#include <core/lethe_grid_tools.h>

#include <fem-dem/gls_vans.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>

// Constructor for class GLS_VANS
template <int dim>
GLSVANSSolver<dim>::GLSVANSSolver(CFDDEMSimulationParameters<dim> &nsparam)
  : GLSNavierStokesSolver<dim>(nsparam.cfd_parameters)
  , cfd_dem_simulation_parameters(nsparam)
  , void_fraction_dof_handler(*this->triangulation)
  , fe_void_fraction(nsparam.cfd_parameters.fem_parameters.void_fraction_order)
  , particle_mapping(1)
  , particle_handler(*this->triangulation,
                     particle_mapping,
                     DEM::get_number_properties())
{
  previous_void_fraction.resize(maximum_number_of_previous_solutions());

  // Check if the simulation has periodic boundaries for fluid.
  std::vector<BoundaryConditions::BoundaryType> boundary_conditions_types =
    this->cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions.type;

  unsigned int n_pbc = 0;
  for (unsigned int i_bc = 0; i_bc < boundary_conditions_types.size(); ++i_bc)
    {
      if (boundary_conditions_types[i_bc] ==
          BoundaryConditions::BoundaryType::periodic)
        {
          if (n_pbc++ > 1)
            {
              throw std::runtime_error(
                "GLS VANS solver does not support more than one periodic boundary condition.");
            }
          else
            {
              has_periodic_boundaries = true;
            }
        }
    }
}

template <int dim>
GLSVANSSolver<dim>::~GLSVANSSolver()
{
  this->dof_handler.clear();
  void_fraction_dof_handler.clear();
}

template <int dim>
void
GLSVANSSolver<dim>::setup_dofs()
{
  GLSNavierStokesSolver<dim>::setup_dofs();

  void_fraction_dof_handler.distribute_dofs(fe_void_fraction);
  locally_owned_dofs_voidfraction =
    void_fraction_dof_handler.locally_owned_dofs();

  DoFTools::extract_locally_relevant_dofs(void_fraction_dof_handler,
                                          locally_relevant_dofs_voidfraction);

  void_fraction_constraints.clear();
  void_fraction_constraints.reinit(locally_relevant_dofs_voidfraction);
  DoFTools::make_hanging_node_constraints(void_fraction_dof_handler,
                                          void_fraction_constraints);

  // Define constraints for periodic boundary conditions
  auto &boundary_conditions =
    this->cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions;
  for (unsigned int i_bc = 0; i_bc < boundary_conditions.size; ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
          BoundaryConditions::BoundaryType::periodic)
        {
          periodic_direction = boundary_conditions.periodic_direction[i_bc];
          DoFTools::make_periodicity_constraints(
            void_fraction_dof_handler,
            boundary_conditions.id[i_bc],
            boundary_conditions.periodic_id[i_bc],
            periodic_direction,
            void_fraction_constraints);

          // Get periodic offset if void fraction method is qcm or spm
          if (this->cfd_dem_simulation_parameters.void_fraction->mode ==
                Parameters::VoidFractionMode::qcm ||
              this->cfd_dem_simulation_parameters.void_fraction->mode ==
                Parameters::VoidFractionMode::spm)
            {
              periodic_offset =
                get_periodic_offset_distance(boundary_conditions.id[i_bc]);
            }
        }
    }

  void_fraction_constraints.close();

  nodal_void_fraction_relevant.reinit(locally_owned_dofs_voidfraction,
                                      locally_relevant_dofs_voidfraction,
                                      this->mpi_communicator);

  for (unsigned int i = 0; i < previous_void_fraction.size(); ++i)
    {
      previous_void_fraction[i].reinit(locally_owned_dofs_voidfraction,
                                       locally_relevant_dofs_voidfraction,
                                       this->mpi_communicator);
    }

  nodal_void_fraction_owned.reinit(locally_owned_dofs_voidfraction,
                                   this->mpi_communicator);

  DynamicSparsityPattern dsp(locally_relevant_dofs_voidfraction);
  DoFTools::make_sparsity_pattern(void_fraction_dof_handler,
                                  dsp,
                                  void_fraction_constraints,
                                  false);
  SparsityTools::distribute_sparsity_pattern(
    dsp,
    locally_owned_dofs_voidfraction,
    this->mpi_communicator,
    locally_relevant_dofs_voidfraction);

  system_matrix_void_fraction.reinit(locally_owned_dofs_voidfraction,
                                     locally_owned_dofs_voidfraction,
                                     dsp,
                                     this->mpi_communicator);

  complete_system_matrix_void_fraction.reinit(locally_owned_dofs_voidfraction,
                                              locally_owned_dofs_voidfraction,
                                              dsp,
                                              this->mpi_communicator);

  system_rhs_void_fraction.reinit(locally_owned_dofs_voidfraction,
                                  this->mpi_communicator);

  complete_system_rhs_void_fraction.reinit(locally_owned_dofs_voidfraction,
                                           this->mpi_communicator);

  active_set.set_size(void_fraction_dof_handler.n_dofs());

  mass_matrix.reinit(locally_owned_dofs_voidfraction,
                     locally_owned_dofs_voidfraction,
                     dsp,
                     this->mpi_communicator);

  assemble_mass_matrix_diagonal(mass_matrix);
}

template <int dim>
void
GLSVANSSolver<dim>::percolate_void_fraction()
{
  for (unsigned int i = previous_void_fraction.size() - 1; i > 0; --i)
    {
      previous_void_fraction[i] = previous_void_fraction[i - 1];
    }
  previous_void_fraction[0] = nodal_void_fraction_relevant;
}

template <int dim>
void
GLSVANSSolver<dim>::finish_time_step_fd()
{
  GLSNavierStokesSolver<dim>::finish_time_step();

  percolate_void_fraction();
}

template <int dim>
void
GLSVANSSolver<dim>::read_dem()
{
  std::string prefix =
    this->cfd_dem_simulation_parameters.void_fraction->dem_file_name;

  // Gather particle serialization information
  std::string   particle_filename = prefix + ".particles";
  std::ifstream input(particle_filename.c_str());
  AssertThrow(input, ExcFileNotOpen(particle_filename));

  std::string buffer;
  std::getline(input, buffer);
  std::istringstream            iss(buffer);
  boost::archive::text_iarchive ia(iss, boost::archive::no_header);

  ia >> particle_handler;

  const std::string filename = prefix + ".triangulation";
  std::ifstream     in(filename.c_str());
  if (!in)
    AssertThrow(false,
                ExcMessage(
                  std::string(
                    "You are trying to restart a previous computation, "
                    "but the restart file <") +
                  filename + "> does not appear to exist!"));

  if (auto parallel_triangulation =
        dynamic_cast<parallel::distributed::Triangulation<dim> *>(
          &*this->triangulation))
    {
      try
        {
          parallel_triangulation->load(filename.c_str());
        }
      catch (...)
        {
          AssertThrow(false,
                      ExcMessage("Cannot open snapshot mesh file or read the"
                                 "triangulation stored there."));
        }
    }
  else
    {
      throw std::runtime_error(
        "VANS equations currently do not support "
        "triangulations other than parallel::distributed");
    }

  // Deserialize particles have the triangulation has been read
  particle_handler.deserialize();
}

template <int dim>
void
GLSVANSSolver<dim>::initialize_void_fraction()
{
  calculate_void_fraction(this->simulation_control->get_current_time(), false);
  for (unsigned int i = 0; i < previous_void_fraction.size(); ++i)
    previous_void_fraction[i] = nodal_void_fraction_relevant;
}

template <int dim>
void
GLSVANSSolver<dim>::calculate_void_fraction(const double time,
                                            bool         load_balance_step)
{
  announce_string(this->pcout, "Void Fraction");

  TimerOutput::Scope t(this->computing_timer, "calculate_void_fraction");
  if (this->cfd_dem_simulation_parameters.void_fraction->mode ==
      Parameters::VoidFractionMode::function)
    {
      const MappingQ<dim> mapping(1);

      this->cfd_dem_simulation_parameters.void_fraction->void_fraction.set_time(
        time);

      VectorTools::interpolate(
        mapping,
        void_fraction_dof_handler,
        this->cfd_dem_simulation_parameters.void_fraction->void_fraction,
        nodal_void_fraction_owned);

      nodal_void_fraction_relevant = nodal_void_fraction_owned;
    }
  else
    {
      if (this->cfd_dem_simulation_parameters.void_fraction->mode ==
          Parameters::VoidFractionMode::pcm)
        particle_centered_method();
      else if (this->cfd_dem_simulation_parameters.void_fraction->mode ==
               Parameters::VoidFractionMode::qcm)
        quadrature_centered_sphere_method(load_balance_step);
      else if (this->cfd_dem_simulation_parameters.void_fraction->mode ==
               Parameters::VoidFractionMode::spm)
        satellite_point_method();

      solve_L2_system_void_fraction();
      if (this->cfd_dem_simulation_parameters.void_fraction
            ->bound_void_fraction == true)
        update_solution_and_constraints();
    }
}

template <int dim>
void
GLSVANSSolver<dim>::assemble_mass_matrix_diagonal(
  TrilinosWrappers::SparseMatrix &mass_matrix)
{
  Assert(fe_void_fraction.degree == 1, ExcNotImplemented());
  QGauss<dim>        quadrature_formula(this->number_quadrature_points);
  FEValues<dim>      fe_void_fraction_values(fe_void_fraction,
                                        quadrature_formula,
                                        update_values | update_JxW_values);
  const unsigned int dofs_per_cell = fe_void_fraction.dofs_per_cell;
  const unsigned int n_qpoints     = quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : void_fraction_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_void_fraction_values.reinit(cell);
          cell_matrix = 0;
          for (unsigned int q = 0; q < n_qpoints; ++q)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_matrix(i, i) += (fe_void_fraction_values.shape_value(i, q) *
                                    fe_void_fraction_values.shape_value(i, q) *
                                    fe_void_fraction_values.JxW(q));
          cell->get_dof_indices(local_dof_indices);
          void_fraction_constraints.distribute_local_to_global(
            cell_matrix, local_dof_indices, mass_matrix);
        }
    }
}

template <int dim>
void
GLSVANSSolver<dim>::update_solution_and_constraints()
{
  const double penalty_parameter = 100;

  TrilinosWrappers::MPI::Vector lambda(locally_owned_dofs_voidfraction);

  nodal_void_fraction_owned = nodal_void_fraction_relevant;

  complete_system_matrix_void_fraction.residual(lambda,
                                                nodal_void_fraction_owned,
                                                system_rhs_void_fraction);

  void_fraction_constraints.clear();

  // reinitialize affine constraints with the correct size
  void_fraction_constraints.reinit(locally_relevant_dofs_voidfraction);

  // Remake hanging node constraints
  DoFTools::make_hanging_node_constraints(void_fraction_dof_handler,
                                          void_fraction_constraints);
  active_set.clear();

  for (const auto &cell : void_fraction_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
               ++v)
            {
              Assert(void_fraction_dof_handler.get_fe().dofs_per_cell ==
                       GeometryInfo<dim>::vertices_per_cell,
                     ExcNotImplemented());
              const unsigned int dof_index = cell->vertex_dof_index(v, 0);
              if (locally_owned_dofs_voidfraction.is_element(dof_index))
                {
                  const double solution_value =
                    nodal_void_fraction_owned(dof_index);
                  if (lambda(dof_index) +
                        penalty_parameter * mass_matrix(dof_index, dof_index) *
                          (solution_value - this->cfd_dem_simulation_parameters
                                              .void_fraction->l2_upper_bound) >
                      0)
                    {
                      active_set.add_index(dof_index);
                      void_fraction_constraints.add_line(dof_index);
                      void_fraction_constraints.set_inhomogeneity(
                        dof_index,
                        this->cfd_dem_simulation_parameters.void_fraction
                          ->l2_upper_bound);
                      nodal_void_fraction_owned(dof_index) =
                        this->cfd_dem_simulation_parameters.void_fraction
                          ->l2_upper_bound;
                      lambda(dof_index) = 0;
                    }
                  else if (lambda(dof_index) +
                             penalty_parameter *
                               mass_matrix(dof_index, dof_index) *
                               (solution_value -
                                this->cfd_dem_simulation_parameters
                                  .void_fraction->l2_lower_bound) <
                           0)
                    {
                      active_set.add_index(dof_index);
                      void_fraction_constraints.add_line(dof_index);
                      void_fraction_constraints.set_inhomogeneity(
                        dof_index,
                        this->cfd_dem_simulation_parameters.void_fraction
                          ->l2_lower_bound);
                      nodal_void_fraction_owned(dof_index) =
                        this->cfd_dem_simulation_parameters.void_fraction
                          ->l2_lower_bound;
                      lambda(dof_index) = 0;
                    }
                }
            }
        }
    }
  active_set.compress();
  nodal_void_fraction_relevant = nodal_void_fraction_owned;
  void_fraction_constraints.close();
}

template <int dim>
void
GLSVANSSolver<dim>::vertices_cell_mapping()
{
  // Find all the cells around each vertex
  TimerOutput::Scope t(this->computing_timer, "vertices_to_cell_map");

  LetheGridTools::vertices_cell_mapping(this->void_fraction_dof_handler,
                                        vertices_to_cell);

  if (has_periodic_boundaries)
    LetheGridTools::vertices_cell_mapping_with_periodic_boundaries(
      this->void_fraction_dof_handler, vertices_to_periodic_cell);
}

template <int dim>
void
GLSVANSSolver<dim>::particle_centered_method()
{
  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(1);

  FEValues<dim> fe_values_void_fraction(mapping,
                                        this->fe_void_fraction,
                                        quadrature_formula,
                                        update_values |
                                          update_quadrature_points |
                                          update_JxW_values | update_gradients);

  const unsigned int dofs_per_cell = this->fe_void_fraction.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  FullMatrix<double> local_matrix_void_fraction(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs_void_fraction(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_vf(dofs_per_cell);
  std::vector<Tensor<1, dim>>          grad_phi_vf(dofs_per_cell);

  system_rhs_void_fraction    = 0;
  system_matrix_void_fraction = 0;

  for (const auto &cell :
       this->void_fraction_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_void_fraction.reinit(cell);

          local_matrix_void_fraction = 0;
          local_rhs_void_fraction    = 0;

          double solid_volume_in_cell = 0;

          // Loop over particles in cell
          // Begin and end iterator for particles in cell
          const auto pic = particle_handler.particles_in_cell(cell);
          for (auto &particle : pic)
            {
              auto particle_properties = particle.get_properties();
              solid_volume_in_cell +=
                M_PI *
                Utilities::fixed_power<dim>(
                  particle_properties[DEM::PropertiesIndex::dp]) /
                (2.0 * dim);
            }
          double cell_volume = cell->measure();

          // Calculate cell void fraction
          double cell_void_fraction =
            (cell_volume - solid_volume_in_cell) / cell_volume;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_vf[k]      = fe_values_void_fraction.shape_value(k, q);
                  grad_phi_vf[k] = fe_values_void_fraction.shape_grad(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Assemble L2 projection
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_void_fraction(i, j) +=
                        (phi_vf[j] * phi_vf[i]) *
                          fe_values_void_fraction.JxW(q) +
                        (this->cfd_dem_simulation_parameters.void_fraction
                           ->l2_smoothing_factor *
                         grad_phi_vf[j] * grad_phi_vf[i] *
                         fe_values_void_fraction.JxW(q));
                    }
                  local_rhs_void_fraction(i) += phi_vf[i] * cell_void_fraction *
                                                fe_values_void_fraction.JxW(q);
                }
            }
          cell->get_dof_indices(local_dof_indices);
          void_fraction_constraints.distribute_local_to_global(
            local_matrix_void_fraction,
            local_rhs_void_fraction,
            local_dof_indices,
            system_matrix_void_fraction,
            system_rhs_void_fraction);
        }
    }
  system_matrix_void_fraction.compress(VectorOperation::add);
  system_rhs_void_fraction.compress(VectorOperation::add);
}

template <int dim>
void
GLSVANSSolver<dim>::quadrature_centered_sphere_method(bool load_balance_step)
{
  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(1);

  FEValues<dim> fe_values_void_fraction(mapping,
                                        this->fe_void_fraction,
                                        quadrature_formula,
                                        update_values |
                                          update_quadrature_points |
                                          update_JxW_values | update_gradients);

  const unsigned int dofs_per_cell = this->fe_void_fraction.dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const unsigned int                   n_q_points = quadrature_formula.size();
  FullMatrix<double>  local_matrix_void_fraction(dofs_per_cell, dofs_per_cell);
  Vector<double>      local_rhs_void_fraction(dofs_per_cell);
  std::vector<double> phi_vf(dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi_vf(dofs_per_cell);

  double r_sphere = 0.0;
  double particles_volume_in_sphere;
  double quadrature_void_fraction;
  double qcm_sphere_diameter =
    this->cfd_dem_simulation_parameters.void_fraction->qcm_sphere_diameter;

  // If the reference sphere diameter is user defined, the radius is calculated
  // from it, otherwise, the value must be calculated while looping over the
  // cells.
  bool calculate_reference_sphere_radius = false;
  if (qcm_sphere_diameter > 1e-16)
    r_sphere = qcm_sphere_diameter / 2.0;
  else
    calculate_reference_sphere_radius = true;

  // Lambda functions for calculating the radius of the reference sphere
  // Calculate the radius by the volume (area in 2D) of sphere:
  // r = (2*dim*V/pi)^(1/dim) / 2
  auto radius_sphere_volume_cell = [](auto cell_measure) {
    return 0.5 * pow(2.0 * dim * cell_measure / M_PI, 1.0 / dim);
  };

  // Calculate the radius is obtained from the volume of sphere based on
  // R_s = h_omega:
  // V_s = pi*(2*V_c^(1/dim))^(dim)/(2*dim) = pi*2^(dim)*V_c/(2*dim)
  auto radius_h_omega = [&radius_sphere_volume_cell](double cell_measure) {
    double reference_sphere_volume =
      M_PI * Utilities::fixed_power<dim>(2.0) * cell_measure / (2.0 * dim);

    return radius_sphere_volume_cell(reference_sphere_volume);
  };

  system_rhs_void_fraction    = 0;
  system_matrix_void_fraction = 0;

  // Clear all contributions of particles from the previous timestep
  for (const auto &cell :
       this->void_fraction_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          const auto pic = particle_handler.particles_in_cell(cell);

          for (auto &particle : pic)
            {
              auto particle_properties = particle.get_properties();

              particle_properties
                [DEM::PropertiesIndex::volumetric_contribution] = 0;
            }
        }
    }

  // Determine the new volumetric contributions of the particles necessary
  // for void fraction calculation
  for (const auto &cell :
       this->void_fraction_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_void_fraction.reinit(cell);

          // Active neighbors include the current cell as well
          auto active_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(vertices_to_cell, cell);

          auto active_periodic_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(
              vertices_to_periodic_cell, cell);

          // Array of real locations for the quadrature points
          std::vector<std::vector<Point<dim>>>
            neighbor_quadrature_point_location(
              active_neighbors.size(), std::vector<Point<dim>>(n_q_points));

          for (unsigned int n = 0; n < active_neighbors.size(); n++)
            {
              fe_values_void_fraction.reinit(active_neighbors[n]);

              neighbor_quadrature_point_location[n] =
                fe_values_void_fraction.get_quadrature_points();
            }

          // Array of real locations for the periodic neighbor quadrature points
          std::vector<std::vector<Point<dim>>>
            periodic_neighbor_quadrature_point_location(
              active_periodic_neighbors.size(),
              std::vector<Point<dim>>(n_q_points));

          for (unsigned int n = 0; n < active_periodic_neighbors.size(); n++)
            {
              fe_values_void_fraction.reinit(active_periodic_neighbors[n]);

              periodic_neighbor_quadrature_point_location[n] =
                fe_values_void_fraction.get_quadrature_points();
            }

          // Loop over the particles in the current cell
          const auto pic = particle_handler.particles_in_cell(cell);
          for (auto &particle : pic)
            {
              auto         particle_properties = particle.get_properties();
              const double r_particle =
                0.5 * particle_properties[DEM::PropertiesIndex::dp];

              // Loop over neighboring cells to determine if a given neighboring
              // particle contributes to the solid volume of the current
              // reference sphere
              //***********************************************************************
              for (unsigned int n = 0; n < active_neighbors.size(); n++)
                {
                  // Define the radius of the reference sphere to be used as the
                  // averaging volume for the QCM, if the reference sphere
                  // diameter was given by the user the value is already defined
                  // since it is not dependent on any measure of the active cell
                  if (calculate_reference_sphere_radius)
                    {
                      if (this->cfd_dem_simulation_parameters.void_fraction
                            ->qcm_sphere_equal_cell_volume == true)
                        {
                          // Get the radius by the volume of sphere which is
                          // equal to the volume of cell
                          r_sphere = radius_sphere_volume_cell(
                            active_neighbors[n]->measure());
                        }
                      else
                        {
                          // The radius is obtained from the volume of sphere
                          // based on R_s = h_omega
                          r_sphere =
                            radius_h_omega(active_neighbors[n]->measure());
                        }
                    }

                  // Loop over quadrature points
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    {
                      // Distance between particle and quadrature point
                      double neighbor_distance =
                        particle.get_location().distance(
                          neighbor_quadrature_point_location[n][k]);

                      // Particle completely in reference sphere
                      if (neighbor_distance <= (r_sphere - r_particle))
                        {
                          particle_properties
                            [DEM::PropertiesIndex::volumetric_contribution] +=
                            M_PI *
                            Utilities::fixed_power<dim>(
                              particle_properties[DEM::PropertiesIndex::dp]) /
                            (2.0 * dim);
                        }

                      // Particle completely outside reference
                      // sphere. Do absolutely nothing.

                      // Particle partially in reference sphere
                      else if ((neighbor_distance > (r_sphere - r_particle)) &&
                               (neighbor_distance < (r_sphere + r_particle)))
                        {
                          if constexpr (dim == 2)
                            particle_properties
                              [DEM::PropertiesIndex::volumetric_contribution] +=
                              particle_circle_intersection_2d(
                                r_particle, r_sphere, neighbor_distance);

                          else if constexpr (dim == 3)
                            particle_properties
                              [DEM::PropertiesIndex::volumetric_contribution] +=
                              particle_sphere_intersection_3d(
                                r_particle, r_sphere, neighbor_distance);
                        }
                    }
                }

              // Loop over periodic neighboring cells to determine if a given
              // neighboring particle contributes to the solid volume of the
              // current reference sphere
              //***********************************************************************
              for (unsigned int n = 0; n < active_periodic_neighbors.size();
                   n++)
                {
                  if (calculate_reference_sphere_radius)
                    {
                      if (this->cfd_dem_simulation_parameters.void_fraction
                            ->qcm_sphere_equal_cell_volume == true)
                        {
                          // Get the radius by the volume of sphere which is
                          // equal to the volume of cell
                          r_sphere = radius_sphere_volume_cell(
                            active_periodic_neighbors[n]->measure());
                        }
                      else
                        {
                          // The radius is obtained from the volume of sphere
                          // based on R_s = h_omega
                          r_sphere = radius_h_omega(
                            active_periodic_neighbors[n]->measure());
                        }
                    }

                  // Loop over quadrature points
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    {
                      // Adjust the location of the particle in the cell to
                      // account for the periodicity. If the position of the
                      // periodic cell if greater than the position of the
                      // current cell, the particle location needs a positive
                      // correction, and vice versa
                      const Point<dim> particle_location =
                        (active_periodic_neighbors[n]
                           ->center()[periodic_direction] >
                         cell->center()[periodic_direction]) ?
                          particle.get_location() + periodic_offset :
                          particle.get_location() - periodic_offset;

                      // Distance between particle and quadrature point
                      double periodic_neighbor_distance =
                        particle_location.distance(
                          periodic_neighbor_quadrature_point_location[n][k]);

                      // Particle completely in reference sphere
                      if (periodic_neighbor_distance <= (r_sphere - r_particle))
                        {
                          particle_properties
                            [DEM::PropertiesIndex::volumetric_contribution] +=
                            M_PI *
                            Utilities::fixed_power<dim>(
                              particle_properties[DEM::PropertiesIndex::dp]) /
                            (2.0 * dim);
                        }

                      // Particle completely outside reference
                      // sphere. Do absolutely nothing.

                      // Particle partially in reference sphere
                      else if ((periodic_neighbor_distance >
                                (r_sphere - r_particle)) &&
                               (periodic_neighbor_distance <
                                (r_sphere + r_particle)))
                        {
                          if constexpr (dim == 2)
                            particle_properties
                              [DEM::PropertiesIndex::volumetric_contribution] +=
                              particle_circle_intersection_2d(
                                r_particle,
                                r_sphere,
                                periodic_neighbor_distance);

                          else if constexpr (dim == 3)
                            particle_properties
                              [DEM::PropertiesIndex::volumetric_contribution] +=
                              particle_sphere_intersection_3d(
                                r_particle,
                                r_sphere,
                                periodic_neighbor_distance);
                        }
                    }
                }

              //*********************************************************************
            }
        }
    }

  // Update ghost particles
  if (load_balance_step)
    {
      particle_handler.sort_particles_into_subdomains_and_cells();
      particle_handler.exchange_ghost_particles(true);
    }
  else
    {
      particle_handler.update_ghost_particles();
    }

  // After the particles' contributions have been determined, calculate and
  // normalize the void fraction
  for (const auto &cell :
       this->void_fraction_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_void_fraction.reinit(cell);

          local_matrix_void_fraction = 0;
          local_rhs_void_fraction    = 0;

          double sum_quadrature_weights = 0;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              sum_quadrature_weights += fe_values_void_fraction.JxW(q);
            }

          // Define the volume of the reference sphere to be used as the
          // averaging volume for the QCM
          if (calculate_reference_sphere_radius)
            {
              if (this->cfd_dem_simulation_parameters.void_fraction
                    ->qcm_sphere_equal_cell_volume == true)
                {
                  // Get the radius by the volume of sphere which is
                  // equal to the volume of cell
                  r_sphere = radius_sphere_volume_cell(cell->measure());
                }
              else
                {
                  // The radius is obtained from the volume of sphere based
                  // on R_s = h_omega
                  r_sphere = radius_h_omega(cell->measure());
                }
            }

          // Array of real locations for the quadrature points
          std::vector<Point<dim>> quadrature_point_location;

          quadrature_point_location =
            fe_values_void_fraction.get_quadrature_points();

          // Active neighbors include the current cell as well
          auto active_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(vertices_to_cell, cell);

          // Periodic neighbors of the current cell
          auto active_periodic_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(
              vertices_to_periodic_cell, cell);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              particles_volume_in_sphere = 0;
              quadrature_void_fraction   = 0;

              for (unsigned int m = 0; m < active_neighbors.size(); m++)
                {
                  // Loop over particles in neighbor cell
                  // Begin and end iterator for particles in neighbor cell
                  const auto pic =
                    particle_handler.particles_in_cell(active_neighbors[m]);
                  for (auto &particle : pic)
                    {
                      double distance            = 0;
                      auto   particle_properties = particle.get_properties();
                      const double r_particle =
                        particle_properties[DEM::PropertiesIndex::dp] * 0.5;
                      double single_particle_volume =
                        M_PI * Utilities::fixed_power<dim>(r_particle * 2.0) /
                        (2 * dim);

                      // Distance between particle and quadrature point
                      // centers
                      distance = particle.get_location().distance(
                        quadrature_point_location[q]);

                      // Particle completely in reference sphere
                      if (distance <= (r_sphere - r_particle))
                        particles_volume_in_sphere +=
                          (M_PI *
                           Utilities::fixed_power<dim>(
                             particle_properties[DEM::PropertiesIndex::dp]) /
                           (2.0 * dim)) *
                          single_particle_volume /
                          particle_properties
                            [DEM::PropertiesIndex::volumetric_contribution];

                      // Particle completely outside reference sphere. Do
                      // absolutely nothing.

                      // Particle partially in reference sphere
                      else if ((distance > (r_sphere - r_particle)) &&
                               (distance < (r_sphere + r_particle)))
                        {
                          if (dim == 2)
                            particles_volume_in_sphere +=
                              particle_circle_intersection_2d(r_particle,
                                                              r_sphere,
                                                              distance) *
                              single_particle_volume /
                              particle_properties
                                [DEM::PropertiesIndex::volumetric_contribution];
                          else if (dim == 3)
                            particles_volume_in_sphere +=
                              particle_sphere_intersection_3d(r_particle,
                                                              r_sphere,
                                                              distance) *
                              single_particle_volume /
                              particle_properties
                                [DEM::PropertiesIndex::volumetric_contribution];
                        }
                    }
                }

              // Execute same operations for periodic neighbors, if the
              // simulation has no periodic boundaries, the container is empty.
              // Also, those operations can not be done in the previous loop
              // because the particles on the periodic side needs a correction
              // with an offset for the distance with the quadrature point
              for (unsigned int m = 0; m < active_periodic_neighbors.size();
                   m++)
                {
                  // Loop over particles in periodic neighbor cell
                  const auto pic = particle_handler.particles_in_cell(
                    active_periodic_neighbors[m]);
                  for (auto &particle : pic)
                    {
                      double distance            = 0;
                      auto   particle_properties = particle.get_properties();
                      const double r_particle =
                        particle_properties[DEM::PropertiesIndex::dp] * 0.5;
                      double single_particle_volume =
                        M_PI * Utilities::fixed_power<dim>(r_particle * 2) /
                        (2 * dim);

                      // Adjust the location of the particle in the cell to
                      // account for the periodicity. If the position of the
                      // periodic cell if greater than the position of the
                      // current cell, the particle location needs a negative
                      // correction, and vice versa. Since the particle is in
                      // the periodic cell, this correction is the inverse of
                      // the correction for the volumetric contribution
                      const Point<dim> particle_location =
                        (active_periodic_neighbors[m]
                           ->center()[periodic_direction] >
                         cell->center()[periodic_direction]) ?
                          particle.get_location() - periodic_offset :
                          particle.get_location() + periodic_offset;

                      // Distance between particle and quadrature point
                      // centers
                      distance = particle_location.distance(
                        quadrature_point_location[q]);

                      // Particle completely in reference sphere
                      if (distance <= (r_sphere - r_particle))
                        particles_volume_in_sphere +=
                          (M_PI *
                           Utilities::fixed_power<dim>(
                             particle_properties[DEM::PropertiesIndex::dp]) /
                           (2.0 * dim)) *
                          single_particle_volume /
                          particle_properties
                            [DEM::PropertiesIndex::volumetric_contribution];

                      // Particle completely outside reference sphere. Do
                      // absolutely nothing.

                      // Particle partially in reference sphere
                      else if ((distance > (r_sphere - r_particle)) &&
                               (distance < (r_sphere + r_particle)))
                        {
                          if (dim == 2)
                            particles_volume_in_sphere +=
                              particle_circle_intersection_2d(r_particle,
                                                              r_sphere,
                                                              distance) *
                              single_particle_volume /
                              particle_properties
                                [DEM::PropertiesIndex::volumetric_contribution];
                          else if (dim == 3)
                            particles_volume_in_sphere +=
                              particle_sphere_intersection_3d(r_particle,
                                                              r_sphere,
                                                              distance) *
                              single_particle_volume /
                              particle_properties
                                [DEM::PropertiesIndex::volumetric_contribution];
                        }
                    }
                }

              // We use the volume of the cell as it is equal to the volume
              // of the sphere
              quadrature_void_fraction =
                ((fe_values_void_fraction.JxW(q) * cell->measure() /
                  sum_quadrature_weights) -
                 particles_volume_in_sphere) /
                (fe_values_void_fraction.JxW(q) * cell->measure() /
                 sum_quadrature_weights);

              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_vf[k]      = fe_values_void_fraction.shape_value(k, q);
                  grad_phi_vf[k] = fe_values_void_fraction.shape_grad(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Assemble L2 projection
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_void_fraction(i, j) +=
                        ((phi_vf[j] * phi_vf[i]) +
                         (this->cfd_dem_simulation_parameters.void_fraction
                            ->l2_smoothing_factor *
                          grad_phi_vf[j] * grad_phi_vf[i])) *
                        fe_values_void_fraction.JxW(q);
                    }

                  local_rhs_void_fraction(i) += phi_vf[i] *
                                                quadrature_void_fraction *
                                                fe_values_void_fraction.JxW(q);
                }
            }

          cell->get_dof_indices(local_dof_indices);
          void_fraction_constraints.distribute_local_to_global(
            local_matrix_void_fraction,
            local_rhs_void_fraction,
            local_dof_indices,
            system_matrix_void_fraction,
            system_rhs_void_fraction);
        }
    }

  system_matrix_void_fraction.compress(VectorOperation::add);
  system_rhs_void_fraction.compress(VectorOperation::add);
}

template <int dim>
void
GLSVANSSolver<dim>::satellite_point_method()
{
  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(1);

  FEValues<dim> fe_values_void_fraction(mapping,
                                        this->fe_void_fraction,
                                        quadrature_formula,
                                        update_values |
                                          update_quadrature_points |
                                          update_JxW_values | update_gradients);

  const unsigned int dofs_per_cell = this->fe_void_fraction.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  FullMatrix<double> local_matrix_void_fraction(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs_void_fraction(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_vf(dofs_per_cell);
  std::vector<Tensor<1, dim>>          grad_phi_vf(dofs_per_cell);

  system_rhs_void_fraction    = 0;
  system_matrix_void_fraction = 0;

  // Creation of reference sphere and components required for mapping into
  // individual particles
  //-------------------------------------------------------------------------
  QGauss<dim>        quadrature_particle(1);
  Triangulation<dim> particle_triangulation;
  Point<dim>         center;

  // Reference particle with radius 1
  GridGenerator::hyper_ball(particle_triangulation, center, 1);
  particle_triangulation.refine_global(
    this->cfd_dem_simulation_parameters.void_fraction
      ->particle_refinement_factor);

  DoFHandler<dim> dof_handler_particle(particle_triangulation);

  FEValues<dim> fe_values_particle(mapping,
                                   this->fe_void_fraction,
                                   quadrature_particle,
                                   update_JxW_values |
                                     update_quadrature_points);

  dof_handler_particle.distribute_dofs(this->fe_void_fraction);

  std::vector<Point<dim>> reference_quadrature_location(
    quadrature_particle.size() *
    particle_triangulation.n_global_active_cells());

  std::vector<double> reference_quadrature_weights(
    quadrature_particle.size() *
    particle_triangulation.n_global_active_cells());

  std::vector<Point<dim>> quadrature_particle_location(
    quadrature_particle.size() *
    particle_triangulation.n_global_active_cells());

  std::vector<double> quadrature_particle_weights(
    quadrature_particle.size() *
    particle_triangulation.n_global_active_cells());

  unsigned int n = 0;
  for (const auto &particle_cell : dof_handler_particle.active_cell_iterators())
    {
      fe_values_particle.reinit(particle_cell);
      for (unsigned int q = 0; q < quadrature_particle.size(); ++q)
        {
          reference_quadrature_weights[n] =
            (M_PI * Utilities::fixed_power<dim>(2.0) / (2.0 * dim)) /
            reference_quadrature_weights.size();
          reference_quadrature_location[n] =
            fe_values_particle.quadrature_point(q);
          n++;
        }
    }
  //-------------------------------------------------------------------------
  for (const auto &cell :
       this->void_fraction_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_void_fraction.reinit(cell);

          local_matrix_void_fraction = 0;
          local_rhs_void_fraction    = 0;

          double solid_volume_in_cell = 0;

          // Active neighbors include the current cell as well
          auto active_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(vertices_to_cell, cell);

          for (unsigned int m = 0; m < active_neighbors.size(); m++)
            {
              // Loop over particles in cell
              // Begin and end iterator for particles in cell
              const auto pic =
                particle_handler.particles_in_cell(active_neighbors[m]);
              for (auto &particle : pic)
                {
                  /*****************************************************************/
                  auto particle_properties = particle.get_properties();
                  auto particle_location   = particle.get_location();

                  // Translation factor used to translate the reference sphere
                  // location and size to those of the particles. Usually, we
                  // take it as the radius of every individual particle. This
                  // makes our method valid for different particle distribution.
                  double translational_factor =
                    particle_properties[DEM::PropertiesIndex::dp] * 0.5;

                  // Resize and translate reference sphere
                  // to the particle size and position according the volume
                  // ratio between sphere and particle.
                  for (unsigned int l = 0;
                       l < reference_quadrature_location.size();
                       ++l)
                    {
                      // For example, in 3D V_particle/V_sphere =
                      // r_particle³/r_sphere³ and since r_sphere is always 1,
                      // then V_particle/V_sphere = r_particle³. Therefore,
                      // V_particle = r_particle³ * V_sphere.
                      quadrature_particle_weights[l] =
                        Utilities::fixed_power<dim>(translational_factor) *
                        reference_quadrature_weights[l];

                      // Here, we translate the position of the reference sphere
                      // into the position of the particle, but for this we have
                      // to shrink or expand the size of the reference sphere to
                      // be equal to the size of the particle as the location of
                      // the quadrature points is affected by the size by
                      // multiplying with the particle's radius. We then
                      // translate by taking the translational vector between
                      // the reference sphere center and the particle's center.
                      // This translates directly into the translational vector
                      // being the particle's position as the reference sphere
                      // is always located at (0,0) in 2D or (0,0,0) in 3D.
                      quadrature_particle_location[l] =
                        (translational_factor *
                         reference_quadrature_location[l]) +
                        particle_location;

                      if (cell->point_inside(quadrature_particle_location[l]))
                        solid_volume_in_cell += quadrature_particle_weights[l];
                    }
                }
            }

          // Same steps for the periodic neighbors with particle location
          // correction
          auto active_periodic_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(
              vertices_to_periodic_cell, cell);

          for (unsigned int m = 0; m < active_periodic_neighbors.size(); m++)
            {
              const auto pic = particle_handler.particles_in_cell(
                active_periodic_neighbors[m]);
              for (auto &particle : pic)
                {
                  /*****************************************************************/
                  auto particle_properties = particle.get_properties();
                  const Point<dim> particle_location =
                    (active_periodic_neighbors[m]
                       ->center()[periodic_direction] >
                     cell->center()[periodic_direction]) ?
                      particle.get_location() - periodic_offset :
                      particle.get_location() + periodic_offset;

                  double translational_factor =
                    particle_properties[DEM::PropertiesIndex::dp] * 0.5;

                  for (unsigned int l = 0;
                       l < reference_quadrature_location.size();
                       ++l)
                    {
                      quadrature_particle_weights[l] =
                        Utilities::fixed_power<dim>(translational_factor) *
                        reference_quadrature_weights[l];

                      quadrature_particle_location[l] =
                        (translational_factor *
                         reference_quadrature_location[l]) +
                        particle_location;

                      if (cell->point_inside(quadrature_particle_location[l]))
                        solid_volume_in_cell += quadrature_particle_weights[l];
                    }
                }
            }

          double cell_volume = cell->measure();

          // Calculate cell void fraction
          double cell_void_fraction =
            (cell_volume - solid_volume_in_cell) / cell_volume;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_vf[k]      = fe_values_void_fraction.shape_value(k, q);
                  grad_phi_vf[k] = fe_values_void_fraction.shape_grad(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_void_fraction(i, j) +=
                        (phi_vf[j] * phi_vf[i]) *
                          fe_values_void_fraction.JxW(q) +
                        (this->cfd_dem_simulation_parameters.void_fraction
                           ->l2_smoothing_factor *
                         grad_phi_vf[j] * grad_phi_vf[i] *
                         fe_values_void_fraction.JxW(q));
                    }
                  local_rhs_void_fraction(i) += phi_vf[i] * cell_void_fraction *
                                                fe_values_void_fraction.JxW(q);
                }
            }
          cell->get_dof_indices(local_dof_indices);
          void_fraction_constraints.distribute_local_to_global(
            local_matrix_void_fraction,
            local_rhs_void_fraction,
            local_dof_indices,
            system_matrix_void_fraction,
            system_rhs_void_fraction);
        }
    }

  system_matrix_void_fraction.compress(VectorOperation::add);
  system_rhs_void_fraction.compress(VectorOperation::add);
}


template <int dim>
void
GLSVANSSolver<dim>::solve_L2_system_void_fraction()
{
  // Solve the L2 projection system
  const double linear_solver_tolerance = 1e-15;

  if (this->cfd_dem_simulation_parameters.cfd_parameters.linear_solver
        .at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const IndexSet locally_owned_dofs_voidfraction =
    void_fraction_dof_handler.locally_owned_dofs();

  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    locally_owned_dofs_voidfraction, this->mpi_communicator);

  SolverControl solver_control(
    this->cfd_dem_simulation_parameters.cfd_parameters.linear_solver
      .at(PhysicsID::fluid_dynamics)
      .max_iterations,
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverCG solver(solver_control);

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill =
    this->cfd_dem_simulation_parameters.cfd_parameters.linear_solver
      .at(PhysicsID::fluid_dynamics)
      .ilu_precond_fill;
  const double ilu_atol =
    this->cfd_dem_simulation_parameters.cfd_parameters.linear_solver
      .at(PhysicsID::fluid_dynamics)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->cfd_dem_simulation_parameters.cfd_parameters.linear_solver
      .at(PhysicsID::fluid_dynamics)
      .ilu_precond_rtol;

  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_matrix_void_fraction,
                                 preconditionerOptions);

  solver.solve(system_matrix_void_fraction,
               completely_distributed_solution,
               system_rhs_void_fraction,
               *ilu_preconditioner);

  if (this->cfd_dem_simulation_parameters.cfd_parameters.linear_solver
        .at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  void_fraction_constraints.distribute(completely_distributed_solution);
  nodal_void_fraction_relevant = completely_distributed_solution;
}

// Do an iteration with the NavierStokes Solver
// Handles the fact that we may or may not be at a first
// iteration with the solver and sets the initial conditions
template <int dim>
void
GLSVANSSolver<dim>::iterate()
{
  announce_string(this->pcout, "Volume-Averaged Fluid Dynamics");
  this->forcing_function->set_time(
    this->simulation_control->get_current_time());

  this->simulation_control->set_assembly_method(
    this->cfd_dem_simulation_parameters.cfd_parameters.simulation_control
      .method);
  PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(false);
}

template <int dim>
void
GLSVANSSolver<dim>::setup_assemblers()
{
  this->assemblers.clear();
  particle_fluid_assemblers.clear();

  if (this->check_existance_of_bc(
        BoundaryConditions::BoundaryType::function_weak))
    {
      this->assemblers.push_back(
        std::make_shared<WeakDirichletBoundaryCondition<dim>>(
          this->simulation_control,
          this->simulation_parameters.boundary_conditions));
    }
  if (this->check_existance_of_bc(
        BoundaryConditions::BoundaryType::partial_slip))
    {
      this->assemblers.push_back(
        std::make_shared<PartialSlipDirichletBoundaryCondition<dim>>(
          this->simulation_control,
          this->simulation_parameters.boundary_conditions));
    }
  if (this->check_existance_of_bc(BoundaryConditions::BoundaryType::outlet))
    {
      this->assemblers.push_back(std::make_shared<OutletBoundaryCondition<dim>>(
        this->simulation_control,
        this->simulation_parameters.boundary_conditions));
    }
  if (this->check_existance_of_bc(BoundaryConditions::BoundaryType::pressure))
    {
      this->assemblers.push_back(
        std::make_shared<PressureBoundaryCondition<dim>>(
          this->simulation_control,
          this->simulation_parameters.boundary_conditions));
    }

  if (this->cfd_dem_simulation_parameters.cfd_dem.drag_force == true)
    {
      // Particle_Fluid Interactions Assembler
      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::difelice)
        {
          // DiFelice Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<GLSVansAssemblerDiFelice<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }

      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::rong)
        {
          // Rong Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<GLSVansAssemblerRong<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }

      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::dallavalle)
        {
          // Dallavalle Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<GLSVansAssemblerDallavalle<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }

      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::kochhill)
        {
          // Koch and Hill Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<GLSVansAssemblerKochHill<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }
      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::beetstra)
        {
          // Beetstra drag model assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<GLSVansAssemblerBeetstra<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }
      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::gidaspow)
        {
          // Gidaspow Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<GLSVansAssemblerGidaspow<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }
    }

  if (this->cfd_dem_simulation_parameters.cfd_dem.saffman_lift_force == true)
    // Saffman Mei Lift Force Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<GLSVansAssemblerSaffmanMei<dim>>(
        this->cfd_dem_simulation_parameters.dem_parameters
          .lagrangian_physical_properties));

  if (this->cfd_dem_simulation_parameters.cfd_dem.magnus_lift_force == true)
    // Magnus Lift Force Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<GLSVansAssemblerMagnus<dim>>(
        this->cfd_dem_simulation_parameters.dem_parameters
          .lagrangian_physical_properties));

  if (this->cfd_dem_simulation_parameters.cfd_dem.rotational_viscous_torque ==
      true)
    // Viscous Torque Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<GLSVansAssemblerViscousTorque<dim>>(
        this->cfd_dem_simulation_parameters.dem_parameters
          .lagrangian_physical_properties));

  if (this->cfd_dem_simulation_parameters.cfd_dem.vortical_viscous_torque ==
      true)
    // Viscous Torque Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<GLSVansAssemblerVorticalTorque<dim>>(
        this->cfd_dem_simulation_parameters.dem_parameters
          .lagrangian_physical_properties));

  if (this->cfd_dem_simulation_parameters.cfd_dem.buoyancy_force == true)
    // Buoyancy Force Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<GLSVansAssemblerBuoyancy<dim>>(
        this->cfd_dem_simulation_parameters.dem_parameters
          .lagrangian_physical_properties));

  if (this->cfd_dem_simulation_parameters.cfd_dem.pressure_force == true)
    // Pressure Force
    particle_fluid_assemblers.push_back(
      std::make_shared<GLSVansAssemblerPressureForce<dim>>(
        this->cfd_dem_simulation_parameters.cfd_dem));

  if (this->cfd_dem_simulation_parameters.cfd_dem.shear_force == true)
    // Shear Force
    particle_fluid_assemblers.push_back(
      std::make_shared<GLSVansAssemblerShearForce<dim>>(
        this->cfd_dem_simulation_parameters.cfd_dem));

  // Time-stepping schemes
  if (is_bdf(this->simulation_control->get_assembly_method()))
    {
      this->assemblers.push_back(std::make_shared<GLSVansAssemblerBDF<dim>>(
        this->simulation_control, this->cfd_dem_simulation_parameters.cfd_dem));
    }

  //  Fluid_Particle Interactions Assembler
  this->assemblers.push_back(std::make_shared<GLSVansAssemblerFPI<dim>>(
    this->cfd_dem_simulation_parameters.cfd_dem));

  // The core assembler should always be the last assembler to be called
  // in the stabilized formulation as to have all strong residual and
  // jacobian stored. Core assembler
  if (this->cfd_dem_simulation_parameters.cfd_dem.vans_model ==
      Parameters::VANSModel::modelA)
    this->assemblers.push_back(
      std::make_shared<GLSVansAssemblerCoreModelA<dim>>(
        this->simulation_control, this->cfd_dem_simulation_parameters.cfd_dem));

  if (this->cfd_dem_simulation_parameters.cfd_dem.vans_model ==
      Parameters::VANSModel::modelB)
    this->assemblers.push_back(
      std::make_shared<GLSVansAssemblerCoreModelB<dim>>(
        this->simulation_control, this->cfd_dem_simulation_parameters.cfd_dem));
}

template <int dim>
void
GLSVANSSolver<dim>::assemble_system_matrix()
{
  {
    TimerOutput::Scope t(this->computing_timer, "Assemble_Matrix");
    this->system_matrix = 0;

    setup_assemblers();

    auto scratch_data = NavierStokesScratchData<dim>(
      this->simulation_parameters.physical_properties_manager,
      *this->fe,
      *this->cell_quadrature,
      *this->mapping,
      *this->face_quadrature);

    scratch_data.enable_void_fraction(fe_void_fraction,
                                      *this->cell_quadrature,
                                      *this->mapping);

    scratch_data.enable_particle_fluid_interactions(
      particle_handler.n_global_max_particles_per_cell(),
      this->cfd_dem_simulation_parameters.cfd_dem.interpolated_void_fraction);

    WorkStream::run(
      this->dof_handler.begin_active(),
      this->dof_handler.end(),
      *this,
      &GLSVANSSolver::assemble_local_system_matrix,
      &GLSVANSSolver::copy_local_matrix_to_global_matrix,
      scratch_data,
      StabilizedMethodsTensorCopyData<dim>(this->fe->n_dofs_per_cell(),
                                           this->cell_quadrature->size()));
    this->system_matrix.compress(VectorOperation::add);
  }
  this->setup_preconditioner();
}

template <int dim>
void
GLSVANSSolver<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim>                         &scratch_data,
  StabilizedMethodsTensorCopyData<dim>                 &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(
    cell,
    this->evaluation_point,
    this->previous_solutions,
    this->forcing_function,
    this->flow_control.get_beta(),
    this->simulation_parameters.stabilization.pressure_scaling_factor);

  typename DoFHandler<dim>::active_cell_iterator void_fraction_cell(
    &(*(this->triangulation)),
    cell->level(),
    cell->index(),
    &this->void_fraction_dof_handler);

  scratch_data.reinit_void_fraction(void_fraction_cell,
                                    nodal_void_fraction_relevant,
                                    previous_void_fraction);

  scratch_data.reinit_particle_fluid_interactions(cell,
                                                  this->evaluation_point,
                                                  this->previous_solutions[0],
                                                  nodal_void_fraction_relevant,
                                                  particle_handler,
                                                  this->dof_handler,
                                                  void_fraction_dof_handler);
  scratch_data.calculate_physical_properties();
  copy_data.reset();

  for (auto &pf_assembler : particle_fluid_assemblers)
    {
      pf_assembler->calculate_particle_fluid_interactions(scratch_data);
    }

  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_matrix(scratch_data, copy_data);
    }

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
GLSVANSSolver<dim>::copy_local_matrix_to_global_matrix(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_matrix,
                                              copy_data.local_dof_indices,
                                              this->system_matrix);
}

template <int dim>
void
GLSVANSSolver<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble_RHS");

  this->system_rhs = 0;

  setup_assemblers();

  auto scratch_data = NavierStokesScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    *this->face_quadrature);

  scratch_data.enable_void_fraction(fe_void_fraction,
                                    *this->cell_quadrature,
                                    *this->mapping);

  scratch_data.enable_particle_fluid_interactions(
    particle_handler.n_global_max_particles_per_cell(),
    this->cfd_dem_simulation_parameters.cfd_dem.interpolated_void_fraction);

  WorkStream::run(
    this->dof_handler.begin_active(),
    this->dof_handler.end(),
    *this,
    &GLSVANSSolver::assemble_local_system_rhs,
    &GLSVANSSolver::copy_local_rhs_to_global_rhs,
    scratch_data,
    StabilizedMethodsTensorCopyData<dim>(this->fe->n_dofs_per_cell(),
                                         this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);

  if (this->simulation_control->is_first_assembly())
    this->simulation_control->provide_residual(this->system_rhs.l2_norm());
}

template <int dim>
void
GLSVANSSolver<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim>                         &scratch_data,
  StabilizedMethodsTensorCopyData<dim>                 &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(
    cell,
    this->evaluation_point,
    this->previous_solutions,
    this->forcing_function,
    this->flow_control.get_beta(),
    this->simulation_parameters.stabilization.pressure_scaling_factor);

  typename DoFHandler<dim>::active_cell_iterator void_fraction_cell(
    &(*(this->triangulation)),
    cell->level(),
    cell->index(),
    &this->void_fraction_dof_handler);

  scratch_data.reinit_void_fraction(void_fraction_cell,
                                    nodal_void_fraction_relevant,
                                    previous_void_fraction);

  scratch_data.reinit_particle_fluid_interactions(cell,
                                                  this->evaluation_point,
                                                  this->previous_solutions[0],
                                                  nodal_void_fraction_relevant,
                                                  particle_handler,
                                                  this->dof_handler,
                                                  void_fraction_dof_handler);

  scratch_data.calculate_physical_properties();
  copy_data.reset();

  for (auto &pf_assembler : particle_fluid_assemblers)
    {
      pf_assembler->calculate_particle_fluid_interactions(scratch_data);
    }

  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_rhs(scratch_data, copy_data);
    }

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
GLSVANSSolver<dim>::copy_local_rhs_to_global_rhs(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_rhs,
                                              copy_data.local_dof_indices,
                                              this->system_rhs);
}

template <int dim>
void
GLSVANSSolver<dim>::output_field_hook(DataOut<dim> &data_out)
{
  data_out.add_data_vector(void_fraction_dof_handler,
                           nodal_void_fraction_relevant,
                           "void_fraction");
}

template <int dim>
void
GLSVANSSolver<dim>::monitor_mass_conservation()
{
  QGauss<dim> quadrature_formula(this->number_quadrature_points);

  FEValues<dim> fe_values(*this->mapping,
                          *this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients |
                            update_hessians);

  FEValues<dim> fe_values_void_fraction(*this->mapping,
                                        this->fe_void_fraction,
                                        quadrature_formula,
                                        update_values |
                                          update_quadrature_points |
                                          update_JxW_values | update_gradients);

  const FEValuesExtractors::Vector velocities(0);

  const unsigned int n_q_points = quadrature_formula.size();

  std::vector<double>         present_void_fraction_values(n_q_points);
  std::vector<Tensor<1, dim>> present_void_fraction_gradients(n_q_points);
  // Values at previous time step for transient schemes for void
  // fraction
  std::vector<double> p1_void_fraction_values(n_q_points);
  std::vector<double> p2_void_fraction_values(n_q_points);
  std::vector<double> p3_void_fraction_values(n_q_points);

  std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);
  std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);

  // Values of the force function, the last component of which is the external
  // mass source
  std::vector<Vector<double>> rhs_force(n_q_points, Vector<double>(dim + 1));

  double continuity           = 0;
  double max_local_continuity = 0;
  double local_mass_source    = 0;

  Vector<double>      bdf_coefs;
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  const auto scheme = this->simulation_control->get_assembly_method();

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf1)
    bdf_coefs = bdf_coefficients(1, time_steps_vector);

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    bdf_coefs = bdf_coefficients(2, time_steps_vector);

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    bdf_coefs = bdf_coefficients(3, time_steps_vector);

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          typename DoFHandler<dim>::active_cell_iterator void_fraction_cell(
            &(*this->triangulation),
            cell->level(),
            cell->index(),
            &this->void_fraction_dof_handler);
          fe_values_void_fraction.reinit(void_fraction_cell);

          // Gather void fraction (values, gradient)
          fe_values_void_fraction.get_function_values(
            nodal_void_fraction_relevant, present_void_fraction_values);
          fe_values_void_fraction.get_function_gradients(
            nodal_void_fraction_relevant, present_void_fraction_gradients);

          fe_values.reinit(cell);

          // Gather velocity (values and gradient)
          auto &evaluation_point = this->evaluation_point;
          fe_values[velocities].get_function_values(evaluation_point,
                                                    present_velocity_values);
          fe_values[velocities].get_function_gradients(
            evaluation_point, present_velocity_gradients);

          // Gather the external mass source at the quadrature point location
          std::vector<Point<dim>> quadrature_points =
            fe_values.get_quadrature_points();

          this->forcing_function->vector_value_list(quadrature_points,
                                                    rhs_force);

          // Gather the previous time steps depending on the number of
          // stages of the time integration scheme for the void fraction

          if (scheme !=
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            {
              fe_values_void_fraction.get_function_values(
                previous_void_fraction[0], p1_void_fraction_values);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                fe_values_void_fraction.get_function_values(
                  previous_void_fraction[1], p2_void_fraction_values);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                fe_values_void_fraction.get_function_values(
                  previous_void_fraction[2], p3_void_fraction_values);
            }

          local_mass_source = 0;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // Calculate external mass source
              const unsigned int component_mass =
                this->fe->system_to_component_index(dim).first;
              double external_mass_source = rhs_force[q](component_mass);

              // Calculate the divergence of the velocity
              const double present_velocity_divergence =
                trace(present_velocity_gradients[q]);


              // Evaluation of global mass conservation
              local_mass_source =
                (present_velocity_values[q] *
                   present_void_fraction_gradients[q] +
                 present_void_fraction_values[q] * present_velocity_divergence -
                 external_mass_source) *
                fe_values_void_fraction.JxW(q);

              if (scheme ==
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
                  scheme == Parameters::SimulationControl::TimeSteppingMethod::
                              steady_bdf)
                local_mass_source +=
                  (bdf_coefs[0] * present_void_fraction_values[q] +
                   bdf_coefs[1] * p1_void_fraction_values[q]) *
                  fe_values_void_fraction.JxW(q);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                local_mass_source +=
                  (bdf_coefs[0] * present_void_fraction_values[q] +
                   bdf_coefs[1] * p1_void_fraction_values[q] +
                   bdf_coefs[2] * p2_void_fraction_values[q]) *
                  fe_values_void_fraction.JxW(q);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                local_mass_source +=
                  (bdf_coefs[0] * present_void_fraction_values[q] +
                   bdf_coefs[1] * p1_void_fraction_values[q] +
                   bdf_coefs[2] * p2_void_fraction_values[q] +
                   bdf_coefs[3] * p3_void_fraction_values[q]) *
                  fe_values_void_fraction.JxW(q);

              continuity += local_mass_source;
            }

          max_local_continuity =
            std::max(max_local_continuity, abs(local_mass_source));
        }
    }

  continuity = Utilities::MPI::sum(continuity, this->mpi_communicator);
  max_local_continuity =
    Utilities::MPI::max(max_local_continuity, this->mpi_communicator);


  this->pcout << "Global continuity equation error: " << continuity << " s^-1"
              << std::endl;
  this->pcout << "Max local continuity error: " << max_local_continuity
              << " s^-1" << std::endl;
}

template <int dim>
void
GLSVANSSolver<dim>::solve()
{
  // This is enforced to 1 right now because it does not provide
  // better speed-up than using MPI. This could be eventually changed...
  MultithreadInfo::set_thread_limit(1);

  read_mesh_and_manifolds(
    *this->triangulation,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh,
    this->cfd_dem_simulation_parameters.cfd_parameters.manifolds_parameters,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
        .restart ||
      this->cfd_dem_simulation_parameters.void_fraction->read_dem == true,
    this->cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions);

  if (this->cfd_dem_simulation_parameters.void_fraction->read_dem == true &&
      this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
          .restart == false)
    read_dem();

  setup_dofs();

  this->set_initial_condition(
    this->cfd_dem_simulation_parameters.cfd_parameters.initial_condition->type,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
      .restart);

  particle_handler.exchange_ghost_particles(true);

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      if ((this->simulation_control->get_step_number() %
               this->simulation_parameters.mesh_adaptation.frequency !=
             0 ||
           this->simulation_parameters.mesh_adaptation.type ==
             Parameters::MeshAdaptation::Type::none ||
           this->simulation_control->is_at_start()) &&
          this->simulation_parameters.boundary_conditions.time_dependent)
        {
          this->update_boundary_conditions();
        }

      this->dynamic_flow_control();

      if (this->simulation_control->is_at_start())
        {
          vertices_cell_mapping();
          initialize_void_fraction();
          this->iterate();
        }
      else
        {
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          vertices_cell_mapping();
          calculate_void_fraction(this->simulation_control->get_current_time(),
                                  false);
          this->iterate();
        }

      this->postprocess(false);
      monitor_mass_conservation();
      finish_time_step_fd();
    }

  this->finish_simulation();
}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the
// library is valid before we actually compile the solver This greatly
// helps with debugging
template class GLSVANSSolver<2>;
template class GLSVANSSolver<3>;
