#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <fem-dem/gls_vans.h>


template <int dim>
double
calculate_pressure_drop(const DoFHandler<dim> &              dof_handler,
                        std::shared_ptr<Mapping<dim>>        mapping,
                        const MPI_Comm &                     mpi_communicator,
                        std::shared_ptr<FESystem<dim>>       fe,
                        const TrilinosWrappers::MPI::Vector &evaluation_point,
                        const unsigned int number_quadrature_points,
                        double             inlet_boundary_id,
                        double             outlet_boundary_id)
{
  QGauss<dim>     quadrature_formula(number_quadrature_points);
  QGauss<dim - 1> face_quadrature_formula(number_quadrature_points);

  FEValues<dim>     fe_values(*mapping,
                          *fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients |
                            update_hessians);
  FEFaceValues<dim> fe_face_values(*fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors | update_JxW_values);

  const FEValuesExtractors::Scalar pressure(dim);

  const unsigned int n_q_points      = quadrature_formula.size();
  const unsigned int face_n_q_points = face_quadrature_formula.size();

  std::vector<double> present_pressure_values(n_q_points);

  double pressure_upper_boundary = 0;
  double upper_surface           = 0;
  double pressure_lower_boundary = 0;
  double lower_surface           = 0;
  double pressure_drop           = 0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          // Gather pressure (values)
          fe_values[pressure].get_function_values(evaluation_point,
                                                  present_pressure_values);

          for (unsigned int q = 0; q < face_n_q_points; ++q)
            {
              for (const auto &face : cell->face_iterators())
                {
                  if (face->at_boundary() &&
                      (face->boundary_id() == outlet_boundary_id))
                    {
                      fe_face_values.reinit(cell, face);
                      pressure_upper_boundary +=
                        present_pressure_values[q] * fe_face_values.JxW(q);

                      upper_surface += fe_face_values.JxW(q);
                    }

                  if (face->at_boundary() &&
                      (face->boundary_id() == inlet_boundary_id))
                    {
                      fe_face_values.reinit(cell, face);
                      pressure_lower_boundary +=
                        present_pressure_values[q] * fe_face_values.JxW(q);

                      lower_surface += fe_face_values.JxW(q);
                    }
                }
            }
        }
    }

  pressure_lower_boundary =
    Utilities::MPI::sum(pressure_lower_boundary, mpi_communicator);
  lower_surface = Utilities::MPI::sum(lower_surface, mpi_communicator);
  pressure_upper_boundary =
    Utilities::MPI::sum(pressure_upper_boundary, mpi_communicator);
  upper_surface = Utilities::MPI::sum(upper_surface, mpi_communicator);

  pressure_upper_boundary = pressure_upper_boundary / upper_surface;
  pressure_lower_boundary = pressure_lower_boundary / lower_surface;
  pressure_drop           = pressure_lower_boundary - pressure_upper_boundary;

  return pressure_drop;
}

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

#if DEAL_II_VERSION_GTE(10, 0, 0)
  fluid_solid_force.resize(particle_handler.get_max_local_particle_index());
#else
  {
    unsigned int max_particle_id = 0;
    for (const auto &particle : particle_handler)
      max_particle_id = std::max(max_particle_id, particle.get_id());
    fluid_solid_force.resize(max_particle_id + 1);
  }
#endif
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
  GLSNavierStokesSolver<dim>::finish_time_step_fd();

  percolate_void_fraction();
}

template <int dim>
void
GLSVANSSolver<dim>::read_dem()
{
  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  std::string prefix =
    this->cfd_dem_simulation_parameters.void_fraction->dem_file_name;

  parallel_triangulation->signals.post_distributed_load.connect(
    std::bind(&Particles::ParticleHandler<dim>::register_load_callback_function,
              &particle_handler,
              true));

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
        "VANS equations currently do not support triangulations other than parallel::distributed");
    }
}


template <int dim>
void
GLSVANSSolver<dim>::initialize_void_fraction()
{
  calculate_void_fraction(this->simulation_control->get_current_time());
  for (unsigned int i = 0; i < previous_void_fraction.size(); ++i)
    previous_void_fraction[i] = nodal_void_fraction_relevant;
}

template <int dim>
void
GLSVANSSolver<dim>::calculate_void_fraction(const double time)
{
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
  else if (this->cfd_dem_simulation_parameters.void_fraction->mode ==
           Parameters::VoidFractionMode::dem)
    {
      assemble_L2_projection_void_fraction();
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
  active_set.clear();
  std::vector<bool> dof_touched(void_fraction_dof_handler.n_dofs(), false);

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
GLSVANSSolver<dim>::assemble_L2_projection_void_fraction()
{
  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(1,
                              this->cfd_dem_simulation_parameters.cfd_parameters
                                .fem_parameters.qmapping_all);

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

          double particles_volume_in_cell = 0;

          // Loop over particles in cell
          // Begin and end iterator for particles in cell
          const auto pic = particle_handler.particles_in_cell(cell);
          for (auto &particle : pic)
            {
              auto particle_properties = particle.get_properties();
              particles_volume_in_cell +=
                M_PI * pow(particle_properties[DEM::PropertiesIndex::dp], dim) /
                (2 * dim);
            }
          double cell_volume = cell->measure();

          // Calculate cell void fraction
          double cell_void_fraction =
            (cell_volume - particles_volume_in_cell) / cell_volume;

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
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const IndexSet locally_owned_dofs_voidfraction =
    void_fraction_dof_handler.locally_owned_dofs();

  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    locally_owned_dofs_voidfraction, this->mpi_communicator);


  SolverControl solver_control(this->cfd_dem_simulation_parameters
                                 .cfd_parameters.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverCG solver(solver_control);

  TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill = this->cfd_dem_simulation_parameters.cfd_parameters
                            .linear_solver.ilu_precond_fill;
  const double ilu_atol = this->cfd_dem_simulation_parameters.cfd_parameters
                            .linear_solver.ilu_precond_atol;
  const double ilu_rtol = this->cfd_dem_simulation_parameters.cfd_parameters
                            .linear_solver.ilu_precond_rtol;

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
GLSVANSSolver<dim>::first_iteration()
{
  // First step if the method is not a multi-step method
  if (!is_bdf_high_order(this->cfd_dem_simulation_parameters.cfd_parameters
                           .simulation_control.method))
    {
      iterate();
    }

  // Taking care of the multi-step methods
  else if (this->cfd_dem_simulation_parameters.cfd_parameters.simulation_control
             .method == Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    {
      Parameters::SimulationControl timeParameters =
        this->cfd_dem_simulation_parameters.cfd_parameters.simulation_control;

      // Start the BDF2 with a single Euler time step with a lower time step
      double time_step =
        timeParameters.dt * timeParameters.startup_timestep_scaling;
      this->simulation_control->set_current_time_step(time_step);

      double intermediate_time =
        this->simulation_control->get_current_time() + time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);


      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf1, false);
      this->percolate_time_vectors_fd();
      percolate_void_fraction();

      // Reset the time step and do a bdf 2 newton iteration using the two
      // steps to complete the full step

      time_step =
        timeParameters.dt * (1. - timeParameters.startup_timestep_scaling);

      this->simulation_control->set_current_time_step(time_step);
      intermediate_time += time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);


      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf2, false);

      this->simulation_control->set_suggested_time_step(timeParameters.dt);
    }

  else if (this->cfd_dem_simulation_parameters.cfd_parameters.simulation_control
             .method == Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    {
      Parameters::SimulationControl timeParameters =
        this->cfd_dem_simulation_parameters.cfd_parameters.simulation_control;

      // Start the BDF3 with a single Euler time step with a lower time step
      double time_step =
        timeParameters.dt * timeParameters.startup_timestep_scaling;

      this->simulation_control->set_current_time_step(time_step);

      double intermediate_time = time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);


      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf1, false);
      this->percolate_time_vectors_fd();
      percolate_void_fraction();

      // Reset the time step and do a bdf 2 newton iteration using the two
      // steps

      this->simulation_control->set_current_time_step(time_step);
      intermediate_time += time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);


      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf1, false);
      this->percolate_time_vectors_fd();
      percolate_void_fraction();

      // Reset the time step and do a bdf 3 newton iteration using the two
      // steps to complete the full step
      time_step =
        timeParameters.dt * (1. - 2. * timeParameters.startup_timestep_scaling);
      this->simulation_control->set_current_time_step(time_step);
      intermediate_time += time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);


      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf3, false);
      this->simulation_control->set_suggested_time_step(timeParameters.dt);
    }
}

// Do an iteration with the NavierStokes Solver
// Handles the fact that we may or may not be at a first
// iteration with the solver and sets the initial conditions
template <int dim>
void
GLSVANSSolver<dim>::iterate()
{
  calculate_void_fraction(this->simulation_control->get_current_time());
  this->forcing_function->set_time(
    this->simulation_control->get_current_time());
  PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
    this->cfd_dem_simulation_parameters.cfd_parameters.simulation_control
      .method,
    false);
}

template <int dim>
void
GLSVANSSolver<dim>::setup_assemblers()
{
  this->assemblers.clear();
  particle_fluid_assemblers.clear();

  // Particle_Fluid Interactions Assembler
  if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
      Parameters::DragModel::difelice)
    {
      // DiFelice Model drag Assembler
      particle_fluid_assemblers.push_back(
        std::make_shared<GLSVansAssemblerDiFelice<dim>>(
          this->simulation_control,
          this->cfd_dem_simulation_parameters.cfd_parameters
            .physical_properties));
    }

  if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
      Parameters::DragModel::rong)
    {
      // Rong Model drag Assembler
      particle_fluid_assemblers.push_back(
        std::make_shared<GLSVansAssemblerRong<dim>>(
          this->simulation_control,
          this->cfd_dem_simulation_parameters.cfd_parameters
            .physical_properties));
    }

  // Buoyancy Force Assembler
  particle_fluid_assemblers.push_back(
    std::make_shared<GLSVansAssemblerBuoyancy<dim>>(
      this->simulation_control,
      this->cfd_dem_simulation_parameters.cfd_parameters.physical_properties,
      this->cfd_dem_simulation_parameters.dem_parameters
        .lagrangian_physical_properties));


  // Time-stepping schemes
  if (is_bdf(this->simulation_control->get_assembly_method()))
    {
      this->assemblers.push_back(std::make_shared<GLSVansAssemblerBDF<dim>>(
        this->simulation_control, this->cfd_dem_simulation_parameters.cfd_dem));
    }

  //  Fluid_Particle Interactions Assembler
  this->assemblers.push_back(std::make_shared<GLSVansAssemblerFPI<dim>>(
    this->simulation_control,
    this->cfd_dem_simulation_parameters.cfd_parameters.physical_properties,
    this->cfd_dem_simulation_parameters.cfd_dem));

  // The core assembler should always be the last assembler to be called in the
  // stabilized formulation as to have all strong residual and jacobian stored.
  // Core assembler
  this->assemblers.push_back(std::make_shared<GLSVansAssemblerCore<dim>>(
    this->simulation_control,
    this->cfd_dem_simulation_parameters.cfd_parameters.physical_properties,
    this->cfd_dem_simulation_parameters.cfd_dem));
}

template <int dim>
void
GLSVANSSolver<dim>::assemble_system_matrix()
{
  this->system_matrix = 0;
  this->simulation_control->set_assembly_method(this->time_stepping_method);

  setup_assemblers();

  auto scratch_data = NavierStokesScratchData<dim>(*this->fe,
                                                   *this->cell_quadrature,
                                                   *this->mapping);

  scratch_data.enable_void_fraction(fe_void_fraction,
                                    *this->cell_quadrature,
                                    *this->mapping);


  scratch_data.enable_particle_fluid_interactions(
    particle_handler.n_global_max_particles_per_cell());

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
  this->setup_preconditioner();
}

template <int dim>
void
GLSVANSSolver<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim> &                        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &                copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      this->solution_stages,
                      this->forcing_function,
                      this->beta);

  typename DoFHandler<dim>::active_cell_iterator void_fraction_cell(
    &(*(this->triangulation)),
    cell->level(),
    cell->index(),
    &this->void_fraction_dof_handler);

  scratch_data.reinit_void_fraction(
    void_fraction_cell,
    nodal_void_fraction_relevant,
    previous_void_fraction,
    std::vector<TrilinosWrappers::MPI::Vector>());


  scratch_data.reinit_particle_fluid_interactions(this->previous_solutions[0],
                                                  nodal_void_fraction_relevant,
                                                  particle_handler,
                                                  this->dof_handler,
                                                  void_fraction_dof_handler);
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
  this->system_rhs = 0;
  this->simulation_control->set_assembly_method(this->time_stepping_method);
  setup_assemblers();

  auto scratch_data = NavierStokesScratchData<dim>(*this->fe,
                                                   *this->cell_quadrature,
                                                   *this->mapping);


  scratch_data.enable_void_fraction(fe_void_fraction,
                                    *this->cell_quadrature,
                                    *this->mapping);


  scratch_data.enable_particle_fluid_interactions(
    particle_handler.n_global_max_particles_per_cell());

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
  NavierStokesScratchData<dim> &                        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &                copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      this->solution_stages,
                      this->forcing_function,
                      this->beta);

  typename DoFHandler<dim>::active_cell_iterator void_fraction_cell(
    &(*(this->triangulation)),
    cell->level(),
    cell->index(),
    &this->void_fraction_dof_handler);


  scratch_data.reinit_void_fraction(
    void_fraction_cell,
    nodal_void_fraction_relevant,
    previous_void_fraction,
    std::vector<TrilinosWrappers::MPI::Vector>());

  scratch_data.reinit_particle_fluid_interactions(this->previous_solutions[0],
                                                  nodal_void_fraction_relevant,
                                                  particle_handler,
                                                  this->dof_handler,
                                                  void_fraction_dof_handler);

  copy_data.reset();

  for (auto &pf_assembler : particle_fluid_assemblers)
    {
      pf_assembler->calculate_particle_fluid_interactions(scratch_data);
    }

  for (unsigned int counter = 0; counter < scratch_data.particle_index;
       ++counter)

    {
      fluid_solid_force[scratch_data.local_particle_id[counter]] =
        scratch_data.fluid_particle_force[counter];
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
GLSVANSSolver<dim>::post_processing()
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

  //  std::vector<double> present_pressure_values(n_q_points);

  double mass_source           = 0;
  double fluid_volume          = 0;
  double bed_volume            = 0;
  double average_void_fraction = 0;
  double pressure_drop         = 0;

  Vector<double>      bdf_coefs;
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

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

          // Gather the previous time steps depending on the number of stages
          // of the time integration scheme for the void fraction

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

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // Calculate the divergence of the velocity
              const double present_velocity_divergence =
                trace(present_velocity_gradients[q]);

              // Evaluation of global mass conservation
              mass_source += (present_velocity_values[q] *
                                present_void_fraction_gradients[q] +
                              present_void_fraction_values[q] *
                                present_velocity_divergence) *
                             fe_values_void_fraction.JxW(q);

              if (scheme ==
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
                  scheme == Parameters::SimulationControl::TimeSteppingMethod::
                              steady_bdf)
                mass_source += (bdf_coefs[0] * present_void_fraction_values[q] +
                                bdf_coefs[1] * p1_void_fraction_values[q]) *
                               fe_values_void_fraction.JxW(q);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                mass_source += (bdf_coefs[0] * present_void_fraction_values[q] +
                                bdf_coefs[1] * p1_void_fraction_values[q] +
                                bdf_coefs[2] * p2_void_fraction_values[q]) *
                               fe_values_void_fraction.JxW(q);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                mass_source += (bdf_coefs[0] * present_void_fraction_values[q] +
                                bdf_coefs[1] * p1_void_fraction_values[q] +
                                bdf_coefs[2] * p2_void_fraction_values[q] +
                                bdf_coefs[3] * p3_void_fraction_values[q]) *
                               fe_values_void_fraction.JxW(q);

              // Calculation of fluid and bed volumes in bed
              if (present_void_fraction_values[q] < 0.6)
                {
                  fluid_volume += present_void_fraction_values[q] *
                                  fe_values_void_fraction.JxW(q);

                  bed_volume += fe_values_void_fraction.JxW(q);
                }
            }
        }
    }

  mass_source  = Utilities::MPI::sum(mass_source, this->mpi_communicator);
  fluid_volume = Utilities::MPI::sum(fluid_volume, this->mpi_communicator);
  bed_volume   = Utilities::MPI::sum(bed_volume, this->mpi_communicator);

  average_void_fraction = fluid_volume / bed_volume;

  pressure_drop = calculate_pressure_drop<dim>(
    this->dof_handler,
    this->mapping,
    this->mpi_communicator,
    this->fe,
    this->evaluation_point,
    this->number_quadrature_points,
    this->cfd_dem_simulation_parameters.cfd_dem.inlet_boundary_id,
    this->cfd_dem_simulation_parameters.cfd_dem.outlet_boundary_id);

  this->pcout << "Mass Source: " << mass_source << " s^-1" << std::endl;
  this->pcout << "Average Void Fraction in Bed: " << average_void_fraction
              << std::endl;
  this->pcout << "Pressure Drop: " << pressure_drop << " m^2.s^-2" << std::endl;
}

template <int dim>
void
GLSVANSSolver<dim>::solve()
{
  // This is enforced to 1 right now because it does not provide
  // better speed-up than using MPI. This could be eventually changed...
  MultithreadInfo::set_thread_limit(1);

  read_mesh_and_manifolds(
    this->triangulation,
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
  calculate_void_fraction(this->simulation_control->get_current_time());
  this->set_initial_condition(
    this->cfd_dem_simulation_parameters.cfd_parameters.initial_condition->type,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
      .restart);

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      if (this->simulation_control->is_at_start())
        {
          this->first_iteration();
        }
      else
        {
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          this->iterate();
        }

      this->postprocess(false);
      this->finish_time_step();

      if (this->cfd_dem_simulation_parameters.cfd_dem.post_processing)
        post_processing();
    }

  this->finish_simulation();
}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the
// library is valid before we actually compile the solver This greatly
// helps with debugging
template class GLSVANSSolver<2>;
template class GLSVANSSolver<3>;
