#include <fem-dem/gls_vans.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>

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
GLSVANSSolver<dim>::GLSVANSSolver(SimulationParameters<dim> &p_nsparam)
  : GLSNavierStokesSolver<dim>(p_nsparam)
  , void_fraction_dof_handler(*this->triangulation)
  , fe_void_fraction(
      this->simulation_parameters.fem_parameters.void_fraction_order)
  , particle_mapping(1)
  , particle_handler(*this->triangulation,
                     particle_mapping,
                     DEM::get_number_properties())
{}

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
  void_fraction_m1.reinit(locally_owned_dofs_voidfraction,
                          locally_relevant_dofs_voidfraction,
                          this->mpi_communicator);
  void_fraction_m2.reinit(locally_owned_dofs_voidfraction,
                          locally_relevant_dofs_voidfraction,
                          this->mpi_communicator);
  void_fraction_m3.reinit(locally_owned_dofs_voidfraction,
                          locally_relevant_dofs_voidfraction,
                          this->mpi_communicator);
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
GLSVANSSolver<dim>::finish_time_step_fd()
{
  GLSNavierStokesSolver<dim>::finish_time_step_fd();

  void_fraction_m3 = void_fraction_m2;
  void_fraction_m2 = void_fraction_m1;
  void_fraction_m1 = nodal_void_fraction_relevant;
}

template <int dim>
void
GLSVANSSolver<dim>::read_dem()
{
  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  std::string prefix = this->simulation_parameters.void_fraction->dem_file_name;

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
  void_fraction_m3 = nodal_void_fraction_relevant;
  void_fraction_m2 = nodal_void_fraction_relevant;
  void_fraction_m1 = nodal_void_fraction_relevant;
}

template <int dim>
void
GLSVANSSolver<dim>::calculate_void_fraction(const double time)
{
  if (this->simulation_parameters.void_fraction->mode ==
      Parameters::VoidFractionMode::function)
    {
      const MappingQ<dim> mapping(1);

      this->simulation_parameters.void_fraction->void_fraction.set_time(time);

      VectorTools::interpolate(
        mapping,
        void_fraction_dof_handler,
        this->simulation_parameters.void_fraction->void_fraction,
        nodal_void_fraction_owned);

      nodal_void_fraction_relevant = nodal_void_fraction_owned;
    }
  else if (this->simulation_parameters.void_fraction->mode ==
           Parameters::VoidFractionMode::dem)
    {
      assemble_L2_projection_void_fraction();
      solve_L2_system_void_fraction();
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
                          (solution_value - this->simulation_parameters
                                              .void_fraction->l2_upper_bound) >
                      0)
                    {
                      active_set.add_index(dof_index);
                      void_fraction_constraints.add_line(dof_index);
                      void_fraction_constraints.set_inhomogeneity(
                        dof_index,
                        this->simulation_parameters.void_fraction
                          ->l2_upper_bound);
                      nodal_void_fraction_owned(dof_index) =
                        this->simulation_parameters.void_fraction
                          ->l2_upper_bound;
                      lambda(dof_index) = 0;
                    }
                  else if (lambda(dof_index) +
                             penalty_parameter *
                               mass_matrix(dof_index, dof_index) *
                               (solution_value -
                                this->simulation_parameters.void_fraction
                                  ->l2_lower_bound) <
                           0)
                    {
                      active_set.add_index(dof_index);
                      void_fraction_constraints.add_line(dof_index);
                      void_fraction_constraints.set_inhomogeneity(
                        dof_index,
                        this->simulation_parameters.void_fraction
                          ->l2_lower_bound);
                      nodal_void_fraction_owned(dof_index) =
                        this->simulation_parameters.void_fraction
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
  const MappingQ<dim> mapping(
    1, this->simulation_parameters.fem_parameters.qmapping_all);

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
                        (this->simulation_parameters.void_fraction
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

  if (this->simulation_parameters.linear_solver.verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const IndexSet locally_owned_dofs_voidfraction =
    void_fraction_dof_handler.locally_owned_dofs();

  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    locally_owned_dofs_voidfraction, this->mpi_communicator);


  SolverControl solver_control(
    this->simulation_parameters.linear_solver.max_iterations,
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverCG solver(solver_control);

  TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill =
    this->simulation_parameters.linear_solver.ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.ilu_precond_rtol;

  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_matrix_void_fraction,
                                 preconditionerOptions);

  solver.solve(system_matrix_void_fraction,
               completely_distributed_solution,
               system_rhs_void_fraction,
               *ilu_preconditioner);

  if (this->simulation_parameters.linear_solver.verbosity !=
      Parameters::Verbosity::quiet)
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
  if (!is_bdf_high_order(this->simulation_parameters.simulation_control.method))
    {
      iterate();
    }

  // Taking care of the multi-step methods
  else if (this->simulation_parameters.simulation_control.method ==
           Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    {
      Parameters::SimulationControl timeParameters =
        this->simulation_parameters.simulation_control;

      // Start the BDF2 with a single Euler time step with a lower time step
      double time_step =
        timeParameters.dt * timeParameters.startup_timestep_scaling;
      this->simulation_control->set_current_time_step(time_step);

      double intermediate_time =
        this->simulation_control->get_current_time() + time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);


      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf1, false, true);
      this->percolate_time_vectors_fd();
      void_fraction_m2 = void_fraction_m1;
      void_fraction_m1 = nodal_void_fraction_relevant;

      // Reset the time step and do a bdf 2 newton iteration using the two
      // steps to complete the full step

      time_step =
        timeParameters.dt * (1. - timeParameters.startup_timestep_scaling);

      this->simulation_control->set_current_time_step(time_step);
      intermediate_time += time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);


      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf2, false, true);

      this->simulation_control->set_suggested_time_step(timeParameters.dt);
    }

  else if (this->simulation_parameters.simulation_control.method ==
           Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    {
      Parameters::SimulationControl timeParameters =
        this->simulation_parameters.simulation_control;

      // Start the BDF3 with a single Euler time step with a lower time step
      double time_step =
        timeParameters.dt * timeParameters.startup_timestep_scaling;

      this->simulation_control->set_current_time_step(time_step);

      double intermediate_time = time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);


      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf1, false, true);
      this->percolate_time_vectors_fd();
      void_fraction_m2 = void_fraction_m1;
      void_fraction_m1 = nodal_void_fraction_relevant;

      // Reset the time step and do a bdf 2 newton iteration using the two
      // steps

      this->simulation_control->set_current_time_step(time_step);
      intermediate_time += time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);


      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf1, false, true);
      this->percolate_time_vectors_fd();
      void_fraction_m3 = void_fraction_m2;
      void_fraction_m2 = void_fraction_m1;
      void_fraction_m1 = nodal_void_fraction_relevant;

      // Reset the time step and do a bdf 3 newton iteration using the two
      // steps to complete the full step
      time_step =
        timeParameters.dt * (1. - 2. * timeParameters.startup_timestep_scaling);
      this->simulation_control->set_current_time_step(time_step);
      intermediate_time += time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);


      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf3, false, true);
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
    this->simulation_parameters.simulation_control.method, false, false);
}


template <int dim>
template <bool                                              assemble_matrix,
          Parameters::SimulationControl::TimeSteppingMethod scheme,
          Parameters::VelocitySource::VelocitySourceType    velocity_source>
void
GLSVANSSolver<dim>::assembleGLS()
{
  if (assemble_matrix)
    this->system_matrix = 0;
  auto &system_rhs = this->system_rhs;
  system_rhs       = 0;

  double viscosity = this->simulation_parameters.physical_properties.viscosity;
  Function<dim> *l_forcing_function = this->forcing_function;

  QGauss<dim>   quadrature_formula(this->number_quadrature_points);
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

  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int void_fraction_dofs_per_cell =
    this->fe_void_fraction.dofs_per_cell;
  const unsigned int               n_q_points = quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  FullMatrix<double>               local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>                   local_rhs(dofs_per_cell);
  std::vector<Vector<double>> rhs_force(n_q_points, Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<Tensor<1, dim>>          present_velocity_values(n_q_points);
  std::vector<Tensor<2, dim>>          present_velocity_gradients(n_q_points);
  std::vector<Tensor<2, dim>> present_velocity_gradients_transpose(n_q_points);
  std::vector<Tensor<2, dim>> stress_tensor_velocity(n_q_points);
  std::vector<double>         present_pressure_values(n_q_points);
  std::vector<Tensor<1, dim>> present_pressure_gradients(n_q_points);
  std::vector<Tensor<1, dim>> present_velocity_laplacians(n_q_points);
  std::vector<Tensor<2, dim>> present_velocity_hess(n_q_points);

  // Data storage vector for the void fraction values/gradients
  std::vector<double>         present_void_fraction_values(n_q_points);
  std::vector<Tensor<1, dim>> present_void_fraction_gradients(n_q_points);

  Tensor<1, dim> force;

  // Drag variables
  double                               beta;
  std::vector<types::global_dof_index> fluid_dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> void_fraction_dof_indices(
    void_fraction_dofs_per_cell);
  Tensor<1, dim> particle_velocity;
  Tensor<1, dim> average_particle_velocity;
  Tensor<1, dim> fluid_velocity_at_particle_location;
  Tensor<1, dim> relative_velocity;
  double         c_d                = 0;
  double         cell_void_fraction = 0;
  int            particles_in_cell;


  // Velocity dependent source term
  //----------------------------------
  // Angular velocity of the rotating frame. This is always a 3D vector even
  // in 2D.
  Tensor<1, dim> omega_vector;

  double omega_z  = this->simulation_parameters.velocity_sources.omega_z;
  omega_vector[0] = this->simulation_parameters.velocity_sources.omega_x;
  omega_vector[1] = this->simulation_parameters.velocity_sources.omega_y;
  if (dim == 3)
    omega_vector[2] = this->simulation_parameters.velocity_sources.omega_z;

  std::vector<double>         div_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<Tensor<3, dim>> hess_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>> laplacian_phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> stress_tensor_phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi_u_transpose(dofs_per_cell);
  std::vector<double>         phi_p(dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi_p(dofs_per_cell);

  // Shape Function tensor for drag calculation
  std::vector<Tensor<1, dim>> shape_function(dofs_per_cell);

  // Values at previous time step for transient schemes
  std::vector<Tensor<1, dim>> p1_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> p2_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> p3_velocity_values(n_q_points);

  // Values at previous time step for transient schemes for void fraction
  std::vector<double> p1_void_fraction_values(n_q_points);
  std::vector<double> p2_void_fraction_values(n_q_points);
  std::vector<double> p3_void_fraction_values(n_q_points);

  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Time steps and inverse time steps which is used for numerous calculations
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Vector for the BDF coefficients
  // The coefficients are stored in the following fashion :
  // 0 - n+1
  // 1 - n
  // 2 - n-1
  // 3 - n-2
  Vector<double> bdf_coefs;

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf1)
    bdf_coefs = bdf_coefficients(1, time_steps_vector);

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    bdf_coefs = bdf_coefficients(2, time_steps_vector);

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    bdf_coefs = bdf_coefficients(3, time_steps_vector);

  // Element size
  double h;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          typename DoFHandler<dim>::active_cell_iterator void_fraction_cell(
            &(*this->triangulation),
            cell->level(),
            cell->index(),
            &this->void_fraction_dof_handler);
          fe_values_void_fraction.reinit(void_fraction_cell);

          if (dim == 2)
            h = std::sqrt(4. * cell->measure() / M_PI) /
                this->velocity_fem_degree;
          else if (dim == 3)
            h = pow(6 * cell->measure() / M_PI, 1. / 3.) /
                this->velocity_fem_degree;

          local_matrix              = 0;
          local_rhs                 = 0;
          beta                      = 0;
          average_particle_velocity = 0;
          particles_in_cell         = 0;

          // Gather velocity (values, gradient and laplacian)
          auto &evaluation_point = this->evaluation_point;
          fe_values[velocities].get_function_values(evaluation_point,
                                                    present_velocity_values);
          fe_values[velocities].get_function_gradients(
            evaluation_point, present_velocity_gradients);
          fe_values[velocities].get_function_laplacians(
            evaluation_point, present_velocity_laplacians);

          // Gather pressure (values, gradient)
          fe_values[pressure].get_function_values(evaluation_point,
                                                  present_pressure_values);
          fe_values[pressure].get_function_gradients(
            evaluation_point, present_pressure_gradients);

          // Gather void fraction (values, gradient)
          fe_values_void_fraction.get_function_values(
            nodal_void_fraction_relevant, present_void_fraction_values);
          fe_values_void_fraction.get_function_gradients(
            nodal_void_fraction_relevant, present_void_fraction_gradients);

          std::vector<Point<dim>> quadrature_points =
            fe_values.get_quadrature_points();

          // Calculate forcing term if there is a forcing function
          if (l_forcing_function)
            l_forcing_function->vector_value_list(quadrature_points, rhs_force);

          // Calculate Transpose of the present_velocity_gradients tensor
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              present_velocity_gradients_transpose[q] =
                transpose(present_velocity_gradients[q]);
            }

          // Calculate the stress tensor from present_velocity_gradients
          if (this->simulation_parameters.cfd_dem.full_stress_tensor == true)
            {
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  stress_tensor_velocity[q] =
                    present_velocity_gradients[q] +
                    present_velocity_gradients_transpose[q];

                  for (int d = 0; d < dim; ++d)
                    {
                      stress_tensor_velocity[q][d][d] +=
                        -(2.0 / 3) *
                        trace(present_velocity_gradients[q]); // Divergence of u
                    }
                }
            }
          else
            stress_tensor_velocity = present_velocity_gradients;


          // Gather the previous time steps depending on the number of stages
          // of the time integration scheme
          if (scheme !=
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            fe_values[velocities].get_function_values(
              this->previous_solutions[0], p1_velocity_values);

          if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf2)
            fe_values[velocities].get_function_values(
              this->previous_solutions[1], p2_velocity_values);

          if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf3)
            fe_values[velocities].get_function_values(
              this->previous_solutions[2], p3_velocity_values);

          // Gather the previous time steps depending on the number of stages
          // of the time integration scheme for the void fraction

          if (scheme !=
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            {
              fe_values_void_fraction.get_function_values(
                void_fraction_m1, p1_void_fraction_values);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                fe_values_void_fraction.get_function_values(
                  void_fraction_m2, p2_void_fraction_values);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                fe_values_void_fraction.get_function_values(
                  void_fraction_m3, p3_void_fraction_values);
            }

          //*******************************************************************
          // Calculation of particles' drag force
          if (this->simulation_parameters.void_fraction->mode ==
              Parameters::VoidFractionMode::dem)
            {
              const auto pic = particle_handler.particles_in_cell(cell);

              double cell_volume = cell->measure();

              // Get the local dof indices for velocity and void fraction field
              // interpolation at particle's position
              const auto &dh_cell =
                typename DoFHandler<dim>::cell_iterator(*cell,
                                                        &this->dof_handler);
              dh_cell->get_dof_indices(fluid_dof_indices);

              const auto &void_fraction_dh_cell =
                typename DoFHandler<dim>::cell_iterator(
                  *void_fraction_cell, &void_fraction_dof_handler);
              void_fraction_dh_cell->get_dof_indices(void_fraction_dof_indices);

              // Loop over particles in cell
              for (auto &particle : pic)
                {
                  auto particle_properties = particle.get_properties();

                  // Stock the values of particle velocity in a
                  // tensor
                  particle_velocity[0] =
                    particle_properties[DEM::PropertiesIndex::v_x];
                  particle_velocity[1] =
                    particle_properties[DEM::PropertiesIndex::v_y];
                  if (dim == 3)
                    particle_velocity[2] =
                      particle_properties[DEM::PropertiesIndex::v_z];

                  // Interpolate velocity and void fraction at particle position
                  // Reference location of the particle
                  auto reference_location = particle.get_reference_location();

                  fluid_velocity_at_particle_location = 0;
                  cell_void_fraction                  = 0;

                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      const auto comp_j =
                        this->fe->system_to_component_index(j).first;
                      if (comp_j < dim)
                        {
                          auto &evaluation_point = this->previous_solutions[0];
                          fluid_velocity_at_particle_location[comp_j] +=
                            evaluation_point[fluid_dof_indices[j]] *
                            this->fe->shape_value(j, reference_location);
                        }
                    }
                  for (unsigned int j = 0; j < void_fraction_dofs_per_cell; ++j)
                    {
                      cell_void_fraction +=
                        nodal_void_fraction_relevant
                          [void_fraction_dof_indices[j]] *
                        this->fe_void_fraction.shape_value(j,
                                                           reference_location);
                    }

                  average_particle_velocity += particle_velocity;
                  particles_in_cell += 1;
                  relative_velocity =
                    fluid_velocity_at_particle_location - particle_velocity;

                  // Particle's Reynolds number
                  double re =
                    1e-1 + relative_velocity.norm() *
                             particle_properties[DEM::PropertiesIndex::dp] /
                             viscosity;

                  // Drag Coefficient Calculation
                  if (this->simulation_parameters.cfd_dem.drag_model ==
                      Parameters::DragModel::difelice)
                    // Di Felice Drag Model CD Calculation
                    {
                      c_d = pow((0.63 + 4.8 / sqrt(re)), 2) *
                            pow(cell_void_fraction,
                                -(3.7 -
                                  0.65 * exp(-pow((1.5 - log10(re)), 2) / 2)));
                    }
                  // Rong Drag Model CD Calculation
                  else if (this->simulation_parameters.cfd_dem.drag_model ==
                           Parameters::DragModel::rong)
                    {
                      c_d = pow((0.63 + 4.8 / sqrt(re)), 2) *
                            pow(cell_void_fraction,
                                -(2.65 * (cell_void_fraction + 1) -
                                  (5.3 - (3.5 * cell_void_fraction)) *
                                    pow(cell_void_fraction, 2) *
                                    exp(-pow(1.5 - log10(re), 2) / 2)));
                    }

                  beta +=
                    (0.5 * c_d * M_PI *
                     pow(particle_properties[DEM::PropertiesIndex::dp], 2) /
                     4) *
                    relative_velocity.norm();
                }

              if (particles_in_cell != 0)
                {
                  average_particle_velocity =
                    average_particle_velocity / particles_in_cell;
                }
              beta = beta / cell_volume;
            }

          //**************************************************************************

          // Loop over the quadrature points
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // Calculation of the magnitude of the velocity for the
              // stabilization parameter
              const double u_mag = std::max(present_velocity_values[q].norm(),
                                            1e-12 * GLS_u_scale);

              // Store JxW in local variable for faster access;
              const double JxW = fe_values.JxW(q);

              // Calculation of the GLS stabilization parameter. The
              // stabilization parameter used is different if the simulation
              // is steady or unsteady. In the unsteady case it includes the
              // value of the time-step
              const double tau =
                scheme ==
                    Parameters::SimulationControl::TimeSteppingMethod::steady ?
                  1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                                 9 * std::pow(4 * viscosity / (h * h), 2)) :
                  1. /
                    std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                              9 * std::pow(4 * viscosity / (h * h), 2));

              // Calcukation of the shock capturing viscosity term
              const double vdcdd =
                (h /
                 (2 * this->simulation_parameters.cfd_dem.reference_velocity)) *
                pow(present_velocity_gradients[q].norm() * h /
                      (this->simulation_parameters.cfd_dem.reference_velocity),
                    this->simulation_parameters.fem_parameters.velocity_order) *
                pow((this->simulation_parameters.cfd_dem.reference_velocity /
                     this->simulation_parameters.void_fraction->l2_lower_bound),
                    2);

              // Grad-div weight factor
              const double gamma = 0.1;


              // Gather the shape functions, their gradient and their
              // laplacian for the velocity and the pressure
              for (unsigned int k = 0; k < dofs_per_cell; ++k)

                {
                  div_phi_u[k]  = fe_values[velocities].divergence(k, q);
                  grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                  phi_u[k]      = fe_values[velocities].value(k, q);
                  hess_phi_u[k] = fe_values[velocities].hessian(k, q);
                  phi_p[k]      = fe_values[pressure].value(k, q);
                  grad_phi_p[k] = fe_values[pressure].gradient(k, q);

                  for (int d = 0; d < dim; ++d)
                    laplacian_phi_u[k][d] = trace(hess_phi_u[k][d]);
                }

              // Calculate Transpose of the grad_phi_u tensor
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  grad_phi_u_transpose[i] = transpose(grad_phi_u[i]);
                }

              // Calculate the stress tensor from grad_phi_u
              if (this->simulation_parameters.cfd_dem.full_stress_tensor ==
                  true)
                {
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      stress_tensor_phi_u[i] =
                        grad_phi_u[i] + grad_phi_u_transpose[i];

                      for (int d = 0; d < dim; ++d)
                        {
                          stress_tensor_phi_u[i][d][d] +=
                            -(2.0 / 3) * div_phi_u[i];
                        }
                    }
                }
              else
                stress_tensor_phi_u = grad_phi_u;

              // Establish the force vector
              for (int i = 0; i < dim; ++i)
                {
                  const unsigned int component_i =
                    this->fe->system_to_component_index(i).first;
                  force[i] = rhs_force[q](component_i);
                }
              const unsigned int component_mass =
                this->fe->system_to_component_index(dim).first;
              double mass_source = rhs_force[q](component_mass);

              // Calculate the divergence of the velocity
              const double present_velocity_divergence =
                trace(present_velocity_gradients[q]);

              // Calculate the strong residual for GLS stabilization
              auto strong_residual =
                present_velocity_gradients[q] * present_velocity_values[q] *
                  present_void_fraction_values[q]
                // Mass source term
                + mass_source * present_velocity_values[q] +
                // Drag Force
                beta *
                  (present_velocity_values[q] - average_particle_velocity) +
                present_pressure_gradients[q] -
                viscosity * present_velocity_laplacians[q] -
                force * present_void_fraction_values[q];


              if (velocity_source ==
                  Parameters::VelocitySource::VelocitySourceType::srf)
                {
                  if (dim == 2)
                    {
                      strong_residual +=
                        2 * omega_z * (-1.) *
                        cross_product_2d(present_velocity_values[q]);
                      auto centrifugal =
                        omega_z * (-1.) *
                        cross_product_2d(
                          omega_z * (-1.) *
                          cross_product_2d(quadrature_points[q]));
                      strong_residual += centrifugal;
                    }
                  else // dim == 3
                    {
                      strong_residual +=
                        2 * cross_product_3d(omega_vector,
                                             present_velocity_values[q]);
                      strong_residual += cross_product_3d(
                        omega_vector,
                        cross_product_3d(omega_vector, quadrature_points[q]));
                    }
                }

              if (this->simulation_parameters.cfd_dem.shock_capturing == true)
                strong_residual += -vdcdd * present_velocity_laplacians[q];

              /* Adjust the strong residual in cases where the scheme is
               transient.
               The BDF schemes require values at previous time steps which
               are stored in the p1, p2 and p3 vectors.
               */

              if (scheme ==
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
                  scheme == Parameters::SimulationControl::TimeSteppingMethod::
                              steady_bdf)
                strong_residual += (bdf_coefs[0] * present_velocity_values[q] +
                                    bdf_coefs[1] * p1_velocity_values[q]) *
                                   present_void_fraction_values[q];


              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                strong_residual += (bdf_coefs[0] * present_velocity_values[q] +
                                    bdf_coefs[1] * p1_velocity_values[q] +
                                    bdf_coefs[2] * p2_velocity_values[q]) *
                                   present_void_fraction_values[q];

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                strong_residual += (bdf_coefs[0] * present_velocity_values[q] +
                                    bdf_coefs[1] * p1_velocity_values[q] +
                                    bdf_coefs[2] * p2_velocity_values[q] +
                                    bdf_coefs[3] * p3_velocity_values[q]) *
                                   present_void_fraction_values[q];



              // Matrix assembly
              if (assemble_matrix)
                {
                  // We loop over the column first to prevent recalculation
                  // of the strong jacobian in the inner loop
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      auto strong_jac =
                        (present_velocity_gradients[q] * phi_u[j] *
                           present_void_fraction_values[q] +
                         grad_phi_u[j] * present_velocity_values[q] *
                           present_void_fraction_values[q]
                         // Mass source term
                         + mass_source * phi_u[j]
                         // Drag Force
                         + beta * phi_u[j] + grad_phi_p[j] -
                         viscosity * laplacian_phi_u[j]);

                      if (is_bdf(scheme))
                        strong_jac += present_void_fraction_values[q] *
                                      phi_u[j] * bdf_coefs[0];

                      if (velocity_source ==
                          Parameters::VelocitySource::VelocitySourceType::srf)
                        {
                          if (dim == 2)
                            strong_jac +=
                              2 * omega_z * (-1.) * cross_product_2d(phi_u[j]);
                          else if (dim == 3)
                            strong_jac +=
                              2 * cross_product_3d(omega_vector, phi_u[j]);
                        }
                      if (this->simulation_parameters.cfd_dem.shock_capturing ==
                          true)
                        strong_jac += -vdcdd * laplacian_phi_u[j];

                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                          local_matrix(i, j) +=
                            (
                              // Momentum terms
                              viscosity * scalar_product(stress_tensor_phi_u[j],
                                                         grad_phi_u[i]) +
                              // Advection terms
                              ((phi_u[j] * present_void_fraction_values[q] *
                                present_velocity_gradients[q] * phi_u[i]) +
                               (grad_phi_u[j] *
                                present_void_fraction_values[q] *
                                present_velocity_values[q] * phi_u[i]))
                              // Mass source term
                              + mass_source * phi_u[j] * phi_u[i]
                              // Pressure
                              - (div_phi_u[i] * phi_p[j]) +
                              // Drag force
                              beta * phi_u[j] * phi_u[i] +
                              // Continuity
                              phi_p[i] *
                                ((present_void_fraction_values[q] *
                                  div_phi_u[j]) +
                                 (phi_u[j] *
                                  present_void_fraction_gradients[q]))) *
                            JxW;

                          // Grad-div stabilization - Bruno test
                          if (this->simulation_parameters.cfd_dem.grad_div ==
                              true)
                            {
                              local_matrix(i, j) +=
                                gamma *
                                (div_phi_u[j] *
                                   present_void_fraction_values[q] +
                                 phi_u[j] *
                                   present_void_fraction_gradients[q]) *
                                div_phi_u[i] * JxW;
                            }

                          // Mass matrix
                          if (is_bdf(scheme))
                            {
                              local_matrix(i, j) +=
                                present_void_fraction_values[q] * phi_u[j] *
                                phi_u[i] * bdf_coefs[0] * JxW;
                            }

                          // PSPG GLS term
                          if (PSPG)
                            local_matrix(i, j) +=
                              tau * strong_jac * grad_phi_p[i] * JxW;

                          if (velocity_source == Parameters::VelocitySource::
                                                   VelocitySourceType::srf)
                            {
                              if (dim == 2)
                                local_matrix(i, j) +=
                                  2 * omega_z * (-1.) *
                                  cross_product_2d(phi_u[j]) * phi_u[i] * JxW;

                              else if (dim == 3)
                                local_matrix(i, j) +=
                                  2 * cross_product_3d(omega_vector, phi_u[j]) *
                                  phi_u[i] * JxW;
                            }


                          // PSPG TAU term is currently disabled because it
                          // does not alter the matrix sufficiently
                          // local_matrix(i, j) +=
                          //  -tau * tau * tau * 4 / h / h *
                          //  (present_velocity_values[q] * phi_u[j]) *
                          //  strong_residual * grad_phi_p[i] *
                          //  fe_values.JxW(q);

                          // Jacobian is currently incomplete
                          if (SUPG)
                            {
                              local_matrix(i, j) +=
                                tau *
                                (strong_jac * (grad_phi_u[i] *
                                               present_velocity_values[q]) +
                                 strong_residual * (grad_phi_u[i] * phi_u[j])) *
                                JxW;

                              // SUPG TAU term is currently disabled because
                              // it does not alter the matrix sufficiently
                              // local_matrix(i, j)
                              // +=
                              //   -strong_residual
                              //   * (grad_phi_u[i]
                              //   *
                              //   present_velocity_values[q])
                              //   * tau * tau *
                              //   tau * 4 / h / h
                              //   *
                              //   (present_velocity_values[q]
                              //   * phi_u[j]) *
                              //   fe_values.JxW(q);
                            }

                          if (this->simulation_parameters.cfd_dem
                                .shock_capturing == true)
                            {
                              local_matrix(i, j) +=
                                vdcdd *
                                scalar_product(grad_phi_u[j], grad_phi_u[i]) *
                                JxW;
                            }
                        }
                    }
                }

              // Assembly of the right-hand side
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Navier-Stokes Residual
                  local_rhs(i) +=
                    (
                      // Momentum
                      -viscosity * scalar_product(stress_tensor_velocity[q],
                                                  grad_phi_u[i]) -
                      // Advection terms
                      (present_velocity_gradients[q] *
                       present_velocity_values[q] *
                       present_void_fraction_values[q] * phi_u[i])
                      // Mass source term
                      - mass_source * present_velocity_values[q] * phi_u[i]
                      // Pressure and force
                      + present_pressure_values[q] * div_phi_u[i] +
                      force * present_void_fraction_values[q] * phi_u[i]
                      // Drag Force
                      - beta *
                          (present_velocity_values[q] -
                           average_particle_velocity) *
                          phi_u[i] -
                      // Continuity
                      (present_velocity_divergence *
                         present_void_fraction_values[q] +
                       present_velocity_values[q] *
                         present_void_fraction_gradients[q] -
                       mass_source) *
                        phi_p[i]) *
                    JxW;


                  // Residual associated with BDF schemes
                  if (scheme == Parameters::SimulationControl::
                                  TimeSteppingMethod::bdf1 ||
                      scheme == Parameters::SimulationControl::
                                  TimeSteppingMethod::steady_bdf)
                    {
                      local_rhs(i) -=
                        (bdf_coefs[0] * present_velocity_values[q] +
                         bdf_coefs[1] * p1_velocity_values[q]) *
                        present_void_fraction_values[q] * phi_u[i] * JxW;

                      local_rhs(i) -=
                        (bdf_coefs[0] * present_void_fraction_values[q] +
                         bdf_coefs[1] * p1_void_fraction_values[q]) *
                        phi_p[i] * JxW;

                      if (this->simulation_parameters.cfd_dem.grad_div == true)
                        {
                          local_rhs(i) -=
                            (bdf_coefs[0] * present_void_fraction_values[q] +
                             bdf_coefs[1] * p1_void_fraction_values[q]) *
                            div_phi_u[i] * JxW;
                        }
                    }

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                    {
                      local_rhs(i) -=
                        (bdf_coefs[0] * present_velocity_values[q] +
                         bdf_coefs[1] * p1_velocity_values[q] +
                         bdf_coefs[2] * p2_velocity_values[q]) *
                        present_void_fraction_values[q] * phi_u[i] * JxW;

                      local_rhs(i) -=
                        (bdf_coefs[0] * present_void_fraction_values[q] +
                         bdf_coefs[1] * p1_void_fraction_values[q] +
                         bdf_coefs[2] * p2_void_fraction_values[q]) *
                        phi_p[i] * JxW;

                      if (this->simulation_parameters.cfd_dem.grad_div == true)
                        {
                          local_rhs(i) -=
                            (bdf_coefs[0] * present_void_fraction_values[q] +
                             bdf_coefs[1] * p1_void_fraction_values[q] +
                             bdf_coefs[2] * p2_void_fraction_values[q]) *
                            div_phi_u[i] * JxW;
                        }
                    }


                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                    {
                      local_rhs(i) -=

                        local_rhs(i) -=
                        (bdf_coefs[0] * present_velocity_values[q] +
                         bdf_coefs[1] * p1_velocity_values[q] +
                         bdf_coefs[2] * p2_velocity_values[q] +
                         bdf_coefs[3] * p3_velocity_values[q]) *
                        present_void_fraction_values[q] * phi_u[i] * JxW;

                      local_rhs(i) -=
                        (bdf_coefs[0] * present_void_fraction_values[q] +
                         bdf_coefs[1] * p1_void_fraction_values[q] +
                         bdf_coefs[2] * p2_void_fraction_values[q] +
                         bdf_coefs[3] * p3_void_fraction_values[q]) *
                        phi_p[i] * JxW;

                      if (this->simulation_parameters.cfd_dem.grad_div == true)
                        {
                          local_rhs(i) -=
                            (bdf_coefs[0] * present_void_fraction_values[q] +
                             bdf_coefs[1] * p1_void_fraction_values[q] +
                             bdf_coefs[2] * p2_void_fraction_values[q] +
                             bdf_coefs[3] * p3_void_fraction_values[q]) *
                            div_phi_u[i] * JxW;
                        }
                    }

                  if (velocity_source ==
                      Parameters::VelocitySource::VelocitySourceType::srf)
                    {
                      if (dim == 2)
                        {
                          local_rhs(i) +=
                            -2 * omega_z * (-1.) *
                            cross_product_2d(present_velocity_values[q]) *
                            phi_u[i] * JxW;
                          auto centrifugal =
                            omega_z * (-1.) *
                            cross_product_2d(
                              omega_z * (-1.) *
                              cross_product_2d(quadrature_points[q]));
                          local_rhs(i) += -centrifugal * phi_u[i] * JxW;
                        }
                      else if (dim == 3)
                        {
                          local_rhs(i) +=
                            -2 *
                            cross_product_3d(omega_vector,
                                             present_velocity_values[q]) *
                            phi_u[i] * JxW;
                          local_rhs(i) +=
                            -cross_product_3d(
                              omega_vector,
                              cross_product_3d(omega_vector,
                                               quadrature_points[q])) *
                            phi_u[i] * JxW;
                        }
                    }

                  // PSPG GLS term
                  if (PSPG)
                    local_rhs(i) +=
                      -tau * (strong_residual * grad_phi_p[i]) * JxW;

                  // SUPG GLS term
                  if (SUPG)
                    {
                      local_rhs(i) +=
                        -tau *
                        (strong_residual *
                         (grad_phi_u[i] * present_velocity_values[q])) *
                        JxW;
                    }

                  if (this->simulation_parameters.cfd_dem.shock_capturing ==
                      true)
                    local_rhs(i) +=
                      -vdcdd *
                      scalar_product(present_velocity_gradients[q],
                                     grad_phi_u[i]) *
                      JxW;

                  // Grad-div stabilization
                  if (this->simulation_parameters.cfd_dem.grad_div == true)
                    {
                      local_rhs(i) -= gamma *
                                      (present_void_fraction_values[q] *
                                         present_velocity_divergence +
                                       present_velocity_values[q] *
                                         present_void_fraction_gradients[q]) *
                                      div_phi_u[i] * JxW;
                    }
                }
            }

          cell->get_dof_indices(local_dof_indices);

          // The non-linear solver assumes that the nonzero
          // constraints have already been applied to the solution
          const AffineConstraints<double> &constraints_used =
            this->zero_constraints;
          // initial_step ? nonzero_constraints : zero_constraints;
          if (assemble_matrix)
            {
              constraints_used.distribute_local_to_global(local_matrix,
                                                          local_rhs,
                                                          local_dof_indices,
                                                          this->system_matrix,
                                                          system_rhs);
            }
          else
            {
              constraints_used.distribute_local_to_global(local_rhs,
                                                          local_dof_indices,
                                                          system_rhs);
            }
        }
    }
  if (assemble_matrix)
    this->system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}



template <int dim>
void
GLSVANSSolver<dim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_system");

  if (this->simulation_parameters.velocity_sources.type ==
      Parameters::VelocitySource::VelocitySourceType::none)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
    }

  else if (this->simulation_parameters.velocity_sources.type ==
           Parameters::VelocitySource::VelocitySourceType::srf)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
    }
}
template <int dim>
void
GLSVANSSolver<dim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_rhs");

  if (this->simulation_parameters.velocity_sources.type ==
      Parameters::VelocitySource::VelocitySourceType::none)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
    }
  if (this->simulation_parameters.velocity_sources.type ==
      Parameters::VelocitySource::VelocitySourceType::srf)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
    }
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
                void_fraction_m1, p1_void_fraction_values);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                fe_values_void_fraction.get_function_values(
                  void_fraction_m2, p2_void_fraction_values);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                fe_values_void_fraction.get_function_values(
                  void_fraction_m3, p3_void_fraction_values);
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
    this->simulation_parameters.cfd_dem.inlet_boundary_id,
    this->simulation_parameters.cfd_dem.outlet_boundary_id);

  this->pcout << "Mass Source: " << mass_source << " s^-1" << std::endl;
  this->pcout << "Average Void Fraction in Bed: " << average_void_fraction
              << std::endl;
  this->pcout << "Pressure Drop: " << pressure_drop << " m^2.s^-2" << std::endl;
}

template <int dim>
void
GLSVANSSolver<dim>::solve()
{
  read_mesh_and_manifolds(
    this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart ||
      this->simulation_parameters.void_fraction->read_dem == true,
    this->simulation_parameters.boundary_conditions);

  if (this->simulation_parameters.void_fraction->read_dem == true)
    read_dem();
  setup_dofs();
  calculate_void_fraction(this->simulation_control->get_current_time());
  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);

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

      if (this->simulation_parameters.cfd_dem.post_processing)
        post_processing();
    }

  this->finish_simulation();
}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the
// library is valid before we actually compile the solver This greatly
// helps with debugging
template class GLSVANSSolver<2>;
template class GLSVANSSolver<3>;
