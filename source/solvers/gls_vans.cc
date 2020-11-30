#include "solvers/gls_vans.h"



// Constructor for class GLS_VANS
template <int dim>
GLSVANSSolver<dim>::GLSVANSSolver(NavierStokesSolverParameters<dim> &p_nsparam,
                                  const unsigned int p_degree_velocity,
                                  const unsigned int p_degree_pressure)
  : GLSNavierStokesSolver<dim>(p_nsparam, p_degree_velocity, p_degree_pressure)
  , void_fraction_dof_handler(*this->triangulation)
  , fe_void_fraction(p_degree_velocity)

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

  const IndexSet locally_owned_dofs_voidfraction =
    void_fraction_dof_handler.locally_owned_dofs();
  IndexSet locally_relevant_dofs_voidfraction;
  DoFTools::extract_locally_relevant_dofs(void_fraction_dof_handler,

                                          locally_relevant_dofs_voidfraction);

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
}

template <int dim>
void
GLSVANSSolver<dim>::finish_time_step()
{
  GLSNavierStokesSolver<dim>::finish_time_step();

  void_fraction_m3 = void_fraction_m2;
  void_fraction_m2 = void_fraction_m1;
  void_fraction_m1 = nodal_void_fraction_relevant;
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
  const MappingQ<dim> mapping(this->velocity_fem_degree);

  this->nsparam.void_fraction->void_fraction.set_time(time);


  VectorTools::interpolate(mapping,
                           void_fraction_dof_handler,
                           this->nsparam.void_fraction->void_fraction,
                           nodal_void_fraction_owned);

  nodal_void_fraction_relevant = nodal_void_fraction_owned;
}

// Do an iteration with the NavierStokes Solver
// Handles the fact that we may or may not be at a first
// iteration with the solver and sets the initial condition
template <int dim>
void
GLSVANSSolver<dim>::first_iteration()
{
  // First step if the method is not a multi-step method
  if (!is_bdf_high_order(this->nsparam.simulation_control.method))
    {
      iterate();
    }

  // Taking care of the multi-step methods
  else if (this->nsparam.simulation_control.method ==
           Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    {
      Parameters::SimulationControl timeParameters =
        this->nsparam.simulation_control;

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
      this->solution_m2 = this->solution_m1;
      this->solution_m1 = this->present_solution;
      void_fraction_m2  = void_fraction_m1;
      void_fraction_m1  = nodal_void_fraction_relevant;

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

  else if (this->nsparam.simulation_control.method ==
           Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    {
      Parameters::SimulationControl timeParameters =
        this->nsparam.simulation_control;

      // Start the BDF3 with a single Euler time step with a lower time step
      double time_step =
        timeParameters.dt * timeParameters.startup_timestep_scaling;

      this->simulation_control->set_current_time_step(time_step);

      double intermediate_time = time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);

      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf1, false, true);
      this->solution_m2 = this->solution_m1;
      this->solution_m1 = this->present_solution;
      void_fraction_m2  = void_fraction_m1;
      void_fraction_m1  = nodal_void_fraction_relevant;

      // Reset the time step and do a bdf 2 newton iteration using the two
      // steps

      this->simulation_control->set_current_time_step(time_step);
      intermediate_time += time_step;
      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);

      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf1, false, true);
      this->solution_m3 = this->solution_m2;
      this->solution_m2 = this->solution_m1;
      this->solution_m1 = this->present_solution;
      void_fraction_m3  = void_fraction_m2;
      void_fraction_m2  = void_fraction_m1;
      void_fraction_m1  = nodal_void_fraction_relevant;

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
// iteration with the solver and sets the initial condition
template <int dim>
void
GLSVANSSolver<dim>::iterate()
{
  if (this->nsparam.simulation_control.method ==
      Parameters::SimulationControl::TimeSteppingMethod::sdirk2)
    {
      const double previous_time =
        this->simulation_control->get_previous_time();

      double stage_step =
        this->simulation_control->get_time_step() *
        time_fraction_at_sdirk_stage(
          Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1);

      double intermediate_time = previous_time + stage_step;

      this->forcing_function->set_time(intermediate_time);
      calculate_void_fraction(intermediate_time);

      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
        false,
        false);
      this->solution_m2 = this->present_solution;
      void_fraction_m2  = nodal_void_fraction_relevant;

      this->forcing_function->set_time(
        this->simulation_control->get_current_time());
      calculate_void_fraction(this->simulation_control->get_current_time());
      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
        false,
        false);
    }

  else if (this->nsparam.simulation_control.method ==
           Parameters::SimulationControl::TimeSteppingMethod::sdirk3)
    {
      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
        false,
        false);

      this->solution_m2 = this->present_solution;
      void_fraction_m2  = nodal_void_fraction_relevant;
      // BB TO FIX THIS IS WRONG
      calculate_void_fraction(this->simulation_control->get_current_time());

      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
        false,
        false);

      this->solution_m3 = this->present_solution;
      void_fraction_m3  = nodal_void_fraction_relevant;
      // BB TO FIX THIS IS WRONG
      calculate_void_fraction(this->simulation_control->get_current_time());

      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
        false,
        false);
    }
  else
    {
      calculate_void_fraction(this->simulation_control->get_current_time());
      this->forcing_function->set_time(
        this->simulation_control->get_current_time());
      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        this->nsparam.simulation_control.method, false, false);
    }
}


template <int dim>
template <bool                                              assemble_matrix,
          Parameters::SimulationControl::TimeSteppingMethod scheme,
          Parameters::VelocitySource::VelocitySourceType    velocity_source>
void
GLSVANSSolver<dim>::assembleGLS()
{
  auto &system_rhs = this->get_system_rhs();

  if (assemble_matrix)
    this->system_matrix = 0;
  this->system_rhs = 0;

  double         viscosity = this->nsparam.physical_properties.viscosity;
  Function<dim> *l_forcing_function = this->forcing_function;

  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(this->velocity_fem_degree,
                              this->nsparam.fem_parameters.qmapping_all);
  FEValues<dim>       fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients |
                            update_hessians);

  FEValues<dim> fe_values_void_fraction(mapping,
                                        this->fe_void_fraction,
                                        quadrature_formula,
                                        update_values |
                                          update_quadrature_points |
                                          update_JxW_values | update_gradients);

  const unsigned int               dofs_per_cell = this->fe.dofs_per_cell;
  const unsigned int               n_q_points    = quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  FullMatrix<double>               local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>                   local_rhs(dofs_per_cell);
  std::vector<Vector<double>> rhs_force(n_q_points, Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<Tensor<1, dim>>          present_velocity_values(n_q_points);
  std::vector<Tensor<2, dim>>          present_velocity_gradients(n_q_points);
  std::vector<double>                  present_pressure_values(n_q_points);
  std::vector<Tensor<1, dim>>          present_pressure_gradients(n_q_points);
  std::vector<Tensor<1, dim>>          present_velocity_laplacians(n_q_points);
  std::vector<Tensor<2, dim>>          present_velocity_hess(n_q_points);

  // Data storage vector for the void fraction values/gradients
  std::vector<double>         present_void_fraction_values(n_q_points);
  std::vector<Tensor<1, dim>> present_void_fraction_gradients(n_q_points);

  Tensor<1, dim> force;


  // Velocity dependent source term
  //----------------------------------
  // Angular velocity of the rotating frame. This is always a 3D vector even
  // in 2D.
  Tensor<1, dim> omega_vector;

  double omega_z  = this->nsparam.velocitySource.omega_z;
  omega_vector[0] = this->nsparam.velocitySource.omega_x;
  omega_vector[1] = this->nsparam.velocitySource.omega_y;
  if (dim == 3)
    omega_vector[2] = this->nsparam.velocitySource.omega_z;

  std::vector<double>         div_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<Tensor<3, dim>> hess_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>> laplacian_phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
  std::vector<double>         phi_p(dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi_p(dofs_per_cell);

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

  // Matrix of coefficients for the SDIRK methods
  // The lines store the information required for each step
  // Column 0 always refer to outcome of the step that is being calculated
  // Column 1 always refer to step n
  // Column 2+ refer to intermediary steps
  FullMatrix<double> sdirk_coefs;
  if (is_sdirk2(scheme))
    sdirk_coefs = sdirk_coefficients(2, dt);

  if (is_sdirk3(scheme))
    sdirk_coefs = sdirk_coefficients(3, dt);

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

          local_matrix = 0;
          local_rhs    = 0;

          // Gather velocity (values, gradient and laplacian)
          auto &evaluation_point = this->get_evaluation_point();
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

          // Gather the previous time steps depending on the number of stages
          // of the time integration scheme
          if (scheme !=
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            fe_values[velocities].get_function_values(this->solution_m1,
                                                      p1_velocity_values);

          if (time_stepping_method_has_two_stages(scheme))
            fe_values[velocities].get_function_values(this->solution_m2,
                                                      p2_velocity_values);

          if (time_stepping_method_has_three_stages(scheme))
            fe_values[velocities].get_function_values(this->solution_m3,
                                                      p3_velocity_values);

          // Gather the previous time steps depending on the number of stages
          // of the time integration scheme for the void fraction

          if (scheme !=
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            {
              fe_values_void_fraction.get_function_values(
                void_fraction_m1, p1_void_fraction_values);

              if (time_stepping_method_has_two_stages(scheme))
                fe_values_void_fraction.get_function_values(
                  void_fraction_m2, p2_void_fraction_values);

              if (time_stepping_method_has_three_stages(scheme))
                fe_values_void_fraction.get_function_values(
                  void_fraction_m3, p3_void_fraction_values);
            }



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


              // Establish the force vector
              for (int i = 0; i < dim; ++i)
                {
                  const unsigned int component_i =
                    this->fe.system_to_component_index(i).first;
                  force[i] = rhs_force[q](component_i);
                }
              const unsigned int component_mass =
                this->fe.system_to_component_index(dim).first;
              double mass_source = rhs_force[q](component_mass);

              // Calculate the divergence of the velocity
              const double present_velocity_divergence =
                trace(present_velocity_gradients[q]);

              // Calculate the strong residual for GLS stabilization
              auto strong_residual =
                present_velocity_gradients[q] * present_velocity_values[q] *
                  present_void_fraction_values[q]
                // Mass source term
                + mass_source * present_velocity_values[q]
                //..........................................
                + present_pressure_gradients[q] -
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

              /* Adjust the strong residual in cases where the scheme is
               transient.
               The BDF schemes require values at previous time steps which are
               stored in the p1, p2 and p3 vectors. The SDIRK scheme require
               the values at the different stages, which are also stored in
               the same arrays.
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


              if (is_sdirk_step1(scheme))
                strong_residual +=
                  (sdirk_coefs[0][0] * present_velocity_values[q] +
                   sdirk_coefs[0][1] * p1_velocity_values[q]) *
                  present_void_fraction_values[q];

              if (is_sdirk_step2(scheme))
                {
                  strong_residual +=
                    (sdirk_coefs[1][0] * present_velocity_values[q] +
                     sdirk_coefs[1][1] * p1_velocity_values[q] +
                     sdirk_coefs[1][2] * p2_velocity_values[q]) *
                    present_void_fraction_values[q];
                }

              if (is_sdirk_step3(scheme))
                {
                  strong_residual +=
                    (sdirk_coefs[2][0] * present_velocity_values[q] +
                     sdirk_coefs[2][1] * p1_velocity_values[q] +
                     sdirk_coefs[2][2] * p2_velocity_values[q] +
                     sdirk_coefs[2][3] * p3_velocity_values[q]) *
                    present_void_fraction_values[q];
                }

              // Matrix assembly
              if (assemble_matrix)
                {
                  // We loop over the column first to prevent recalculation of
                  // the strong jacobian in the inner loop
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      auto strong_jac =
                        (present_velocity_gradients[q] * phi_u[j] *
                           present_void_fraction_values[q] +
                         grad_phi_u[j] * present_velocity_values[q] *
                           present_void_fraction_values[q]
                         // Mass source term
                         + mass_source * phi_u[j]
                         //..........................................
                         + grad_phi_p[j] - viscosity * laplacian_phi_u[j]);

                      if (is_bdf(scheme))
                        strong_jac += present_void_fraction_values[q] *
                                      phi_u[j] * bdf_coefs[0];

                      if (is_sdirk(scheme))
                        strong_jac += present_void_fraction_values[q] *
                                      phi_u[j] * sdirk_coefs[0][0];

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

                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                          local_matrix(i, j) +=
                            (
                              // Momentum terms
                              viscosity *
                                scalar_product(grad_phi_u[j], grad_phi_u[i]) +
                              // TODO VERIFY ADVECTION TERMS
                              // Advection terms (eps.u.del(u))'
                              ((phi_u[j] * present_void_fraction_values[q] *
                                present_velocity_gradients[q] * phi_u[i]) +
                               (grad_phi_u[j] *
                                present_void_fraction_values[q] *
                                present_velocity_values[q] * phi_u[i]))
                              // Mass source term
                              + mass_source * phi_u[j] * phi_u[i]
                              //.........................................
                              // Pressure
                              - (div_phi_u[i] * phi_p[j]) +
                              // Continuity
                              phi_p[i] *
                                ((present_void_fraction_values[q] *
                                  div_phi_u[j]) +
                                 (phi_u[j] *
                                  present_void_fraction_gradients[q]))) *
                            JxW;

                          // Mass matrix
                          if (is_bdf(scheme))
                            local_matrix(i, j) +=
                              present_void_fraction_values[q] * phi_u[j] *
                              phi_u[i] * bdf_coefs[0] * JxW;

                          if (is_sdirk(scheme))
                            local_matrix(i, j) +=
                              present_void_fraction_values[q] * phi_u[j] *
                              phi_u[i] * sdirk_coefs[0][0] * JxW;

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
                      -viscosity * scalar_product(present_velocity_gradients[q],
                                                  grad_phi_u[i]) -
                      // Advection terms eps.udel(u)
                      (present_velocity_gradients[q] *
                       present_velocity_values[q] *
                       present_void_fraction_values[q] * phi_u[i])
                      // Mass source term
                      - mass_source * present_velocity_values[q] * phi_u[i]
                      //.............................................
                      // Pressure and force
                      + present_pressure_values[q] * div_phi_u[i] +
                      force * present_void_fraction_values[q] * phi_u[i] -
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
                    }

                  // std::cout << "mass source " << mass_source << std::endl;
                  // std::cout << "mass derivative "
                  //          << (bdf_coefs[0] * present_void_fraction_values[q]
                  //          +
                  //              bdf_coefs[1] * p1_void_fraction_values[q])
                  //          << std::endl;
                  // std::cout << "velocity divergence "
                  //          << present_velocity_divergence *
                  //               present_void_fraction_values[q]
                  //          << std::endl;
                  // std::cout << "void gradient "
                  //          << present_velocity_values[q] *
                  //               present_void_fraction_gradients[q]
                  //          << std::endl;


                  // std::cout << "first term : "
                  //          << (bdf_coefs[0] *
                  //          present_velocity_values[q] +
                  //              bdf_coefs[1] *
                  //              p1_velocity_values[q])
                  //          << std::endl;
                  // std::cout << "first term : "
                  //          << (bdf_coefs[0] *
                  //          present_void_fraction_values[q]
                  //          +
                  //              bdf_coefs[1] *
                  //              p1_void_fraction_values[q])
                  //          << std::endl;

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
                    }

                  // Residuals associated with SDIRK schemes
                  if (is_sdirk_step1(scheme))
                    {
                      local_rhs(i) -=
                        (sdirk_coefs[0][0] * present_velocity_values[q] +
                         sdirk_coefs[0][1] * p1_velocity_values[q]) *
                        present_void_fraction_values[q] * phi_u[i] * JxW;

                      local_rhs(i) -=
                        (sdirk_coefs[0][0] * present_void_fraction_values[q] +
                         sdirk_coefs[0][1] * p1_void_fraction_values[q]) *
                        phi_p[i] * JxW;
                    }

                  if (is_sdirk_step2(scheme))
                    {
                      local_rhs(i) -=
                        (sdirk_coefs[1][0] * present_velocity_values[q] +
                         sdirk_coefs[1][1] * p1_velocity_values[q] +
                         sdirk_coefs[1][2] * p2_velocity_values[q]) *
                        present_void_fraction_values[q] * phi_u[i] * JxW;

                      local_rhs(i) -=
                        (sdirk_coefs[1][0] * present_void_fraction_values[q] +
                         sdirk_coefs[1][1] * p1_void_fraction_values[q] +
                         sdirk_coefs[1][2] * p2_void_fraction_values[q]) *
                        phi_p[i] * JxW;
                    }

                  if (is_sdirk_step3(scheme))
                    {
                      local_rhs(i) -=
                        (sdirk_coefs[2][0] * present_velocity_values[q] +
                         sdirk_coefs[2][1] * p1_velocity_values[q] +
                         sdirk_coefs[2][2] * p2_velocity_values[q] +
                         sdirk_coefs[2][3] * p3_velocity_values[q]) *
                        present_void_fraction_values[q] * phi_u[i] * JxW;

                      local_rhs(i) -=
                        (sdirk_coefs[2][0] * present_void_fraction_values[q] +
                         sdirk_coefs[2][1] * p1_void_fraction_values[q] +
                         sdirk_coefs[2][2] * p2_void_fraction_values[q] +
                         sdirk_coefs[2][3] * p3_void_fraction_values[q]) *
                        phi_p[i] * JxW;
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
                }
            }

          cell->get_dof_indices(local_dof_indices);

          // The non-linear solver assumes that the nonzero constraints have
          // already been applied to the solution
          const AffineConstraints<double> &constraints_used =
            this->zero_constraints;
          // initial_step ? nonzero_constraints : zero_constraints;
          if (assemble_matrix)
            {
              constraints_used.distribute_local_to_global(local_matrix,
                                                          local_rhs,
                                                          local_dof_indices,
                                                          this->system_matrix,
                                                          this->system_rhs);
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
  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
GLSVANSSolver<dim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_system");

  if (this->nsparam.velocitySource.type ==
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
    }

  else if (this->nsparam.velocitySource.type ==
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
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

  if (this->nsparam.velocitySource.type ==
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
    }
  if (this->nsparam.velocitySource.type ==
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
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
                           nodal_void_fraction_owned,
                           "void_fraction");
}

template <int dim>
void
GLSVANSSolver<dim>::solve()
{
  read_mesh_and_manifolds(this->triangulation,
                          this->nsparam.mesh,
                          this->nsparam.manifolds_parameters,
                          this->nsparam.boundary_conditions);

  setup_dofs();
  calculate_void_fraction(this->simulation_control->get_current_time());
  this->set_initial_condition(this->nsparam.initial_condition->type,
                              this->nsparam.restart_parameters.restart);

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
    }

  this->finish_simulation();
}


// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library
// is valid before we actually compile the solver This greatly helps with
// debugging
template class GLSVANSSolver<2>;
