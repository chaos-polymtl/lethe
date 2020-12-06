/**
 * @brief This code checks the read and write of simulationcontrol.
 */

// Deal.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/types.h>

// Lethe
#include <core/parameters.h>
#include <solvers/gls_navier_stokes.h>
#include <solvers/navier_stokes_solver_parameters.h>

// Tests
#include <../tests/tests.h>

template <int dim>
class ExactSolutionMMS : public Function<dim>
{
public:
  ExactSolutionMMS()
    : Function<dim>(3)
  {}
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const;
};
template <int dim>
void
ExactSolutionMMS<dim>::vector_value(const Point<dim> &p,
                                    Vector<double> &  values) const
{
  assert(dim == 2);
  const double a = M_PI;
  double       x = p[0];
  double       y = p[1];
  values(0)      = sin(a * x) * sin(a * x) * cos(a * y) * sin(a * y);
  values(1)      = -cos(a * x) * sin(a * x) * sin(a * y) * sin(a * y);
}


template <int dim>
class MMSSineForcingFunction : public Function<dim>
{
public:
  MMSSineForcingFunction()
    : Function<dim>(3)
  {}
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const;
};
template <int dim>
void
MMSSineForcingFunction<dim>::vector_value(const Point<dim> &p,
                                          Vector<double> &  values) const
{
  assert(dim == 2);
  const double a = M_PI;
  const double x = p[0];
  const double y = p[1];
  values(0) =
    (2 * a * a * (-sin(a * x) * sin(a * x) + cos(a * x) * (cos(a * x))) *
       sin(a * y) * cos(a * y) -
     4 * a * a * sin(a * x) * sin(a * x) * sin(a * y) * cos(a * y) - 2.0 * x) *
      (-1.) +
    a * std::pow(sin(a * x), 3.) * std::pow(sin(a * y), 2.) * std::cos(a * x);
  values(1) =
    (2 * a * a * (sin(a * y) * (sin(a * y)) - cos(a * y) * cos(a * y)) *
       sin(a * x) * cos(a * x) +
     4 * a * a * sin(a * x) * sin(a * y) * sin(a * y) * cos(a * x) - 2.0 * y) *
      (-1) +
    a * std::pow(sin(a * x), 2.) * std::pow(sin(a * y), 3.) * std::cos(a * y);
}


template <int dim>
class RestartNavierStokes : public GLSNavierStokesSolver<dim>
{
public:
  RestartNavierStokes(NavierStokesSolverParameters<dim> nsparam)
    : GLSNavierStokesSolver<dim>(nsparam)
  {}
  void
  run();
};

template <int dim>
void
RestartNavierStokes<dim>::run()
{
  const int initialSize = 4;
  GridGenerator::hyper_cube(*this->triangulation, -1, 1);
  this->triangulation->refine_global(initialSize);
  this->setup_dofs();
  this->exact_solution                        = new ExactSolutionMMS<dim>;
  this->forcing_function                      = new MMSSineForcingFunction<dim>;
  this->nsparam.physical_properties.viscosity = 1.;

  this->simulation_control->print_progression(this->pcout);
  this->first_iteration();
  this->postprocess(false);
  auto & present_solution = this->get_present_solution();
  auto   errors_p1        = calculate_L2_error(this->dof_handler,
                                      present_solution,
                                      this->exact_solution,
                                      this->nsparam.fem_parameters,
                                      this->mpi_communicator);
  double error1           = errors_p1.first;
  deallog << "Error after first simulation : " << error1 << std::endl;
  this->finish_time_step();

  this->set_solution_vector(0.);
  auto errors_p2 = calculate_L2_error(this->dof_handler,
                                      present_solution,
                                      this->exact_solution,
                                      this->nsparam.fem_parameters,
                                      this->mpi_communicator);

  double error2 = errors_p2.first;

  deallog << "Error after zeroing the solution: " << error2 << std::endl;

  this->simulation_control->print_progression(this->pcout);
  this->triangulation->clear();
  GridGenerator::hyper_cube(*this->triangulation, -1, 1);
  this->triangulation->refine_global(0);

  this->set_initial_condition(this->nsparam.initial_condition->type, true);
  auto errors_p3 = calculate_L2_error(this->dof_handler,
                                      present_solution,
                                      this->exact_solution,
                                      this->nsparam.fem_parameters,
                                      this->mpi_communicator);

  double error3 = errors_p3.first;
  deallog << "Error after restarting the simulation: " << error3 << std::endl;
}

void
test()
{
  ParameterHandler                prm;
  NavierStokesSolverParameters<2> NSparam;
  NSparam.declare(prm);
  NSparam.parse(prm);

  // Manually alter some of the default parameters of the solver
  NSparam.restart_parameters.checkpoint = true;
  NSparam.restart_parameters.frequency  = 1;
  NSparam.non_linear_solver.verbosity   = Parameters::Verbosity::quiet;
  NSparam.linear_solver.verbosity       = Parameters::Verbosity::quiet;
  NSparam.boundary_conditions.createDefaultNoSlip();

  RestartNavierStokes<2> problem_2d(NSparam);
  problem_2d.run();
}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
      test();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
