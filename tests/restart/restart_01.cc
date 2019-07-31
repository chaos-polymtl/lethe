// check the read and write of simulationcontrol

#include "simulationcontrol.h"
#include "parameters.h"
#include "navierstokessolverparameters.h"
#include "../tests.h"

#include "glsNS.h"

template <int dim>
class RestartNavierStokes : public GLSNavierStokesSolver<dim>
{
public:
  RestartNavierStokes(NavierStokesSolverParameters<dim> nsparam, const unsigned int degreeVelocity, const unsigned int degreePressure):
    GLSNavierStokesSolver<dim>(nsparam, degreeVelocity,degreePressure)
  {
  }
  void run();
};

template<int dim>
void RestartNavierStokes<dim>::run()
{
  const int initialSize=4;
  this->make_cube_grid(initialSize);
  this->setup_dofs();
  this->exact_solution = new ExactSolutionMMS<dim>;
  this->forcing_function = new MMSSineForcingFunction<dim>;
  this->viscosity_=1;

  printTime(this->pcout,this->simulationControl);
  this->iterate(false);
  this->postprocess();
  double error1 = this->calculate_L2_error();
  deallog << "Error 1 : " << error1 << std::endl;
  this->finish_time_step();

  this->set_solution_vector(0.);
  double error2 = this->calculate_L2_error();

  deallog << "Error 2 : " << error2 << std::endl;
  if(error2<error1) std::runtime_error("Zeroed solution has a lower error");

  printTime(this->pcout,this->simulationControl);
  this->triangulation.clear();
  this->make_cube_grid(0);

  this->set_initial_condition(this->nsparam.initialCondition->type,true);

  double error3 = this->calculate_L2_error();
  deallog << "Error 3 : " << error3 << std::endl;
  if(!approximatelyEqual(error1,error3,1e-10)) std::runtime_error("Reloaded solution is not the same");
}

int
main(int argc, char* argv[])
{
  try
  {
    initlog();
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);

    ParameterHandler prm;
    NavierStokesSolverParameters<2> NSparam;
    NSparam.declare(prm);
    NSparam.parse(prm);

    // Manually alter some of the default parameters of the solver
    NSparam.restartParameters.checkpoint=true;
    NSparam.restartParameters.frequency=1;
    NSparam.nonLinearSolver.verbosity=Parameters::NonLinearSolver::quiet;
    NSparam.linearSolver.verbosity=Parameters::LinearSolver::quiet;
    NSparam.boundaryConditions.createDefaultNoSlip();

    RestartNavierStokes<2> problem_2d(NSparam,NSparam.femParameters.velocityOrder,NSparam.femParameters.pressureOrder);
    problem_2d.run();
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what()  << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
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
