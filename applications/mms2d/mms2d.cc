#include "glsNS.h"
#include "parameters.h"

template <int dim>
class MMSNavierStokes : public GLSNavierStokesSolver<dim>
{
public:
  MMSNavierStokes(NavierStokesSolverParameters<dim> nsparam, const unsigned int degreeVelocity, const unsigned int degreePressure):
    GLSNavierStokesSolver<dim>(nsparam, degreeVelocity,degreePressure)
  {
  }
  void run();
};

template<int dim>
void MMSNavierStokes<dim>::run()
{
  std::vector<double>                   ErrorLog;
  std::vector<double>                   wallTime;
  const int initialSize=this->meshParameters.initialRefinement;
  this->make_cube_grid(initialSize);
  this->setup_dofs();
  this->exact_solution = new ExactSolutionMMS<dim>;
  this->forcing_function = new MMSSineForcingFunction<dim>;
  this->viscosity_=this->physicalProperties.viscosity;

  Timer timer;
  this->set_initial_condition(this->nsparam.initialCondition->type, this->restartParameters.restart);
  while(this->simulationControl.integrate())
  {
    printTime(this->pcout,this->simulationControl);
    timer.start();
    if (this->simulationControl.getIter() !=1) this->refine_mesh();
    this->iterate(this->simulationControl.firstIter());
    this->postprocess();
    {
      double L2Error= this->calculate_L2_error();
      this->pcout << "L2Error U is : " << std::setprecision(this->analyticalSolutionParameters.errorPrecision) << L2Error << std::endl;
      ErrorLog.push_back(L2Error);
      wallTime.push_back((timer.wall_time()));
    }
    this->finish_time_step();
  }

  if(this->this_mpi_process==0)
  {
    assert (wallTime.size()==ErrorLog.size());
    std::ofstream output_file("./L2Error-2D.dat");
    for (unsigned int i=0 ; i < ErrorLog.size() ; ++i)
    {
      output_file << i+initialSize << " " << ErrorLog[i] << " " << wallTime[i] << std::endl;
    }
    output_file.close();
  }
}

int main (int argc, char *argv[])
{
    try
    {
    if (argc < 2)
      {
        std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
        std::exit(1);
      }
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);

        ParameterHandler prm;
        NavierStokesSolverParameters<2> NSparam;
        NSparam.declare(prm);
        prm.parse_input (argv[1]);
        NSparam.parse(prm);

        MMSNavierStokes<2> problem_2d(NSparam,NSparam.femParameters.velocityOrder,NSparam.femParameters.pressureOrder);
        problem_2d.run();
    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl << std::endl
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
