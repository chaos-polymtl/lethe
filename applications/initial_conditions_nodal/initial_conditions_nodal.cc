#include "glsNS.h"
#include "parameters.h"


template<int dim>
class ExactInitialSolution : public Function<dim>
{
public:
    ExactInitialSolution() : Function<dim>(3)
    {}
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;
};
template<int dim>
void ExactInitialSolution<dim>::vector_value(const Point<dim> &p,
                                                    Vector<double> &values) const
{
    double x = p[0];
    double y = p[1];
    values(0) = x;
    values(1) = y;
    values(2) = 0.;
}

template <int dim>
class InitialConditionsNavierStokes : public GLSNavierStokesSolver<dim>
{
public:
  InitialConditionsNavierStokes(NavierStokesSolverParameters<dim> nsparam, const unsigned int degreeVelocity, const unsigned int degreePressure):
    GLSNavierStokesSolver<dim>(nsparam, degreeVelocity,degreePressure)
  {
  }
  void run();
};



template<int dim>
void InitialConditionsNavierStokes<dim>::run()
{
  const int initialSize=this->meshParameters.initialRefinement;
  this->make_cube_grid(initialSize);
  this->setup_dofs();
  this->exact_solution = new ExactInitialSolution<dim>;
  this->setInitialCondition(Parameters::nodal);
  double error_nodal= this->calculateL2Error();
  if (error_nodal<1e-9)
  {
    this->pcout << "The nodal initial condition is OK"<< std::endl;
  }
  else
  {
    throw "Nodal initial condition error";
  }
}

int main (int argc, char *argv[])
{
    try
    {
    if (argc != 2)
      {
        std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
        std::exit(1);
      }
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);
        ParameterHandler prm;
        NavierStokesSolverParameters<2> NSparam;
        NSparam.declare(prm);
        // Parsing of the file
        prm.parse_input (argv[1]);
        NSparam.parse(prm);
        InitialConditionsNavierStokes<2> problem_2d(NSparam,NSparam.femParameters.velocityOrder,NSparam.femParameters.pressureOrder);
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
