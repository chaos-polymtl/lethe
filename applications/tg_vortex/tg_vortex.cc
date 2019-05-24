#include "glsNS.h"
#include "parameters.h"


template<int dim>
class ExactSolutionTGV : public Function<dim>
{
public:
    ExactSolutionTGV(double p_viscosity, double p_time) : Function<dim>(3),viscosity(p_viscosity),time(p_time) {}
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;
private:
    double viscosity;
    double time;
};

template<int dim>
void ExactSolutionTGV<dim>::vector_value(const Point<dim> &p,
                                         Vector<double> &values) const
{
    assert(dim==2);
    double x = p[0];
    double y = p[1];
    double factor = std::exp(-2.*viscosity*time);
    values(0) = cos(x)*sin(y)*factor;
    values(1) = -sin(x)*cos(y)*factor;
}

template <int dim>
class TaylorGreenVortex : public GLSNavierStokesSolver<dim>
{
public:
  TaylorGreenVortex(NavierStokesSolverParameters<dim> nsparam, const unsigned int degreeVelocity, const unsigned int degreePressure):
    GLSNavierStokesSolver<dim>(nsparam, degreeVelocity,degreePressure)
  {
  }
  void run2DTGV();
};

template<int dim>
void TaylorGreenVortex<dim>::run2DTGV()
{
  const int initialSize=this->meshParameters.initialRefinement;
  GridGenerator::hyper_cube (this->triangulation, 0, 2.*M_PI,true);
  this->setPeriodicity();
  this->triangulation.refine_global (initialSize);
  this->setup_dofs();
  this->forcing_function = new NoForce<dim>;
  this->viscosity_=this->physicalProperties.viscosity;
  this->exact_solution = new ExactSolutionTGV<dim>(this->viscosity_,0.);
  this->setInitialCondition(this->initialConditionParameters->type,this->restartParameters.restart);

  Timer timer;
  while(this->simulationControl.integrate())
  {
    printTime(this->pcout,this->simulationControl);
    timer.start ();
    this->refine_mesh();
    this->newton_iteration(false);
    this->postprocess();
    {
      delete this->exact_solution;
      this->exact_solution = new ExactSolutionTGV<dim>(this->viscosity_, this->simulationControl.getTime());
      double error = this->calculateL2Error() / std::exp(-2*this->viscosity_*this->simulationControl.getTime());
      this->pcout << "L2 error : " << std::setprecision(this->analyticalSolutionParameters.errorPrecision) << error << std::endl;
    }
    this->finishTimeStep();
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
        // Parsing of the file
        prm.parse_input (argv[1]);
        NSparam.parse(prm);

        if (NSparam.femParameters.dimension==2)
        {
          TaylorGreenVortex<2> problem_2d(NSparam,NSparam.femParameters.velocityOrder,NSparam.femParameters.pressureOrder);
          problem_2d.run2DTGV();
        }
        else
        {
          throw std::runtime_error("3D has not been implemented yet for this case");
        }
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
