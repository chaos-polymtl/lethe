#include "glsNS.h"

template<int dim>
class ExactSolutionMMS3D : public Function<dim>
{
public:
    ExactSolutionMMS3D() : Function<dim>(3) {}
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;
};
template<int dim>
void ExactSolutionMMS3D<dim>::vector_value(const Point<dim> &p,
                                                    Vector<double> &values) const
{
    const double a = M_PI;
    double x = p[0];
    double y = p[1];
    const double z = p[2];
    values(0) = sin(a*x)*sin(a*x)*cos(a*y)*sin(a*y)*cos(a*z)*sin(a*z);
    values(1) = cos(a*x)*sin(a*x)*sin(a*y)*sin(a*y)*cos(a*z)*sin(a*z);
    values(2) = -2*cos(a*x)*sin(a*x)*cos(a*y)*sin(a*y)*sin(a*z)*sin(a*z);
}

template <int dim>
class MMS3DNavierStokes : public GLSNavierStokesSolver<dim>
{
public:
  MMS3DNavierStokes(const std::string input_filename, const unsigned int degreeVelocity, const unsigned int degreePressure);
  void runMMS_3D();

private:
  std::vector<double>          wallTime_;
};

template <int dim>
MMS3DNavierStokes<dim>::MMS3DNavierStokes(const std::string input_filename, const unsigned int degreeVelocity, const unsigned int degreePressure):
  GLSNavierStokesSolver<dim>(input_filename,degreeVelocity,degreePressure)
{

}

template<int dim>
void MMS3DNavierStokes<dim>::runMMS_3D()
{
    assert(dim==3);
    const int initialSize=this->meshParameters.initialRefinement;

    this->make_cube_grid(initialSize);
    this->setup_dofs();

    this->exact_solution = new ExactSolutionMMS3D<dim>;
    this->forcing_function = new MMS3DSineForcingFunction<dim>;
    this->viscosity_=this->physicalProperties.viscosity;

    Timer timer;
    while(this->simulationControl.integrate())
    {
      printTime(this->pcout,this->simulationControl);
      timer.start ();
      if (this->simulationControl.getIter() !=1)
      {
        this->refine_mesh();
      }
      this->setSolutionVector(0.);
      this->newton_iteration(false);
      this->postprocess();
      this->oldCalculateL2Error();

      this->wallTime_.push_back((timer.wall_time()));
      this->finishTimeStep();
    }

    if(this->this_mpi_process==0)
    {
      assert (this->wallTime_.size()==this->L2ErrorU_.size());
      std::ofstream output_file("./L2Error-3D.dat");
      for (unsigned int i=0 ; i < this->L2ErrorU_.size() ; ++i)
      {
        output_file << i+initialSize << " " << this->L2ErrorU_[i] << " " << this->wallTime_[i] << std::endl;
      }
      output_file.close();
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

        Parameters::FEM              fem;
        fem=Parameters::getFEMParameters2D(argv[1]);
        MMS3DNavierStokes<3> problem(argv[1],fem.velocityOrder,fem.pressureOrder);
        problem.runMMS_3D();
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
