#include "glsNS.h"

template <int dim>
class VonKarmanNavierStokes : public GLSNavierStokesSolver<dim>
{
public:
  VonKarmanNavierStokes(NavierStokesSolverParameters<dim> nsparam, const unsigned int degreeVelocity, const unsigned int degreePressure):
    GLSNavierStokesSolver<dim>(nsparam, degreeVelocity,degreePressure){}
  void run();
  void runTest();
};

template<int dim>
void VonKarmanNavierStokes<dim>::run()
{
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (this->triangulation);
  std::ifstream input_file(this->meshParameters.fileName);
  grid_in.read_msh(input_file);


  Point<dim,double> circleCenter(8,8);
  static const SphericalManifold<dim> boundary(circleCenter);
  this->triangulation.set_all_manifold_ids_on_boundary(0,0);
  this->triangulation.set_manifold (0, boundary);
  this->setup_dofs();
  this->forcing_function = new NoForce<dim>;

  this->setInitialCondition(this->initialConditionParameters->type);
  while(this->simulationControl.integrate())
    {
      printTime(this->pcout,this->simulationControl);
      this->newton_iteration(false);
      this->postprocess();
      this->refine_mesh();
      this->finishTimeStep();
    }
}

template<int dim>
void VonKarmanNavierStokes<dim>::runTest()
{
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (this->triangulation);
  std::ifstream input_file(this->meshParameters.fileName);
  grid_in.read_msh(input_file);


  Point<dim,double> circleCenter(8,8);
  static const SphericalManifold<dim> boundary(circleCenter);
  this->triangulation.set_all_manifold_ids_on_boundary(0,0);
  this->triangulation.set_manifold (0, boundary);
  this->setup_dofs();
  this->forcing_function = new NoForce<dim>;

  this->setInitialCondition(this->initialConditionParameters->type);
  while(this->simulationControl.integrate())
    {
      printTime(this->pcout,this->simulationControl);
      this->newton_iteration(false);
      this->postprocess();
      this->refine_mesh();
      this->finishTimeStep();
      for (unsigned int i =0 ; i < this->forces_.size() ; ++i)
      {
        this->pcout << " fx : " << this->forces_[i][0] << std::endl;
      }
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
        NavierStokesSolverParameters<2> nsparam;
        nsparam.declare(prm);
        // Parsing of the file
        prm.parse_input (argv[1]);
        nsparam.parse(prm);

        VonKarmanNavierStokes<2> problem_2d(nsparam,nsparam.femParameters.velocityOrder,nsparam.femParameters.pressureOrder);
        if (nsparam.test.enabled) problem_2d.runTest();
        else problem_2d.run();
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
