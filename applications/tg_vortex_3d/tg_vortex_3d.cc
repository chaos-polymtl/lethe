#include "glsNS.h"
#include "parameters.h"

template <int dim>
class TaylorGreenVortex : public GLSNavierStokesSolver<dim>
{
public:
  TaylorGreenVortex(NavierStokesSolverParameters<dim> nsparam,
                    const unsigned int                degreeVelocity,
                    const unsigned int                degreePressure)
    : GLSNavierStokesSolver<dim>(nsparam, degreeVelocity, degreePressure)
  {
    if (this->nsparam.restartParameters.restart && this->this_mpi_process == 0)
      {
        read_ke();
        read_enstrophy();
        if (enstrophy_values.size() != ke_values.size())
          {
            throw std::runtime_error("Unalligned KE and enstrophy files");
          }
      }
  }
  void
  run3DTGV();

private:
  std::vector<double> time;
  std::vector<double> ke_values;
  std::vector<double> enstrophy_values;

  // Read previously stored values of enstrophy and Ke
  void
  read_enstrophy();
  void
  read_ke();
};

template <int dim>
void
TaylorGreenVortex<dim>::read_ke()
{
  std::ifstream input("./KE.dat");
  if (!input)
    {
      throw("Unable to open file");
    }
  std::string buffer;
  std::getline(input, buffer);

  std::string lineData;
  while (getline(input, lineData))
    {
      double              t, ke;
      std::vector<double> row;
      std::stringstream   lineStream(lineData);
      lineStream >> t;
      lineStream >> ke;
      time.push_back(t);
      ke_values.push_back(ke);
    }
}

template <int dim>
void
TaylorGreenVortex<dim>::read_enstrophy()
{
  std::ifstream input("./enstrophy.dat");
  if (!input)
    {
      throw("Unable to open file");
    }
  std::string buffer;
  std::getline(input, buffer);

  std::string lineData;
  while (getline(input, lineData))
    {
      double              t, ke;
      std::vector<double> row;
      std::stringstream   lineStream(lineData);
      lineStream >> t;
      lineStream >> ke;
      enstrophy_values.push_back(ke);
    }
}

// 3D TGV case
template <int dim>
void
TaylorGreenVortex<dim>::run3DTGV()
{
  const int initialSize = this->nsparam.mesh.initialRefinement;
  GridGenerator::hyper_cube(this->triangulation, 0, 2. * M_PI, true);
  this->set_periodicity();
  this->triangulation.refine_global(initialSize);
  this->setup_dofs();
  this->forcing_function = new NoForce<dim>;
  this->set_initial_condition(this->nsparam.initialCondition->type,
                              this->nsparam.restartParameters.restart);

  double kE        = this->calculate_average_KE();
  double enstrophy = this->calculate_average_enstrophy();
  if (this->this_mpi_process == 0 && !this->nsparam.restartParameters.restart)
    {
      time.push_back(0.);
      ke_values.push_back(kE);
      enstrophy_values.push_back(enstrophy);
    }

  Timer timer;
  while (this->simulationControl.integrate())
    {
      printTime(this->pcout, this->simulationControl);
      timer.start();
      this->refine_mesh();
      this->iterate(this->simulationControl.firstIter());
      this->postprocess(false);

      // Post-processing for Kinetic enegery and Enstrophy
      {
        double kE        = this->calculate_average_KE();
        double enstrophy = this->calculate_average_enstrophy();
        if (this->this_mpi_process == 0)
          {
            time.push_back((this->simulationControl.getTime()));
            ke_values.push_back(kE);
            enstrophy_values.push_back(enstrophy);

            this->pcout << "Kinetic energy : " << kE << std::endl;
            {
              std::ofstream output_file("./KE.dat");
              output_file << "Time Kinetic_Energy" << std::endl;
              for (unsigned int i = 0; i < ke_values.size(); ++i)
                output_file << std::setprecision(12) << time[i] << " "
                            << std::setprecision(12) << ke_values[i]
                            << std::endl;
              output_file.close();
            }

            this->pcout << "Enstrophy  : " << enstrophy << std::endl;
            {
              std::ofstream output_file("./enstrophy.dat");
              output_file << "Time Enstrophy" << std::endl;
              for (unsigned int i = 0; i < enstrophy_values.size(); ++i)
                output_file << std::setprecision(12) << time[i] << " "
                            << std::setprecision(12) << enstrophy_values[i]
                            << std::endl;
              output_file.close();
            }
          }
      }
      this->finish_time_step();
    }
}
int
main(int argc, char *argv[])
{
  try
    {
      if (argc < 2)
        {
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
          std::exit(1);
        }
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
      ParameterHandler                prm;
      NavierStokesSolverParameters<3> NSparam;
      NSparam.declare(prm);
      prm.parse_input(argv[1]);
      NSparam.parse(prm);
      TaylorGreenVortex<3> problem_3d(NSparam,
                                      NSparam.femParameters.velocityOrder,
                                      NSparam.femParameters.pressureOrder);
      problem_3d.run3DTGV();
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
