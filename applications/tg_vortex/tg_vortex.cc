#include "glsNS.h"
#include "parameters.h"

template <int dim> class ExactSolutionTGV : public Function<dim> {
public:
  ExactSolutionTGV(double p_viscosity, double p_time)
      : Function<dim>(dim + 1), viscosity(p_viscosity), time(p_time) {}
  virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;

private:
  double viscosity;
  double time;
};

template <int dim>
void ExactSolutionTGV<dim>::vector_value(const Point<dim> &p,
                                         Vector<double> &values) const {
  assert(dim == 2);
  double x = p[0];
  double y = p[1];
  double factor = std::exp(-2. * viscosity * time);
  values(0) = cos(x) * sin(y) * factor;
  values(1) = -sin(x) * cos(y) * factor;
}

template <int dim> class TaylorGreenVortex : public GLSNavierStokesSolver<dim> {
public:
  TaylorGreenVortex(NavierStokesSolverParameters<dim> nsparam,
                    const unsigned int degreeVelocity,
                    const unsigned int degreePressure)
      : GLSNavierStokesSolver<dim>(nsparam, degreeVelocity, degreePressure) {}
  void run2DTGV();
};

template <int dim> void TaylorGreenVortex<dim>::run2DTGV() {
  std::vector<double> ke_values;
  std::vector<double> timeTaken;
  std::vector<double> enstrophy_values;
  this->read_mesh();
  this->setup_dofs();
  this->forcing_function = new NoForce<dim>;
  this->exact_solution =
      new ExactSolutionTGV<dim>(this->nsparam.physicalProperties.viscosity, 0.);
  this->set_initial_condition(this->nsparam.initialCondition->type,
                              this->nsparam.restartParameters.restart);

  Timer timer;
  while (this->simulationControl.integrate()) {
    printTime(this->pcout, this->simulationControl);
    timer.start();
    this->refine_mesh();
    this->iterate(this->simulationControl.firstIter());
    this->postprocess();

    {
      delete this->exact_solution;
      this->exact_solution =
          new ExactSolutionTGV<dim>(this->nsparam.physicalProperties.viscosity,
                                    this->simulationControl.getTime());
      double error = this->calculate_L2_error() /
                     std::exp(-2 * this->nsparam.physicalProperties.viscosity *
                              this->simulationControl.getTime());
      this->pcout << "L2 error : "
                  << std::setprecision(
                         this->nsparam.analyticalSolution.errorPrecision)
                  << error << std::endl;
      double kE = this->calculate_average_KE();
      this->pcout << "Kinetic energy : "
                  << std::setprecision(
                         this->nsparam.analyticalSolution.errorPrecision)
                  << kE << std::endl;
      ke_values.push_back(kE);
      timeTaken.push_back((this->simulationControl.getTime()));
      double enstrophy = this->calculate_average_enstrophy();
      this->pcout << "Enstrophy  : "
                  << std::setprecision(
                         this->nsparam.analyticalSolution.errorPrecision)
                  << enstrophy << std::endl;
      enstrophy_values.push_back(enstrophy);
    }

    this->finish_time_step();
  }
  if (this->this_mpi_process == 0) {
    assert(timeTaken.size() == ke_values.size());
    std::ofstream output_file("./KE-2D.dat");
    for (unsigned int i = 0; i < ke_values.size(); ++i) {

      output_file << timeTaken[i] << " " << ke_values[i] << std::endl;
    }
    output_file.close();
  }
  if (this->this_mpi_process == 0) {
    assert(timeTaken.size() == enstrophy_values.size());
    std::ofstream output_file("./Enstrophy-2D.dat");
    for (unsigned int i = 0; i < enstrophy_values.size(); ++i) {

      output_file << timeTaken[i] << " " << enstrophy_values[i] << std::endl;
    }
    output_file.close();
  }
}

int main(int argc, char *argv[]) {
  try {
    if (argc < 2) {
      std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
      std::exit(1);
    }
    Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);

    ParameterHandler prm;
    NavierStokesSolverParameters<2> NSparam;
    NSparam.declare(prm);
    prm.parse_input(argv[1]);
    NSparam.parse(prm);

    TaylorGreenVortex<2> problem_2d(NSparam,
                                    NSparam.femParameters.velocityOrder,
                                    NSparam.femParameters.pressureOrder);
    problem_2d.run2DTGV();
  } catch (std::exception &exc) {
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
  } catch (...) {
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
