#include <deal.II/base/convergence_table.h>

#include "glsNS.h"

template <int dim>
class ExactSolutionTaylorCouette : public Function<dim>
{
public:
  ExactSolutionTaylorCouette()
    : Function<dim>(3)
  {
    eta_ = 0.25;
    ri_  = 0.25;
  }
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const;

private:
  double eta_;
  double ri_ = 0.25;
};
template <int dim>
void
ExactSolutionTaylorCouette<dim>::vector_value(const Point<dim> &p,
                                              Vector<double> &  values) const
{
  double x      = p[0];
  double y      = p[1];
  double r      = std::sqrt(x * x + y * y);
  double theta  = std::atan2(y, x);
  double A      = -(eta_ * eta_) / (1. - eta_ * eta_);
  double B      = ri_ * ri_ / (1. - eta_ * eta_);
  double utheta = A * r + B / r;
  values(0)     = -std::sin(theta) * utheta;
  values(1)     = std::cos(theta) * utheta;
  values(2)     = 0.;
}

template <int dim>
class TaylorCouetteNavierStokes : public GLSNavierStokesSolver<dim>
{
public:
  TaylorCouetteNavierStokes(NavierStokesSolverParameters<dim> nsparam,
                            const unsigned int                degreeVelocity,
                            const unsigned int                degreePressure);
  void
  run();
};

template <int dim>
TaylorCouetteNavierStokes<dim>::TaylorCouetteNavierStokes(
  NavierStokesSolverParameters<dim> nsparam,
  const unsigned int                degreeVelocity,
  const unsigned int                degreePressure)
  : GLSNavierStokesSolver<dim>(nsparam, degreeVelocity, degreePressure)
{}

template <int dim>
void
TaylorCouetteNavierStokes<dim>::run()
{
  Point<dim> circleCenter;
  if (dim == 2)
    circleCenter = Point<dim>(0, 0);
  this->read_mesh();

  this->setup_dofs();
  this->forcing_function = new NoForce<dim>;
  this->exact_solution   = new ExactSolutionTaylorCouette<dim>;

  ConvergenceTable table;

  this->iterate(this->simulationControl.firstIter());
  while (this->simulationControl.integrate())
    {
      printTime(this->pcout, this->simulationControl);
      if (this->simulationControl.firstIter())
        {
          this->iterate(this->simulationControl.firstIter());
        }
      else
        {
          this->refine_mesh();
          this->iterate(this->simulationControl.firstIter());
        }
      this->postprocess();
      table.add_value("cells", this->triangulation.n_global_active_cells());

      const double error = this->calculate_L2_error();
      table.add_value("error", error);

      this->finish_time_step();
    }
  table.omit_column_from_convergence_rate_evaluation("cells");
  table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
  table.set_scientific("error", true);
  if (this->this_mpi_process == 0)
    {
      table.write_text(std::cout);
    }
}

int
main(int argc, char *argv[])
{
  try
    {
      if (argc != 2)
        {
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
          std::exit(1);
        }
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);

      ParameterHandler                prm;
      NavierStokesSolverParameters<2> NSparam;
      NSparam.declare(prm);
      // Parsing of the file
      prm.parse_input(argv[1]);
      NSparam.parse(prm);

      TaylorCouetteNavierStokes<2> problem_2d(
        NSparam,
        NSparam.femParameters.velocityOrder,
        NSparam.femParameters.pressureOrder);
      problem_2d.run();
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
