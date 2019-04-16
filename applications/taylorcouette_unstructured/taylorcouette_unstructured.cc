#include "glsNS.h"
#include <deal.II/base/convergence_table.h>

template<int dim>
class ExactSolutionTaylorCouette : public Function<dim>
{
public:
    ExactSolutionTaylorCouette() : Function<dim>(3)
    {
        eta_=0.25;
        ri_=0.25;
    }
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;

private:
    double eta_;
    double ri_=0.25;
};
template<int dim>
void ExactSolutionTaylorCouette<dim>::vector_value(const Point<dim> &p,
                                                    Vector<double> &values) const
{
    double x = p[0];
    double y = p[1];
    double r= std::sqrt(x*x+y*y);
    double theta= std::atan2(y,x);
    double A= -(eta_*eta_)/(1.-eta_*eta_);
    double B= ri_ * ri_ / (1.-eta_*eta_);
    double utheta= A*r + B/r;
    values(0) = -std::sin(theta)*utheta;
    values(1) = std::cos(theta)*utheta;
    values(2) = 0.;
}

template <int dim>
class TaylorCouetteNavierStokes : public GLSNavierStokesSolver<dim>
{
public:
  TaylorCouetteNavierStokes(const std::string input_filename, const unsigned int degreeVelocity, const unsigned int degreePressure);
  void run();

private:

  std::vector<double>          wallTime_;
};

template <int dim>
TaylorCouetteNavierStokes<dim>::TaylorCouetteNavierStokes(const std::string input_filename, const unsigned int degreeVelocity, const unsigned int degreePressure):
  GLSNavierStokesSolver<dim>(input_filename,degreeVelocity,degreePressure)
{

}


template<int dim>
void TaylorCouetteNavierStokes<dim>::run()
{
  Point<dim>           circleCenter;

  if (dim==2) circleCenter = Point<dim>(0,0);

  GridIn<dim> grid_in;
  grid_in.attach_triangulation (this->triangulation);
  std::ifstream input_file(this->meshParameters.fileName);
  grid_in.read_msh(input_file);

  static const SphericalManifold<dim> manifold_description(circleCenter);
  this->triangulation.set_manifold (0, manifold_description);
  this->triangulation.set_all_manifold_ids_on_boundary(0);

  this->triangulation.refine_global(this->meshParameters.initialRefinement);

  this->setup_dofs();
  this->forcing_function = new NoForce<dim>;
  this->exact_solution = new ExactSolutionTaylorCouette<dim>;

  ConvergenceTable table;

  while(this->simulationControl.integrate())
    {
      printTime(this->pcout,this->simulationControl);
      if (this->simulationControl.firstIter())
      {
        this->newton_iteration(true);
      }
      else
      {
        this->refine_mesh();
        this->newton_iteration(false);
      }
      this->postprocess();
      table.add_value("cells", this->triangulation.n_global_active_cells());

      const double error = this->calculateL2Error();
      table.add_value("error",   error);

      this->finishTimeStep();
    }
  table.omit_column_from_convergence_rate_evaluation("cells");
  table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
  table.set_scientific("error", true);
  if(this->this_mpi_process==0)
  {
    table.write_text(std::cout);
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
        TaylorCouetteNavierStokes<2> problem_2d(argv[1],fem.velocityOrder,fem.pressureOrder);
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
