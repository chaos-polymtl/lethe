#include "glsNS.h"
#include <deal.II/base/convergence_table.h>


template<int dim>
class ConstantXForce : public Function<dim>
{
public:
    ConstantXForce() : Function<dim>(3) {};
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;

};
template<int dim>
void ConstantXForce<dim>::vector_value(const Point<dim> &/*p*/,
                                Vector<double> &values) const
{
    assert(dim==2);
    values(0) = 1.;
    values(1) = 0.;

}

template<int dim>
class ExactSolutionPoiseuille : public Function<dim>
{
public:
    ExactSolutionPoiseuille() : Function<dim>(3) {}
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;
};
template<int dim>
void ExactSolutionPoiseuille<dim>::vector_value(const Point<dim> &p,
                                         Vector<double> &values) const
{
    assert(dim==2);
    const double H=1.;
    const double G=1.;
    const double mu=1.;

    //double x = p[0];
    double y = p[1];
    values(0) = 0.5*G/mu * (y)*(H-y);
    values(1) = 0.;
}

template <int dim>
class PeriodicPoiseuille : public GLSNavierStokesSolver<dim>
{
public:
  PeriodicPoiseuille(NavierStokesSolverParameters<dim> nsparam, const unsigned int degreeVelocity, const unsigned int degreePressure):
    GLSNavierStokesSolver<dim>(nsparam, degreeVelocity,degreePressure){}
  void run();
};

template<int dim>
void PeriodicPoiseuille<dim>::run()
{
  this->read_mesh();
  this->setup_dofs();
  this->forcing_function = new ConstantXForce<dim>;
  this->exact_solution = new ExactSolutionPoiseuille<dim>;

  ConvergenceTable table;
  this->set_initial_condition(this->nsparam.initialCondition->type);
  while(this->simulationControl.integrate())
  {
    printTime(this->pcout,this->simulationControl);
    if (this->simulationControl.getIter() !=1) this->refine_mesh();
    this->iterate(this->simulationControl.firstIter());
    this->postprocess();
    table.add_value("cells", this->triangulation.n_global_active_cells());
    const double error = this->calculate_L2_error();
    table.add_value("error",   error);
    this->finish_time_step();
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
        ParameterHandler prm;
        NavierStokesSolverParameters<2> nsparam;
        nsparam.declare(prm);
        // Parsing of the file
        prm.parse_input (argv[1]);
        nsparam.parse(prm);

        PeriodicPoiseuille<2> problem_2d(nsparam,nsparam.femParameters.velocityOrder,nsparam.femParameters.pressureOrder);
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
