#include "solvers/gls_vans.h"

// Constructor for class GLS_VANS
template <int dim>
GLSVANSSolver<dim>::GLSVANSSolver(NavierStokesSolverParameters<dim> &p_nsparam,
                                  const unsigned int p_degree_velocity,
                                  const unsigned int p_degree_pressure)
  : GLSNavierStokesSolver<dim>(p_nsparam, p_degree_velocity, p_degree_pressure)
  , void_fraction_dof_handler(*this->triangulation)
  , fe_void_fraction(p_degree_velocity)

{}


template <int dim>
void
GLSVANSSolver<dim>::setup_dofs()
{
  GLSNavierStokesSolver<dim>::setup_dofs();
  void_fraction_dof_handler.distribute_dofs(fe_void_fraction);

  const IndexSet locally_owned_dofs_voidfraction =
    void_fraction_dof_handler.locally_owned_dofs();
  IndexSet locally_relevant_dofs_voidfraction;
  DoFTools::extract_locally_relevant_dofs(void_fraction_dof_handler,
                                          locally_relevant_dofs_voidfraction);

  nodal_void_fraction_relevant.reinit(locally_owned_dofs_voidfraction,
                                      locally_relevant_dofs_voidfraction,
                                      this->mpi_communicator);
  nodal_void_fraction_owned.reinit(locally_owned_dofs_voidfraction,
                                   this->mpi_communicator);
}

template <int dim>
void
GLSVANSSolver<dim>::calculate_void_fraction()
{
  const MappingQ<dim> mapping(this->velocity_fem_degree);

  const double t = this->simulationControl->get_current_time();
  this->nsparam.void_fraction->void_fraction.set_time(t);

  VectorTools::interpolate(mapping,
                           void_fraction_dof_handler,
                           this->nsparam.void_fraction->void_fraction,
                           nodal_void_fraction_owned);

  nodal_void_fraction_relevant = nodal_void_fraction_owned;
}


template <int dim>
void
GLSVANSSolver<dim>::output_field_hook(DataOut<dim> &data_out)
{
  data_out.add_data_vector(void_fraction_dof_handler,
                           nodal_void_fraction_owned,
                           "void_fraction");
}

template <int dim>
void
GLSVANSSolver<dim>::solve()
{
  read_mesh_and_manifolds(this->triangulation,
                          this->nsparam.mesh,
                          this->nsparam.manifolds_parameters,
                          this->nsparam.boundary_conditions);

  setup_dofs();
  this->set_initial_condition(this->nsparam.initial_condition->type,
                              this->nsparam.restart_parameters.restart);

  while (this->simulationControl->integrate())
    {
      this->simulationControl->print_progression(this->pcout);
      if (this->simulationControl->is_at_start())
        this->first_iteration();
      else
        {
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          calculate_void_fraction();
          this->iterate();
        }
      this->postprocess(false);
      this->finish_time_step();
    }


  this->finish_simulation();
}


// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class GLSVANSSolver<2>;
