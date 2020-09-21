#include "solvers/gls_vans.h"
// Constructor for class GLS_VANS
template <int dim>
GLSVANSSolver<dim>::GLSVANSSolver(NavierStokesSolverParameters<dim> &p_nsparam,
                                  const unsigned int p_degree_velocity,
                                  const unsigned int p_degree_pressure)
  : GLSNavierStokesSolver<dim>(p_nsparam, p_degree_velocity, p_degree_pressure)

{}

template <int dim>
GLSVANSSolver<dim>::~GLSVANSSolver()
{}

template <int dim>
void
GLSVANSSolver<dim>::calculate_void_fraction()
{
    QGauss<2> quadrature_formula(this->fe.degree + 1);
    FEValues<2> fe_values(this->fe,
                          quadrature_formula,
                          update_values | update_gradients | update_JxW_values);

    const unsigned int  dofs_per_cell = this->fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const unsigned int n_dofs = this->dof_handler.n_dofs();
    v_fraction.reinit(dofs_per_cell, this->mpi_communicator);
    solution.reinit(n_dofs, this->mpi_communicator);

    for (const auto &cell : this->dof_handler.active_cell_iterators())
      {
        v_fraction = 0;

        for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
             {
               vertex = cell->vertex(i);
               //v_fraction = this->nsparam.void_fraction;
               v_fraction (i) = 0;
             }
        cell->get_dof_indices(local_dof_indices);

        for (const unsigned int i : fe_values.dof_indices())
            {
               solution(local_dof_indices[i]) = v_fraction(i);
            }
        }

}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class GLSVANSSolver<2>;
