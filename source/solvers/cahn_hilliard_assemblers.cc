#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/cahn_hilliard_assemblers.h>
#include <solvers/copy_data.h>


template <int dim>
void
CahnHilliardAssemblerCore<dim>::assemble_matrix()
{
  return;
} // end loop on quadrature points



template <int dim>
void
CahnHilliardAssemblerCore<dim>::assemble_rhs()
{
  return;
}

template class CahnHilliardAssemblerCore<2>;
template class CahnHilliardAssemblerCore<3>;

template <int dim>
void
CahnHilliardAssemblerBDF<dim>::assemble_matrix()
{
  return;
}

template <int dim>
void
CahnHilliardAssemblerBDF<dim>::assemble_rhs()
{
  return;
}

template class CahnHilliardAssemblerBDF<2>;
template class CahnHilliardAssemblerBDF<3>;
