#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/copy_data.h>
#include <solvers/cahn_hilliard_assemblers.h>


template <int dim>
void
CahnHilliardAssemblerCore<dim>::assemble_matrix(CahnHilliardScratchData<dim> &scratch_data,
                                          StabilizedMethodsCopyData &copy_data)
{
  return;
} // end loop on quadrature points



template <int dim>
void
CahnHilliardAssemblerCore<dim>::assemble_rhs(CahnHilliardScratchData<dim> &   scratch_data,
                                       StabilizedMethodsCopyData &copy_data)
{
  return;
}

template class CahnHilliardAssemblerCore<2>;
template class CahnHilliardAssemblerCore<3>;

template <int dim>
void
CahnHilliardAssemblerBDF<dim>::assemble_matrix(CahnHilliardScratchData<dim> &scratch_data,
                                         StabilizedMethodsCopyData &copy_data)
{
  return;
}

template <int dim>
void
CahnHilliardAssemblerBDF<dim>::assemble_rhs(CahnHilliardScratchData<dim> &   scratch_data,
                                      StabilizedMethodsCopyData &copy_data)
{
  return;
}

template class CahnHilliardAssemblerBDF<2>;
template class CahnHilliardAssemblerBDF<3>;
