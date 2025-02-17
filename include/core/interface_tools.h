// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_interface_tools_h
#define lethe_interface_tools_h

#include <deal.II/base/function.h>

#include <deal.II/fe/fe.h>

using namespace dealii;
namespace InterfaceTools
{
  /**
   * @brief Scalar function defined by the DOF values of a single cell. Based on the CellWiseFunction and RefSpaceFEFieldFunction of dealii.
   */
  template <int dim, typename VectorType = Vector<double>>
  class CellWiseFunction : public Function<dim>
  {
  public:
    /**
     * @brief Constructor.
     *
     * @param[in] p_fe Finite element discretizing the field we want to convert to
     * a CellWiseFunction.
     *
     */
    CellWiseFunction(const FiniteElement<dim, dim> &p_fe);
  
    /**
     * @brief Set the cell that the function should be evaluated on. 
     *
     * @param[in] in_local_dof_values Cell's DOF values
     *
     */
      void
    set_active_cell(const VectorType &in_local_dof_values);
  private:
    /// Finite element discretizing the field of interest
    FiniteElement<dim> fe;
    
    /// Number of dofs per element
    unsigned int n_cell_wise_dofs;
    
  };
  
  template <int dim, typename VectorType>
  CellWiseFunction<dim, VectorType>::CellWiseFunction(
    const FiniteElement<dim, dim> &p_fe)
    : fe(p_fe)
  {
    n_cell_wise_dofs = fe.dofs_per_cell;
  }


} // namespace InterfaceTools

#endif
