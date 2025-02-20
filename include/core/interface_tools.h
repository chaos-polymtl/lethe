// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_interface_tools_h
#define lethe_interface_tools_h

#include <deal.II/base/function.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/non_matching/quadrature_generator.h>

using namespace dealii;

namespace InterfaceTools
{
  /**
   * @brief Scalar function defined by the DOF values of a single cell. Based on the CellWiseFunction and RefSpaceFEFieldFunction of dealii.
   */
  template <int dim,
            typename VectorType = Vector<double>,
            typename FEType     = FE_Q<dim>>
  class CellWiseFunction : public Function<dim>
  {
  public:
    /**
     * @brief Constructor.
     *
     * @param[in] p_fe Finite element discretizing the field we want to convert
     * to a CellWiseFunction.
     *
     */
    CellWiseFunction(const unsigned int p_fe_degree);

    /**
     * @brief Set the cell that the function should be evaluated on.
     *
     * @param[in] in_local_dof_values Cell's DOF values
     *
     */
    void
    set_active_cell(const VectorType &in_local_dof_values);

    /**
     * @brief Return the value of the function at the given point in the reference cell.
     *
     * @param[in] point Coordinates of the point in the reference cell
     *
     * @param[in] component Index of the component of the function for which the
     * value is computed
     *
     * @return Value of the function (or the component of interest) at the given point
     *
     */
    double
    value(const Point<dim>  &point,
          const unsigned int component = 0) const override;

    /**
     * @brief Return the function gradient of the specified component  at the given point in the reference cell.
     *
     * @param[in] point Coordinates of the point in the reference cell
     *
     * @param[in] component Index of the component of the function for which the
     * gradient is computed
     *
     * @return Value of the gradient of the specified component at the given point
     *
     */
    Tensor<1, dim>
    gradient(const Point<dim>  &point,
             const unsigned int component = 0) const override;

    /**
     * @brief Return the function Hessian of the specified component  at the given point in the reference cell.
     *
     * @param[in] point Coordinates of the point in the reference cell
     *
     * @param[in] component Index of the component of the function for which the
     * Hessian is computed
     *
     * @return Value of the Hessian of the specified component at the given point
     *
     */
    SymmetricTensor<2, dim>
    hessian(const Point<dim>  &point,
            const unsigned int component = 0) const override;

  private:
    /// Finite element discretizing the field of interest
    FEType fe;

    /// Number of dofs per element
    unsigned int n_cell_wise_dofs;

    ///
    VectorType cell_dof_values;
  };

  /**
   * @brief
   * Compute the volume enclosed by the 0 level of a level-set field
   * inside a cell. The inside volume is computed, defined by negative values of
   * the level-set field.
   *
   * @param[in] fe_point_evaluation FePointEvaluation
   *
   * @param[in] cell Cell for which the volume is computed
   *
   * @param[in] cell_dof_values cell DOFs value of the level set field
   *
   * @param[in] corr correction to apply to the DOF values (constant for all
   * DOFs)
   *
   * @param[in] n_quad_points number of quadrature points for the volume
   * integration faces
   *
   * @return cell-wise volume
   */
  template <int dim>
  double
  compute_cell_wise_volume(
    FEPointEvaluation<1, dim>                            &fe_point_evaluation,
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    Vector<double>     cell_dof_level_set_values,
    const double       corr,
    const unsigned int n_quad_points);

  /**
   * @brief
   * Compute the volume enclosed by the 0 level of a level set field
   * in the domain. The inside volume is computed, defined by negative values of
   * the level-set field.
   *
   * @param[in] mapping Mapping of the domain
   *
   * @param[in] dof_handler DofHandler associated to the triangulation on which
   * the volume is computed
   *
   * @param[in] level_set_vector Level-set vector
   *
   * @param[in] mpi_communicator MPI communicator
   *
   * @return Volume enclosed by the 0 level
   */
  template <int dim, typename VectorType>
  double
  compute_volume(const Mapping<dim>       &mapping,
                 const DoFHandler<dim>    &dof_handler,
                 const FiniteElement<dim> &fe,
                 const VectorType         &level_set_vector,
                 const MPI_Comm           &mpi_communicator);

  /**
   * @brief
   * Reconstruct the interface defined by the 0 level of a level set field
   * in the domain.
   *
   * @param[in] mapping Mapping of the domain
   *
   * @param[in] dof_handler DofHandler associated to the triangulation for which
   * the interface is reconstructed
   *
   * @param[in] fe Finite element
   *
   * @param[in] level_set_vector Level-set vector
   *
   * @param[in,out] interface_reconstruction_vertices Cell-wise map of the
   * reconstructed surface vertices. The map contains vectors storing the
   * vertices of the reconstructed surface for each intersected volume cell
   * (dim).
   *
   * @param[in,out] interface_reconstruction_cells Cell-wise map of the
   * reconstructed surface cells. The map contains vectors storing the cell
   * (dim-1) of the reconstructed surface for each intersected volume cell
   * (dim).
   *
   * @param[in,out] intersected_dofs Set of DOFs that belong to intersected
   * volume cell (dim).
   *
   */
  template <int dim, typename VectorType>
  void
  reconstruct_interface(
    const Mapping<dim>       &mapping,
    const DoFHandler<dim>    &dof_handler,
    const FiniteElement<dim> &fe,
    const VectorType         &level_set_vector,
    std::map<types::global_cell_index, std::vector<Point<dim>>>
      &interface_reconstruction_vertices,
    std::map<types::global_cell_index, std::vector<CellData<dim - 1>>>
                                      &interface_reconstruction_cells,
    std::set<types::global_dof_index> &intersected_dofs);

} // namespace InterfaceTools

#endif
