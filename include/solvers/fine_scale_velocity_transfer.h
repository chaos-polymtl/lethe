// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fine_scale_velocity_transfer_h
#define lethe_fine_scale_velocity_transfer_h

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>

#include <deal.II/distributed/cell_data_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/tria.h>

#include <memory>
#include <string>
#include <vector>

using namespace dealii;

/**
 * @brief Interpolate the history of the dynamic RBVMS fine-scale velocity of a
 * cell to its children. The values of a cell are ordered by history entry,
 * quadrature point and velocity component. Since the quadrature points of a
 * cell are the support points of @p fe, the values of every history entry and
 * component are the nodal values of a polynomial, which is evaluated at the
 * quadrature points of the children with the prolongation matrices of @p fe.
 *
 * @tparam dim Number of spatial dimensions.
 *
 * @param[in] fe Scalar discontinuous finite element whose support points are
 * the quadrature points of a cell.
 *
 * @param[in] parent_values Values of the parent cell.
 *
 * @return Values of each child of an isotropically refined cell.
 */
template <int dim>
std::vector<std::vector<double>>
refine_fine_scale_velocity_cell_data(const FiniteElement<dim>  &fe,
                                     const std::vector<double> &parent_values);

/**
 * @brief Project the history of the dynamic RBVMS fine-scale velocity of the
 * children of a cell to this cell. This is the counterpart of
 * refine_fine_scale_velocity_cell_data(): the values of every history entry
 * and component are projected (L2 projection over the parent cell) with the
 * restriction matrices of @p fe.
 *
 * @tparam dim Number of spatial dimensions.
 *
 * @param[in] fe Scalar discontinuous finite element whose support points are
 * the quadrature points of a cell.
 *
 * @param[in] children_values Values of each child of an isotropically refined
 * cell.
 *
 * @return Values of the parent cell.
 */
template <int dim>
std::vector<double>
coarsen_fine_scale_velocity_cell_data(
  const FiniteElement<dim>               &fe,
  const std::vector<std::vector<double>> &children_values);

/**
 * @brief Transfer of the history of the dynamic RBVMS fine-scale velocity,
 * which is stored at the quadrature points of the cells, to a checkpoint and
 * across mesh adaptation.
 *
 * The history is exchanged as cell-wise data: one vector of values per active
 * cell, ordered by history entry, quadrature point and velocity component. The
 * data is attached to the cells of a parallel::distributed::Triangulation, so
 * it follows them when the mesh is repartitioned, refined or coarsened.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class FineScaleVelocityTransfer
{
public:
  /**
   * @brief Cell-wise data. Entry cell->active_cell_index() holds the values of
   * a locally owned cell, ordered by history entry, quadrature point and
   * velocity component. The entries of the other cells are empty.
   */
  using CellData = std::vector<std::vector<double>>;

  /**
   * @brief Constructor.
   *
   * @param[in] triangulation Triangulation holding the cells. It must be a
   * parallel::distributed::Triangulation for the data to be transferred.
   *
   * @param[in] mpi_communicator MPI communicator of the triangulation.
   *
   * @param[in] n_q_points_1d Number of Gauss quadrature points per direction
   * at which the fine-scale velocity is stored.
   */
  FineScaleVelocityTransfer(const Triangulation<dim> &triangulation,
                            const MPI_Comm            mpi_communicator,
                            const unsigned int        n_q_points_1d);

  /**
   * @brief Attach the cell-wise data to the triangulation, so that it is saved
   * with it in a checkpoint. The number of values per cell is always written
   * to the file <prefix>.fine_scale_velocity, also when it is zero, in which
   * case nothing is attached. The data is kept by this object, and clear()
   * must be called once the triangulation has been saved.
   *
   * @param[in] checkpoint_file_prefix Prefix of the checkpoint files.
   *
   * @param[in] n_values_per_cell Number of values per cell, which must be the
   * same on every process.
   *
   * @param[in] cell_data Cell-wise data to attach.
   */
  void
  attach_to_checkpoint(const std::string &checkpoint_file_prefix,
                       const unsigned int n_values_per_cell,
                       CellData           cell_data);

  /**
   * @brief Read the cell-wise data attached to the triangulation loaded from a
   * checkpoint. If the checkpoint holds no data (zero values per cell, or
   * checkpoint written without the file <prefix>.fine_scale_velocity), or if
   * the simulation does not expect any, the returned data is empty. Data that
   * is not expected is still read and discarded, since all the data attached
   * to a triangulation must be read.
   *
   * @param[in] checkpoint_file_prefix Prefix of the checkpoint files.
   *
   * @param[in] n_values_per_cell Number of values per cell expected by the
   * simulation, which must be the same on every process.
   *
   * @param[in] pcout Stream used to report how the restart is carried out.
   *
   * @return Cell-wise data of the checkpoint, or an empty container if there
   * is nothing to restart from.
   */
  CellData
  read_from_checkpoint(const std::string        &checkpoint_file_prefix,
                       const unsigned int        n_values_per_cell,
                       const ConditionalOStream &pcout);

  /**
   * @brief Register the cell-wise data for the transfer across the mesh
   * adaptation that is about to be executed. The data is kept by this object
   * until unpack_after_coarsening_and_refinement() is called.
   *
   * @param[in] cell_data Cell-wise data on the mesh before adaptation.
   */
  void
  prepare_for_coarsening_and_refinement(CellData cell_data);

  /**
   * @brief Return the cell-wise data on the adapted mesh. The data of a
   * refined cell is interpolated to its children, and the data of coarsened
   * cells is projected to their parent.
   *
   * @return Cell-wise data on the mesh after adaptation.
   */
  CellData
  unpack_after_coarsening_and_refinement();

  /**
   * @brief Release the cell-wise data kept by this object and the object that
   * attached it to the triangulation.
   */
  void
  clear();

private:
  /**
   * @brief Create the object that attaches the cell-wise data to the cells of
   * the triangulation.
   *
   * @return Object transferring the cell-wise data.
   */
  std::unique_ptr<parallel::distributed::CellDataTransfer<dim, dim, CellData>>
  create_cell_data_transfer() const;

  /// Triangulation holding the cells, or nullptr if it is not a
  /// parallel::distributed::Triangulation.
  const parallel::distributed::Triangulation<dim> *triangulation;

  /// MPI communicator of the triangulation.
  const MPI_Comm mpi_communicator;

  /// Scalar discontinuous element whose support points are the quadrature
  /// points of a cell.
  const FE_DGQArbitraryNodes<dim> fe;

  /// Cell-wise data attached to the triangulation.
  CellData attached_cell_data;

  /// Object attaching attached_cell_data to the triangulation.
  std::unique_ptr<parallel::distributed::CellDataTransfer<dim, dim, CellData>>
    cell_data_transfer;
};

#endif
