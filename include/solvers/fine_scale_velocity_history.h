// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fine_scale_velocity_history_h
#define lethe_fine_scale_velocity_history_h

#include <deal.II/base/tensor.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <vector>

using namespace dealii;

/**
 * @brief History of the fine-scale velocity of the dynamic RBVMS velocity
 * subscale at the quadrature points of the cells, for the matrix-based
 * solvers.
 *
 * The fine-scale velocity u' is advanced in time with the BDF scheme of the
 * velocity, u'_{n+1} = -tau_t (r_M + m'), where the memory
 * m' = sum_{j>=1} beta_j u'_{n+1-j} gathers the fine-scale velocities of the
 * previous time steps. This class stores these previous fine-scale velocities
 * cell by cell. The values of a cell are ordered by history entry (entry 0 is
 * the most recent time step), quadrature point and velocity component, which
 * is the layout used by FineScaleVelocityTransfer to write the history to a
 * checkpoint and to transfer it across mesh adaptation.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class FineScaleVelocityHistory
{
public:
  /**
   * @brief Cell-wise data. Entry cell->active_cell_index() holds the values of
   * a locally owned cell, and the entries of the other cells are empty.
   */
  using CellData = std::vector<std::vector<double>>;

  /**
   * @brief Allocate the history, with a fine-scale velocity equal to zero.
   *
   * @param[in] triangulation Triangulation holding the cells.
   *
   * @param[in] n_history_entries Number of previous fine-scale velocities
   * stored, that is the number of previous solutions of the time-stepping
   * method.
   *
   * @param[in] n_q_points Number of quadrature points per cell.
   */
  void
  reinit(const Triangulation<dim> &triangulation,
         const unsigned int        n_history_entries,
         const unsigned int        n_q_points);

  /**
   * @brief Return the number of values stored per cell, that is the number of
   * history entries times the number of quadrature points times dim.
   *
   * @return Number of values of the history per cell.
   */
  unsigned int
  n_values_per_cell() const
  {
    return n_history_entries * n_q_points * dim;
  }

  /**
   * @brief Calculate the memory m' = sum_{j>=1} beta_j u'_{n+1-j} of the
   * fine-scale velocity at a quadrature point.
   *
   * @param[in] active_cell_index Active cell index of a locally owned cell.
   *
   * @param[in] q Index of the quadrature point.
   *
   * @param[in] bdf_coefficients Coefficients beta_j of the BDF scheme.
   *
   * @param[in] n_previous_solutions Number of previous solutions used by the
   * BDF scheme of the current time step.
   *
   * @return Memory of the fine-scale velocity.
   */
  Tensor<1, dim>
  calculate_memory(const unsigned int    active_cell_index,
                   const unsigned int    q,
                   const Vector<double> &bdf_coefficients,
                   const unsigned int    n_previous_solutions) const
  {
    AssertIndexRange(active_cell_index, cell_data.size());
    AssertIndexRange(n_previous_solutions, n_history_entries + 1);
    const std::vector<double> &values = cell_data[active_cell_index];
    Tensor<1, dim>             memory;
    for (unsigned int p = 0; p < n_previous_solutions; ++p)
      for (int d = 0; d < dim; ++d)
        memory[d] +=
          bdf_coefficients[p + 1] * values[(p * n_q_points + q) * dim + d];
    return memory;
  }

  /**
   * @brief Store the fine-scale velocity of the time step that just converged
   * at a quadrature point. It only enters the history when advance() is
   * called. Different cells can be written concurrently.
   *
   * @param[in] active_cell_index Active cell index of a locally owned cell.
   *
   * @param[in] q Index of the quadrature point.
   *
   * @param[in] fine_scale_velocity Fine-scale velocity u'_{n+1}.
   */
  void
  store_new_value(const unsigned int    active_cell_index,
                  const unsigned int    q,
                  const Tensor<1, dim> &fine_scale_velocity)
  {
    AssertIndexRange(active_cell_index, new_values.size());
    for (int d = 0; d < dim; ++d)
      new_values[active_cell_index][q * dim + d] = fine_scale_velocity[d];
  }

  /**
   * @brief Shift the history by one time step: the values given to
   * store_new_value() become entry 0 and the oldest entry is discarded.
   */
  void
  advance();

  /**
   * @brief Tell whether the fine-scale velocity computed by the assembler
   * must be stored with store_new_value().
   *
   * @return True if the fine-scale velocity must be stored.
   */
  bool
  is_recording() const
  {
    return recording;
  }

  /**
   * @brief Set whether the fine-scale velocity computed by the assembler must
   * be stored with store_new_value().
   *
   * @param[in] is_recording True if the fine-scale velocity must be stored.
   */
  void
  set_recording(const bool is_recording)
  {
    recording = is_recording;
  }

  /**
   * @brief Return the cell-wise data of the history.
   *
   * @return Cell-wise data of the history.
   */
  const CellData &
  get_cell_data() const
  {
    return cell_data;
  }

  /**
   * @brief Set the history from cell-wise data, for instance when a simulation
   * is restarted or after the mesh has been adapted. The history must have
   * been allocated by reinit() for the same triangulation and the same number
   * of values per cell.
   *
   * @param[in] new_cell_data Cell-wise data of the history.
   */
  void
  set_cell_data(const CellData &new_cell_data);

private:
  /// Number of previous fine-scale velocities stored.
  unsigned int n_history_entries = 0;

  /// Number of quadrature points per cell.
  unsigned int n_q_points = 0;

  /// History of the fine-scale velocity.
  CellData cell_data;

  /// Fine-scale velocity of the time step that just converged.
  CellData new_values;

  /// Whether the assembler must store the fine-scale velocity.
  bool recording = false;
};

#endif
