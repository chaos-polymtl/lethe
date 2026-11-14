// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/fine_scale_velocity_transfer.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/cell_data_transfer.templates.h>

#include <deal.II/lac/full_matrix.h>

#include <filesystem>
#include <fstream>
#include <utility>

template <int dim>
std::vector<std::vector<double>>
refine_fine_scale_velocity_cell_data(const FiniteElement<dim>  &fe,
                                     const std::vector<double> &parent_values)
{
  const unsigned int n_q_points = fe.n_dofs_per_cell();
  const unsigned int n_children = GeometryInfo<dim>::max_children_per_cell;
  AssertDimension(parent_values.size() % (n_q_points * dim), 0);
  const unsigned int n_history_entries =
    parent_values.size() / (n_q_points * dim);

  std::vector<std::vector<double>> children_values(
    n_children, std::vector<double>(parent_values.size(), 0.));

  for (unsigned int child = 0; child < n_children; ++child)
    {
      const FullMatrix<double> &prolongation =
        fe.get_prolongation_matrix(child);
      for (unsigned int history = 0; history < n_history_entries; ++history)
        for (unsigned int i = 0; i < n_q_points; ++i)
          for (unsigned int j = 0; j < n_q_points; ++j)
            for (unsigned int d = 0; d < dim; ++d)
              children_values[child][(history * n_q_points + i) * dim + d] +=
                prolongation(i, j) *
                parent_values[(history * n_q_points + j) * dim + d];
    }

  return children_values;
}

template <int dim>
std::vector<double>
coarsen_fine_scale_velocity_cell_data(
  const FiniteElement<dim>               &fe,
  const std::vector<std::vector<double>> &children_values)
{
  const unsigned int n_q_points = fe.n_dofs_per_cell();
  AssertDimension(children_values.size(),
                  GeometryInfo<dim>::max_children_per_cell);
  const unsigned int n_values = children_values[0].size();
  AssertDimension(n_values % (n_q_points * dim), 0);
  const unsigned int n_history_entries = n_values / (n_q_points * dim);

  std::vector<double> parent_values(n_values, 0.);

  for (unsigned int child = 0; child < children_values.size(); ++child)
    {
      AssertDimension(children_values[child].size(), n_values);
      const FullMatrix<double> &restriction = fe.get_restriction_matrix(child);
      for (unsigned int history = 0; history < n_history_entries; ++history)
        for (unsigned int i = 0; i < n_q_points; ++i)
          for (unsigned int j = 0; j < n_q_points; ++j)
            for (unsigned int d = 0; d < dim; ++d)
              parent_values[(history * n_q_points + i) * dim + d] +=
                restriction(i, j) *
                children_values[child][(history * n_q_points + j) * dim + d];
    }

  return parent_values;
}

template <int dim>
FineScaleVelocityTransfer<dim>::FineScaleVelocityTransfer(
  const Triangulation<dim> &triangulation,
  const MPI_Comm            mpi_communicator,
  const unsigned int        n_q_points_1d)
  : triangulation(
      dynamic_cast<const parallel::distributed::Triangulation<dim> *>(
        &triangulation))
  , mpi_communicator(mpi_communicator)
  , fe(QGauss<1>(n_q_points_1d))
{
  // The transfer across mesh adaptation relies on the support points of the
  // element being the quadrature points of a cell, in the same order.
  const QGauss<dim> quadrature(n_q_points_1d);
  const double      tolerance = 1e-12;
  bool points_match = fe.get_unit_support_points().size() == quadrature.size();
  for (unsigned int q = 0; points_match && q < quadrature.size(); ++q)
    points_match =
      fe.get_unit_support_points()[q].distance(quadrature.point(q)) < tolerance;
  AssertThrow(points_match,
              ExcMessage("The support points of the element used to transfer "
                         "the history of the dynamic RBVMS fine-scale "
                         "velocity do not match the quadrature points of the "
                         "cells."));
}

template <int dim>
std::unique_ptr<parallel::distributed::CellDataTransfer<
  dim,
  dim,
  typename FineScaleVelocityTransfer<dim>::CellData>>
FineScaleVelocityTransfer<dim>::create_cell_data_transfer() const
{
  AssertThrow(this->triangulation != nullptr,
              ExcMessage("The history of the dynamic RBVMS fine-scale "
                         "velocity can only be transferred with a "
                         "parallel::distributed::Triangulation."));

  // Every locally owned cell holds the same number of values, so the data
  // has a fixed size.
  return std::make_unique<
    parallel::distributed::CellDataTransfer<dim, dim, CellData>>(
    *this->triangulation,
    false,
    [this](const typename Triangulation<dim>::cell_iterator &,
           const std::vector<double> parent_values) {
      return refine_fine_scale_velocity_cell_data(this->fe, parent_values);
    },
    [this](const typename Triangulation<dim>::cell_iterator &,
           const std::vector<std::vector<double>> &children_values) {
      return coarsen_fine_scale_velocity_cell_data(this->fe, children_values);
    });
}

template <int dim>
void
FineScaleVelocityTransfer<dim>::attach_to_checkpoint(
  const std::string &checkpoint_file_prefix,
  const unsigned int n_values_per_cell,
  CellData           cell_data)
{
  // The number of values per cell is always written, also when it is zero, so
  // that the file of a previous checkpoint with the same prefix is overwritten
  // and the restart never looks for a history that is not attached.
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      std::ofstream output(checkpoint_file_prefix + ".fine_scale_velocity");
      output << n_values_per_cell << std::endl;
    }

  if (n_values_per_cell == 0)
    return;

  // The cell-wise data and the object attaching it must stay alive until the
  // triangulation is saved, so they are members released by clear().
  this->attached_cell_data = std::move(cell_data);
  this->cell_data_transfer = this->create_cell_data_transfer();
  this->cell_data_transfer->prepare_for_serialization(this->attached_cell_data);
}

template <int dim>
typename FineScaleVelocityTransfer<dim>::CellData
FineScaleVelocityTransfer<dim>::read_from_checkpoint(
  const std::string        &checkpoint_file_prefix,
  const unsigned int        n_values_per_cell,
  const ConditionalOStream &pcout)
{
  // Number of values per cell of the history attached to the checkpoint. It is
  // read by the first process and broadcast, so that every process takes the
  // same decisions below. A checkpoint written without this file holds no
  // history.
  const std::string filename = checkpoint_file_prefix + ".fine_scale_velocity";
  unsigned int      n_values_per_cell_in_checkpoint = 0;
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0 &&
      std::filesystem::exists(filename))
    {
      std::ifstream input(filename);
      if (!(input >> n_values_per_cell_in_checkpoint))
        n_values_per_cell_in_checkpoint = numbers::invalid_unsigned_int;
    }
  n_values_per_cell_in_checkpoint =
    Utilities::MPI::broadcast(this->mpi_communicator,
                              n_values_per_cell_in_checkpoint);
  AssertThrow(n_values_per_cell_in_checkpoint != numbers::invalid_unsigned_int,
              ExcMessage("The number of values per cell of the history of the "
                         "dynamic RBVMS fine-scale velocity could not be read "
                         "from the checkpoint file " +
                         filename + "."));

  if (n_values_per_cell_in_checkpoint == 0)
    {
      if (n_values_per_cell > 0)
        pcout << "The checkpoint holds no history of the dynamic RBVMS "
                 "fine-scale velocity, which restarts from zero."
              << std::endl;
      return {};
    }

  // The history attached to the triangulation is read even if the simulation
  // does not use it, since all the data attached to a triangulation must be
  // read before it is saved again.
  auto     transfer = this->create_cell_data_transfer();
  CellData cell_data(this->triangulation->n_active_cells());
  transfer->deserialize(cell_data);

  if (n_values_per_cell == 0)
    {
      pcout << "The history of the dynamic RBVMS fine-scale velocity "
               "held by the checkpoint is not used, since the "
               "simulation does not use the dynamic velocity subscale "
               "in a transient simulation."
            << std::endl;
      return {};
    }

  AssertThrow(
    n_values_per_cell_in_checkpoint == n_values_per_cell,
    ExcMessage(
      "The checkpoint holds " +
      std::to_string(n_values_per_cell_in_checkpoint) +
      " values per cell of the history of the dynamic RBVMS fine-scale "
      "velocity, whereas the simulation expects " +
      std::to_string(n_values_per_cell) +
      " (number of previous solutions of the time-stepping method times "
      "number of quadrature points times dimension). Restart with the "
      "time-stepping method and the velocity degree of the simulation that "
      "wrote the checkpoint."));

  return cell_data;
}

template <int dim>
void
FineScaleVelocityTransfer<dim>::prepare_for_coarsening_and_refinement(
  CellData cell_data)
{
  this->attached_cell_data = std::move(cell_data);
  this->cell_data_transfer = this->create_cell_data_transfer();
  this->cell_data_transfer->prepare_for_coarsening_and_refinement(
    this->attached_cell_data);
}

template <int dim>
typename FineScaleVelocityTransfer<dim>::CellData
FineScaleVelocityTransfer<dim>::unpack_after_coarsening_and_refinement()
{
  AssertThrow(this->cell_data_transfer != nullptr,
              ExcMessage("The history of the dynamic RBVMS fine-scale "
                         "velocity was not prepared for the mesh adaptation. "
                         "Call prepare_for_coarsening_and_refinement() before "
                         "the mesh is adapted."));

  CellData cell_data(this->triangulation->n_active_cells());
  this->cell_data_transfer->unpack(cell_data);
  this->clear();

  return cell_data;
}

template <int dim>
void
FineScaleVelocityTransfer<dim>::clear()
{
  this->cell_data_transfer.reset();
  this->attached_cell_data.clear();
  this->attached_cell_data.shrink_to_fit();
}

template std::vector<std::vector<double>>
refine_fine_scale_velocity_cell_data<2>(const FiniteElement<2> &,
                                        const std::vector<double> &);
template std::vector<std::vector<double>>
refine_fine_scale_velocity_cell_data<3>(const FiniteElement<3> &,
                                        const std::vector<double> &);
template std::vector<double>
coarsen_fine_scale_velocity_cell_data<2>(
  const FiniteElement<2> &,
  const std::vector<std::vector<double>> &);
template std::vector<double>
coarsen_fine_scale_velocity_cell_data<3>(
  const FiniteElement<3> &,
  const std::vector<std::vector<double>> &);

template class FineScaleVelocityTransfer<2>;
template class FineScaleVelocityTransfer<3>;
