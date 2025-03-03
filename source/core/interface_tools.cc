// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/interface_tools.h>

template <int dim>
double
InterfaceTools::compute_cell_wise_volume(
  FEPointEvaluation<1, dim>                            &fe_point_evaluation,
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  Vector<double>     cell_dof_level_set_values,
  const double       corr,
  const unsigned int n_quad_points)
{
  const unsigned int n_dofs = cell_dof_level_set_values.size();

  // Initialize required variables to compute local volume
  const BoundingBox<dim> unit_box = create_unit_bounding_box<dim>();
  CellWiseFunction<dim>  signed_distance_function =
    CellWiseFunction<dim>(cell->get_fe().degree);

  hp::QCollection<1> q_collection;
  q_collection.push_back(QGauss<1>(n_quad_points));

  NonMatching::QuadratureGenerator<dim> quadrature_generator =
    NonMatching::QuadratureGenerator<dim>(q_collection);

  for (unsigned int j = 0; j < n_dofs; j++)
    {
      cell_dof_level_set_values[j] += corr;
    }

  signed_distance_function.set_active_cell(cell_dof_level_set_values);
  quadrature_generator.generate(signed_distance_function, unit_box);

  const Quadrature<dim> inside_quadrature =
    quadrature_generator.get_inside_quadrature();

  if (inside_quadrature.size() == 0)
    return 0.0;

  fe_point_evaluation.reinit(cell, inside_quadrature.get_points());

  double inside_cell_volume = 0.0;
  for (unsigned int q = 0; q < inside_quadrature.size(); q++)
    {
      /* Compute the volume int 1*JxW*dOmega. FEPointEvaluation.JxW() does not
      return the right thing.*/
      inside_cell_volume += fe_point_evaluation.jacobian(q).determinant() *
                            inside_quadrature.weight(q);
    }

  return inside_cell_volume;
}

template <int dim, typename VectorType>
double
InterfaceTools::compute_volume(const Mapping<dim>       &mapping,
                               const DoFHandler<dim>    &dof_handler,
                               const FiniteElement<dim> &fe,
                               const VectorType         &level_set_vector,
                               const double              iso_level,
                               const MPI_Comm           &mpi_communicator)
{
  FEPointEvaluation<1, dim> fe_point_evaluation(
    mapping, fe, update_jacobians | update_JxW_values);

  double volume = 0.0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          const unsigned int n_dofs_per_cell = cell->get_fe().n_dofs_per_cell();
          Vector<double>     cell_dof_level_set_values(n_dofs_per_cell);

          cell->get_dof_values(level_set_vector,
                               cell_dof_level_set_values.begin(),
                               cell_dof_level_set_values.end());

          const double level_set_correction = -iso_level;
          volume += compute_cell_wise_volume(fe_point_evaluation,
                                             cell,
                                             cell_dof_level_set_values,
                                             level_set_correction,
                                             cell->get_fe().degree + 1);
        }
    }
  volume = Utilities::MPI::sum(volume, mpi_communicator);

  return volume;
}

template double
InterfaceTools::compute_volume(const Mapping<2>       &mapping,
                               const DoFHandler<2>    &dof_handler,
                               const FiniteElement<2> &fe,
                               const Vector<double>   &level_set_vector,
                               const double            iso_level,
                               const MPI_Comm         &mpi_communicator);
template double
InterfaceTools::compute_volume(const Mapping<3>       &mapping,
                               const DoFHandler<3>    &dof_handler,
                               const FiniteElement<3> &fe,
                               const Vector<double>   &level_set_vector,
                               const double            iso_level,
                               const MPI_Comm         &mpi_communicator);

template <int dim, typename VectorType>
void
InterfaceTools::reconstruct_interface(
  const Mapping<dim>       &mapping,
  const DoFHandler<dim>    &dof_handler,
  const FiniteElement<dim> &fe,
  const VectorType         &level_set_vector,
  const double              iso_level,
  std::map<types::global_cell_index, std::vector<Point<dim>>>
    &interface_reconstruction_vertices,
  std::map<types::global_cell_index, std::vector<CellData<dim - 1>>>
                                    &interface_reconstruction_cells,
  std::set<types::global_dof_index> &intersected_dofs)
{
  GridTools::MarchingCubeAlgorithm<dim, VectorType> marching_cube(mapping,
                                                                  fe,
                                                                  1,
                                                                  1e-10);
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          std::vector<Point<dim>>        surface_vertices;
          std::vector<CellData<dim - 1>> surface_cells;

          marching_cube.process_cell(
            cell, level_set_vector, iso_level, surface_vertices, surface_cells);

          // If the cell is intersected, reconstruct the interface in it
          if (surface_vertices.size() != 0)
            {
              const unsigned int cell_index = cell->global_active_cell_index();

              // Store the interface reconstruction vertices and cells
              interface_reconstruction_vertices[cell_index] = surface_vertices;
              interface_reconstruction_cells[cell_index]    = surface_cells;

              // Store the dofs of the intersected volume cell
              std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
              cell->get_dof_indices(dof_indices);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  intersected_dofs.insert(dof_indices[i]);
                }
            }
        }
    }
}

template void
InterfaceTools::reconstruct_interface(
  const Mapping<2>       &mapping,
  const DoFHandler<2>    &dof_handler,
  const FiniteElement<2> &fe,
  const Vector<double>   &level_set_vector,
  const double            iso_level,
  std::map<types::global_cell_index, std::vector<Point<2>>>
    &interface_reconstruction_vertices,
  std::map<types::global_cell_index, std::vector<CellData<1>>>
                                    &interface_reconstruction_cells,
  std::set<types::global_dof_index> &intersected_dofs);
template void
InterfaceTools::reconstruct_interface(
  const Mapping<3>       &mapping,
  const DoFHandler<3>    &dof_handler,
  const FiniteElement<3> &fe,
  const Vector<double>   &level_set_vector,
  const double            iso_level,
  std::map<types::global_cell_index, std::vector<Point<3>>>
    &interface_reconstruction_vertices,
  std::map<types::global_cell_index, std::vector<CellData<2>>>
                                    &interface_reconstruction_cells,
  std::set<types::global_dof_index> &intersected_dofs);
