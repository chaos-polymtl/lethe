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


template <int dim, typename VectorType>
void
InterfaceTools::SignedDistanceSolver<dim, VectorType>::setup_dofs()
{
  const MPI_Comm mpi_communicator = dof_handler.get_communicator();

  dof_handler.distribute_dofs(*this->fe);

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
  locally_active_dofs   = DoFTools::extract_locally_active_dofs(dof_handler);

  level_set.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  signed_distance.reinit(locally_owned_dofs,
                         locally_active_dofs,
                         mpi_communicator);
  signed_distance_with_ghost.reinit(locally_owned_dofs,
                                    locally_active_dofs,
                                    mpi_communicator);

  distance.reinit(locally_owned_dofs, locally_active_dofs, mpi_communicator);
  distance_with_ghost.reinit(locally_owned_dofs,
                             locally_active_dofs,
                             mpi_communicator);

  volume_correction.reinit(locally_owned_dofs,
                           locally_active_dofs,
                           mpi_communicator);

  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();
}

template <int dim, typename VectorType>
void
InterfaceTools::SignedDistanceSolver<dim, VectorType>::
  set_level_set_from_background_mesh(
    const DoFHandler<dim> &background_dof_handler,
    const VectorType      &background_level_set_vector)
{
  const MPI_Comm mpi_communicator = dof_handler.get_communicator();

  VectorType tmp_local_level_set(this->locally_owned_dofs, mpi_communicator);

  FETools::interpolate(background_dof_handler,
                       background_level_set_vector,
                       dof_handler,
                       tmp_local_level_set);

  level_set = tmp_local_level_set;
}

template <int dim, typename VectorType>
void
InterfaceTools::SignedDistanceSolver<dim, VectorType>::solve()
{
  // Gain the writing right.
  zero_out_ghost_values();

  // Clear maps and sets
  interface_reconstruction_vertices.clear();
  interface_reconstruction_cells.clear();
  intersected_dofs.clear();

  // Initialize local distance vetors.
  initialize_distance();

  // Identify intersected cells and compute the interface reconstruction.
  InterfaceTools::reconstruct_interface(*mapping,
                                        dof_handler,
                                        *fe,
                                        level_set,
                                        iso_level,
                                        interface_reconstruction_vertices,
                                        interface_reconstruction_cells,
                                        intersected_dofs);

  /* Compute the distance for the dofs of the intersected cells (the ones in
  the intersected_dofs set). They correspond to the first neighbor dofs.*/
  compute_first_neighbors_distance();

  /* Compute signed distance from distance (only first neighbors have an
  updated value)*/
  compute_signed_distance_from_distance();

  // Conserve local and global volume
  // compute_cell_wise_volume_correction();
  // conserve_global_volume();

  /* Compute the distance for the dofs of the rest of the mesh. They
  correspond to the second neighbors dofs. */
  compute_second_neighbors_distance();

  // Compute signed distance from distance (all DOFs have updated value)
  compute_signed_distance_from_distance();

  // Update ghost values to regain reading ability.
  update_ghost_values();
}

template <int dim, typename VectorType>
void
InterfaceTools::SignedDistanceSolver<dim, VectorType>::initialize_distance()
{
  /* Initialization of the active dofs to the max distance we want to
   redistanciate. It requires a loop on the active dofs to initialize also
   the distance value of the ghost dofs. Otherwise, they are set to 0.0 by
   default. */
  for (auto p : this->locally_active_dofs)
    {
      distance(p)            = max_distance;
      distance_with_ghost(p) = max_distance;

      const double level_set_value  = level_set(p);
      signed_distance(p)            = max_distance * sgn(level_set_value);
      signed_distance_with_ghost(p) = max_distance * sgn(level_set_value);
    }
}

template <int dim, typename VectorType>
VectorType &
InterfaceTools::SignedDistanceSolver<dim, VectorType>::get_signed_distance()
{
  const MPI_Comm mpi_communicator = dof_handler.get_communicator();

  VectorType tmp_local_level_set(this->locally_owned_dofs, mpi_communicator);

  // Loop on the DoFs to be compatible with the difference in vector type
  // between level_set and signed_distance
  for (auto p : this->locally_owned_dofs)
    {
      tmp_local_level_set(p) = signed_distance(p);
    }

  level_set = tmp_local_level_set;

  return level_set;
}

template <int dim, typename VectorType>
void
InterfaceTools::SignedDistanceSolver<dim, VectorType>::zero_out_ghost_values()
{
  /* To have the right to write in a LinearAlgebra::distributed::Vector, we
  have to zero out the ghost values.*/
  signed_distance.zero_out_ghost_values();
  signed_distance_with_ghost.zero_out_ghost_values();
  distance.zero_out_ghost_values();
  distance_with_ghost.zero_out_ghost_values();
  volume_correction.zero_out_ghost_values();
}

template <int dim, typename VectorType>
void
InterfaceTools::SignedDistanceSolver<dim, VectorType>::update_ghost_values()
{
  /* To have the right to read a LinearAlgebra::distributed::Vector, we have
  to update the ghost values.*/
  signed_distance.update_ghost_values();
  signed_distance_with_ghost.update_ghost_values();
  distance.update_ghost_values();
  distance_with_ghost.update_ghost_values();
  volume_correction.update_ghost_values();
}

template <int dim, typename VectorType>
void
InterfaceTools::SignedDistanceSolver<dim, VectorType>::exchange_distance()
{
  // Exchange and "select" min value between the processes
  distance.compress(VectorOperation::min);

  // Update local ghost (distance becomes read only)
  distance.update_ghost_values();

  /* Copy distance to distance_with_ghost to keep the knowledge of local ghost
  values and to have a read only version of the vector*/
  distance_with_ghost = distance;

  /* Zero out ghost DOFs to regain write functionalities in distance (it
  becomes write only, that is why we need distance_with_ghost - to read the
  ghost values in it).*/
  distance.zero_out_ghost_values();

  /* Copy the ghost values back in distance (zero_out_ghost_values() puts
  zeros in ghost DOFs)*/
  for (auto p : this->locally_active_dofs)
    {
      /* We need to have the ghost values in distance for future
      compress(VectorOperation::min) operation*/
      distance(p) = distance_with_ghost(p);
    }
}

template <int dim, typename VectorType>
void
InterfaceTools::SignedDistanceSolver<dim, VectorType>::
  compute_first_neighbors_distance()
{
  /* The signed distance for the first neighbors (the Dofs belonging to the
   * cells intersected by the reconstructed interface. This is a brute force
   * distance computation, meaning the distance is computed geometrically, as
   * presented by Ausas (2010). */

  // DoF coordinates
  std::map<types::global_dof_index, Point<dim>> dof_support_points =
    DoFTools::map_dofs_to_support_points(*mapping, dof_handler);

  // Loop over the intersected cells (volume cells)
  for (auto &intersected_cell : interface_reconstruction_cells)
    {
      const unsigned int cell_index = intersected_cell.first;

      // Create interface recontruction triangulation (surface triangulation)
      // in the intersected volume cell
      std::vector<Point<dim>> surface_vertices =
        interface_reconstruction_vertices.at(cell_index);
      std::vector<CellData<dim - 1>> surface_cells = intersected_cell.second;

      Triangulation<dim - 1, dim> surface_triangulation;
      surface_triangulation.create_triangulation(surface_vertices,
                                                 surface_cells,
                                                 {});

      /* Loop over all DoFs of the volume mesh belonging to a intersected
       * volume cell. This is more expensive, but it is required to have the
       * the right signed distance approximation for the first neighbors.*/
      for (const unsigned int &intersected_dof : intersected_dofs)
        {
          const Point<dim> y = dof_support_points.at(intersected_dof);

          /* Loop over the surface cells of the interface reconstruction in
           * the volume cell. In 2D, there is only 1 surface cell (line),
           * while in 3D, it can vary from 1 to 4 or 5 (triangles), depending
           * on the marching cube algorithm.*/
          for (const auto &surface_cell :
               surface_triangulation.active_cell_iterators())
            {
              // Store the current surface cell vertex coordinates
              unsigned int surface_cell_n_vertices = surface_cell->n_vertices();
              std::vector<Point<dim>> surface_cell_vertices(
                surface_cell_n_vertices);
              for (unsigned int p = 0; p < surface_cell_n_vertices; p++)
                {
                  surface_cell_vertices[p] = surface_cell->vertex(p);
                }

              // Compute the geometrical distance between the surface cell
              // (line in 2D, triangle in 3D) and the DoF
              double D = LetheGridTools::find_point_triangle_distance(
                surface_cell_vertices, y);

              // Select the minimum distance
              distance(intersected_dof) =
                std::min(std::abs(distance(intersected_dof)), std::abs(D));
            }
        }
    }
  exchange_distance();
}

template <int dim, typename VectorType>
void
InterfaceTools::SignedDistanceSolver<dim, VectorType>::
  compute_second_neighbors_distance()
{
  /* The signed distance for the second neighbors (the cells not intersected
   * by the interface is resolved according the minimization problem
   * presented by Ausas (2010). The method looks for the point in the opposite
   * faces of each second neighbor DoFs that minimizes the distance to the
   * interface. It works in a similar manner as a marching algorithm from the
   * knowledge of the signed distance for the interface first neighbors. */

  const MPI_Comm mpi_communicator = dof_handler.get_communicator();

  const unsigned int n_opposite_faces_per_dofs = dim;
  const unsigned int dofs_per_cell             = fe->n_dofs_per_cell();

  std::map<types::global_dof_index, Point<dim>> dof_support_points =
    DoFTools::map_dofs_to_support_points(*mapping, dof_handler);

  Point<dim - 1> ref_face_center_point = Point<dim - 1>();
  ref_face_center_point(0)             = 0.5;
  if constexpr (dim == 3)
    ref_face_center_point(1) = 0.5;

  FEPointEvaluation<1, dim> fe_point_evaluation(
    *mapping, *fe, update_values | update_gradients | update_jacobians);

  /* The method is iterative, hence, we solve as long as the distance
  approximation changes for at least one dof. We use the flag change to
  track this change. */
  bool change = true;

  /* The count corresponds to how many times we iterate. In fact, it
  correspond to the number of cell layers (starting from the interface)
  that the approximation of the distance is known. */
  int count = 0;
  while (change)
    {
      pcout << "Redistanciating layer " << count << std::endl;
      change = false;

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              const unsigned int cell_index = cell->global_active_cell_index();

              // If the cell is intersected, the distance is already computed.
              if (interface_reconstruction_vertices.find(cell_index) !=
                  interface_reconstruction_vertices.end())
                {
                  continue;
                }

              std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

              cell->get_dof_indices(dof_indices);
              std::vector<double> cell_dof_values(dofs_per_cell);
              cell->get_dof_values(distance_with_ghost,
                                   cell_dof_values.begin(),
                                   cell_dof_values.end());

              // Loop over the cell's Dofs
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  /* If the dof belongs to an intersected cell, the distance
                  is already computed */
                  if (intersected_dofs.find(dof_indices[i]) !=
                      intersected_dofs.end())
                    {
                      continue;
                    }

                  // Get opposite faces
                  std::vector<unsigned int> dof_opposite_faces(
                    n_opposite_faces_per_dofs);
                  get_dof_opposite_faces(i, dof_opposite_faces);

                  // Get the real coordinates of the current DoF I
                  const Point<dim> x_I_real =
                    dof_support_points.at(dof_indices[i]);

                  // Loop on opposite faces F_J
                  for (unsigned int j = 0; j < n_opposite_faces_per_dofs; ++j)
                    {
                      /* The minimization problem is: Find x in the face F_J
                      (opposite to the DoF of interest I) such that:

                        |d|_{x_I} = min(phi(x) +|x_I - x|)

                      where x_I is the coord of the DoF I, phi(x) is the
                      distance (not signed) at the point x, belonging to the
                      face F_J. Here, we solve the problem in the reference
                      space (dim - 1).
                      */

                      // Initialize required variables
                      Point<dim> x_n_ref = transform_ref_face_point_to_ref_cell(
                        ref_face_center_point, dof_opposite_faces[j]);
                      Point<dim> x_n_real;

                      double correction_norm = 1.0;
                      int    newton_it       = 0;

                      // Check to constrain the solution in the face F_J
                      int outside_check = 0;

                      // Solve the minimization problem with Newton method
                      // using a numerical jacobian
                      while (correction_norm > 1e-10 && outside_check < 3 &&
                             newton_it < 100)
                        {
                          /* Set stencil for numerical jacobian computation.
                           The entries of the vector are the following:
                                    4

                               1    0    2

                                    3
                          The entry 0 is the current evaluation point. */

                          const double            perturbation = 0.01;
                          std::vector<Point<dim>> stencil_ref =
                            compute_numerical_jacobian_stencil(
                              x_n_ref, dof_opposite_faces[j], perturbation);

                          std::vector<Point<dim>>     stencil_real(2 * dim - 1);
                          std::vector<Tensor<1, dim>> distance_gradients(
                            2 * dim - 1);
                          std::vector<DerivativeForm<1, dim, dim>>
                            cell_transformation_jacobians(2 * dim - 1);
                          std::vector<DerivativeForm<1, dim - 1, dim>>
                            face_transformation_jacobians(2 * dim - 1);


                          /* Prepare FEPointEvaluation to compute value and
                          gradient at the stencil points*/
                          fe_point_evaluation.reinit(cell, stencil_ref);
                          fe_point_evaluation.evaluate(
                            cell_dof_values, EvaluationFlags::gradients);

                          // Get the required values at each stencil point
                          for (unsigned int k = 0; k < 2 * dim - 1; k++)
                            {
                              stencil_real[k] =
                                fe_point_evaluation.quadrature_point(k);
                              distance_gradients[k] =
                                fe_point_evaluation.get_gradient(k);
                              cell_transformation_jacobians[k] =
                                fe_point_evaluation.jacobian(k);
                              get_face_transformation_jacobian(
                                cell_transformation_jacobians[k],
                                dof_opposite_faces[j],
                                face_transformation_jacobians[k]);
                            }

                          /* Compute the jacobian matrix. The Ax=b system is
                          formulated as the dim-1 system. We solve for the
                          correction in the reference face. */
                          LAPACKFullMatrix<double> jacobian_matrix(dim - 1,
                                                                   dim - 1);
                          compute_numerical_jacobian(
                            stencil_real,
                            x_I_real,
                            distance_gradients,
                            face_transformation_jacobians,
                            perturbation,
                            jacobian_matrix);

                          const Tensor<1, dim> x_n_to_x_I_real =
                            x_I_real - stencil_real[0];

                          // Compute the right hand side.
                          Tensor<1, dim - 1> residual_n;
                          compute_residual(x_n_to_x_I_real,
                                           distance_gradients[0],
                                           face_transformation_jacobians[0],
                                           residual_n);

                          // Convert the right hand side to the right format
                          // for the linear solver
                          Vector<double> residual_n_vec(dim - 1);
                          residual_n.unroll(residual_n_vec.begin(),
                                            residual_n_vec.end());
                          residual_n_vec *= -1.0;

                          jacobian_matrix.set_property(LAPACKSupport::general);

                          /* Factorize and solve the matrix. The correction is
                          put back in residual_n_vec. */
                          jacobian_matrix.compute_lu_factorization();
                          jacobian_matrix.solve(residual_n_vec);

                          // Compute the norm of the correction
                          correction_norm = residual_n_vec.l2_norm();

                          /* Transform the dim-1 correction (in the reference
                          face) to dim (in the reference cell) */
                          Tensor<1, dim> correction =
                            transform_ref_face_correction_to_ref_cell(
                              residual_n_vec, dof_opposite_faces[j]);

                          /* Compute the solution (the point x_n_ref on the
                          face minimizing the distance)*/
                          Point<dim> x_n_p1_ref = stencil_ref[0] + correction;

                          /* Relax the correction if it brings us outside of
                          the cell */
                          double relaxation = 1.0;

                          /* Check if the Newton method results in a solution
                          outside the face. For example in 3D, we could have:
                               _____________
                              |             |     solution
                              |             |    *
                              |             |
                              |             |
                              |             |
                              |_____________|

                          Each time it does, we relax the scheme to bring
                          back the estimation of the solution in the face:
                               _____________
                              |             | relaxed solution
                              |           * |
                              |             |
                              |             |
                              |             |
                              |_____________|

                          If the solution is outside the face more than three
                          times, we constraint the solution on the right
                          boundary of the face:
                               _____________
                              |             |         real
                              |  constraint *     * solution
                              |   solution  |
                              |             |
                              |             |
                              |_____________|

                          */

                          /* Flag indicating if the correction brings us
                          outside of the cell.*/
                          bool check = false;
                          for (unsigned int k = 0; k < dim; ++k)
                            {
                              if (x_n_p1_ref[k] > 1.0 + 1e-12 ||
                                  x_n_p1_ref[k] < 0.0 - 1e-12)
                                {
                                  check = true;

                                  /* Set the correction to put the solution on
                                  the face boundary. Select the minimum
                                  relaxation of the all directions to ensure
                                  the solution stays inside the face.*/
                                  if (correction[k] > 1e-12)
                                    {
                                      relaxation =
                                        std::min((1.0 - x_n_ref[k]) /
                                                   (correction[k] + 1e-12),
                                                 relaxation);
                                    }
                                  else if (correction[k] < -1e-12)
                                    {
                                      relaxation =
                                        std::min((0.0 - x_n_ref[k]) /
                                                   (correction[k] + 1e-12),
                                                 relaxation);
                                    }
                                }
                            }

                          // Increment the outside_check if the correction
                          // brought us outside the face
                          if (check)
                            outside_check += 1;

                          // Re-compute the solution with the relaxation
                          x_n_p1_ref = stencil_ref[0] + relaxation * correction;

                          // Transform the solution from reference to the real
                          // cell. This could be improved to not call
                          // fe_point_evaluation.
                          std::vector<Point<dim>> x_n_p1_ref_vec = {x_n_p1_ref};
                          fe_point_evaluation.reinit(cell, x_n_p1_ref_vec);
                          Point<dim> x_n_p1_real =
                            fe_point_evaluation.quadrature_point(0);

                          // Update the solution.
                          x_n_ref  = x_n_p1_ref;
                          x_n_real = x_n_p1_real;

                          newton_it += 1;
                        } // End of the Newton solver.

                      // Compute the distance approximation: distance(x_I) =
                      // distance(x_n) + |x_n - x_I|
                      const Tensor<1, dim> x_n_to_x_I_real =
                        x_I_real - x_n_real;
                      fe_point_evaluation.evaluate(cell_dof_values,
                                                   EvaluationFlags::values);

                      double distance_value_at_x_n =
                        fe_point_evaluation.get_value(0);

                      double approx_distance =
                        compute_distance(x_n_to_x_I_real,
                                         distance_value_at_x_n);

                      // If the new distance is smaller than the previous,
                      // update the value and flag the change
                      if (distance(dof_indices[i]) > (approx_distance + 1e-8))
                        {
                          change                   = true;
                          distance(dof_indices[i]) = approx_distance;
                        }
                    } // End of the loop on the opposite faces
                }     // End of the loop on the dofs
            }
        } // End of the loop on the cells

      exchange_distance();

      // Track the change flag across the processes
      change = Utilities::MPI::logical_or(change, mpi_communicator);

      count += 1;
    } // End of the iterative while loop

  // Update the hagging node values
  constraints.distribute(distance);
}

template <int dim, typename VectorType>
void
InterfaceTools::SignedDistanceSolver<dim, VectorType>::
  compute_signed_distance_from_distance()
{
  for (auto p : this->locally_active_dofs)
    {
      signed_distance(p) = distance(p) * sgn(signed_distance_with_ghost(p));
    }

  // Update local ghost (signed_distance becomes read only)
  signed_distance.update_ghost_values();

  /* Copy distance to signed_distance_with_ghost to keep the knowledge of
  local ghost values and to have a read-only version of the vector. Ghost
  values of the signed distance are needed for volume computations.*/
  signed_distance_with_ghost = signed_distance;

  /* Zero out ghost DOFs to regain write functionalities in signed_distance
  (it becomes write-only, that is why we need signed_distance_with_ghost -
  to read the ghost values in it).*/
  signed_distance.zero_out_ghost_values();

  /* Copy the ghost values back in signed_distance (zero_out_ghost_values()
  puts zeros in ghost DOFs)*/
  for (auto p : this->locally_active_dofs)
    {
      signed_distance(p) = signed_distance_with_ghost(p);
    }
}

template class InterfaceTools::SignedDistanceSolver<2, GlobalVectorType>;
template class InterfaceTools::SignedDistanceSolver<3, GlobalVectorType>;

#ifndef LETHE_USE_LDV
template class InterfaceTools::
  SignedDistanceSolver<2, LinearAlgebra::distributed::Vector<double>>;
template class InterfaceTools::
  SignedDistanceSolver<3, LinearAlgebra::distributed::Vector<double>>;
#endif
