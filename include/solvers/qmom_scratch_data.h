// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file qmom_scratch_data.h
 * @brief Scratch and copy data used to assemble the discontinuous Galerkin
 * transport of the QMOM components.
 *
 * All the QMOM components (scalars and moments) are discretized with the same
 * scalar discontinuous element and advected by the same fluid velocity. The
 * matrix of the transport is therefore the same for every component, while
 * each component has its own right-hand side.
 */

#ifndef lethe_qmom_scratch_data_h
#define lethe_qmom_scratch_data_h

#include <core/vector.h>

#include <solvers/physics_scratch_data.h>
#include <solvers/qmom_field_function.h>

#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <memory>
#include <vector>

using namespace dealii;

/**
 * @brief Local contributions of one cell to the transport matrix, shared by
 * all the QMOM components, and to the right-hand side of every component. The
 * contributions of the interior faces to the matrix are stored separately,
 * since they couple the degrees of freedom of two cells.
 */
struct QMOMCopyData
{
  /**
   * @brief Allocate the cell contributions.
   *
   * @param[in] n_dofs Number of degrees of freedom of a cell.
   * @param[in] n_components Number of QMOM components (scalars and moments).
   */
  QMOMCopyData(const unsigned int n_dofs, const unsigned int n_components)
    : local_matrix(n_dofs, n_dofs)
    , local_rhs(n_components, Vector<double>(n_dofs))
    , local_dof_indices(n_dofs)
  {}

  /// Contributions of the cell and of its boundary faces to the matrix
  FullMatrix<double> local_matrix;

  /// Contributions of the cell and of its boundary faces to the right-hand
  /// side of every component: local_rhs[component]
  std::vector<Vector<double>> local_rhs;

  /// Global indices of the degrees of freedom of the cell
  std::vector<types::global_dof_index> local_dof_indices;

  /// Whether the cell is locally owned, i.e. whether local_matrix and
  /// local_rhs are valid
  bool cell_is_local = false;

  /// Contributions of an interior face to the degrees of freedom of the two
  /// cells that share it
  struct FaceData
  {
    /// Face contributions to the matrix, ordered as the interface degrees of
    /// freedom
    FullMatrix<double> face_matrix;

    /// Global indices of the interface degrees of freedom
    std::vector<types::global_dof_index> joint_dof_indices;
  };

  /// Contributions of the interior faces assembled with the cell
  std::vector<FaceData> face_data;
};

/**
 * @brief Scratch data of the transport of the QMOM components.
 *
 * The transport is linear and is solved directly for the solution at the end
 * of the time step (linear solver strategy). Unlike the tracer physics, which
 * assembles the residual of a Newton iteration, the assembly therefore never
 * needs the present solution: the matrix only needs the shape functions and
 * the fluid velocity, and the right-hand side the previous solutions of the
 * BDF scheme, the source term and the Dirichlet values.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class QMOMScratchData : public PhysicsScratchDataBase
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] fe_qmom Scalar discontinuous element shared by all the QMOM
   * components.
   * @param[in] quadrature Cell quadrature.
   * @param[in] face_quadrature Face quadrature.
   * @param[in] mapping Mapping of the QMOM physics, also used to evaluate the
   * fluid velocity at the same quadrature points.
   * @param[in] fe_fluid FESystem of the fluid dynamics (velocity and
   * pressure).
   * @param[in] n_components Number of QMOM components (scalars and moments).
   */
  QMOMScratchData(const FiniteElement<dim>  &fe_qmom,
                  const Quadrature<dim>     &quadrature,
                  const Quadrature<dim - 1> &face_quadrature,
                  const Mapping<dim>        &mapping,
                  const FiniteElement<dim>  &fe_fluid,
                  const unsigned int         n_components);

  /**
   * @brief Copy constructor, required by MeshWorker::mesh_loop. The FEValues
   * objects cannot be copied and are rebuilt.
   *
   * @param[in] sd Scratch data to copy.
   */
  QMOMScratchData(const QMOMScratchData<dim> &sd);

  /**
   * @brief Allocate the tables of shape functions and the vectors of values at
   * the quadrature points.
   */
  void
  allocate() override;

  /**
   * @brief Evaluate the shape functions, their gradients, the quadrature
   * points and the JxW values of a cell.
   *
   * @param[in] cell Cell of the QMOM DoFHandler.
   */
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell);

  /**
   * @brief Evaluate the previous solutions of every component at the
   * quadrature points of the cell of the last call to reinit().
   *
   * @param[in] previous_solutions Ghosted previous solutions:
   * previous_solutions[component][p] is the solution p + 1 time steps ago.
   */
  void
  reinit_previous_solutions(
    const std::vector<std::vector<GlobalVectorType>> &previous_solutions);

  /**
   * @brief Evaluate the source term of every component at the quadrature
   * points of the cell of the last call to reinit().
   *
   * @param[in] source_function Source term of all the components, at the
   * current time.
   */
  void
  reinit_source(const QMOMFieldFunction<dim> &source_function);

  /**
   * @brief Gather the fluid velocity at the quadrature points of a cell.
   *
   * @tparam VectorType Vector type of the fluid dynamics solution.
   * @param[in] fluid_cell Cell of the fluid dynamics DoFHandler.
   * @param[in] fluid_solution Present (ghosted) fluid dynamics solution. The
   * fluid dynamics is solved before the QMOM physics, so this is the velocity
   * at the end of the time step.
   */
  template <typename VectorType>
  void
  reinit_velocity(
    const typename DoFHandler<dim>::active_cell_iterator &fluid_cell,
    const VectorType                                     &fluid_solution);

  /**
   * @brief Evaluate the shape functions of the two cells that share an
   * interior face (the "here" and "there" sides), and their jump, at the face
   * quadrature points.
   *
   * @param[in] cell Cell on the "here" side of the face.
   * @param[in] face_no Face number in cell.
   * @param[in] sub_face_no Subface number in cell, or
   * numbers::invalid_unsigned_int if the face is not refined from this side.
   * @param[in] neighbor_cell Cell on the "there" side of the face.
   * @param[in] neighbor_face_no Face number in neighbor_cell.
   * @param[in] neighbor_sub_face_no Subface number in neighbor_cell.
   */
  void
  reinit_internal_face(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const unsigned int                                    face_no,
    const unsigned int                                    sub_face_no,
    const typename DoFHandler<dim>::active_cell_iterator &neighbor_cell,
    const unsigned int                                    neighbor_face_no,
    const unsigned int                                    neighbor_sub_face_no);

  /**
   * @brief Reinitialize the interface values on a boundary face.
   *
   * @param[in] cell Cell of the QMOM DoFHandler.
   * @param[in] face_no Face number in cell.
   */
  void
  reinit_boundary_face(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const unsigned int                                    face_no);

  /**
   * @brief Gather the fluid velocity at the quadrature points of a face, as
   * seen from the "here" side of the face. These quadrature points are the
   * ones of the last call to reinit_internal_face() or
   * reinit_boundary_face().
   *
   * @tparam VectorType Vector type of the fluid dynamics solution.
   * @param[in] fluid_cell Cell of the fluid dynamics DoFHandler on the "here"
   * side of the face.
   * @param[in] face_no Face number in fluid_cell.
   * @param[in] sub_face_no Subface number in fluid_cell, or
   * numbers::invalid_unsigned_int if the face is not refined from this side.
   * @param[in] fluid_solution Present (ghosted) fluid dynamics solution.
   */
  template <typename VectorType>
  void
  reinit_face_velocity(
    const typename DoFHandler<dim>::active_cell_iterator &fluid_cell,
    const unsigned int                                    face_no,
    const unsigned int                                    sub_face_no,
    const VectorType                                     &fluid_solution);

  /// FEValues of the QMOM element on cells
  FEValues<dim> fe_values;

  /// FEInterfaceValues of the QMOM element on interior and boundary faces
  FEInterfaceValues<dim> fe_interface_values;

  /// FEValues of the fluid dynamics on cells
  FEValues<dim> fe_values_fluid;

  /// FEFaceValues of the fluid dynamics on faces
  FEFaceValues<dim> fe_face_values_fluid;

  /// FESubfaceValues of the fluid dynamics on refined faces. Only allocated
  /// for hypercube meshes, since simplex meshes have no hanging nodes.
  std::unique_ptr<FESubfaceValues<dim>> fe_subface_values_fluid;

  /// Extractor of the velocity in the fluid dynamics FESystem
  const FEValuesExtractors::Vector velocities;

  /// Number of QMOM components (scalars and moments)
  const unsigned int n_components;

  /// Number of degrees of freedom of a cell
  unsigned int n_dofs;

  /// Number of cell quadrature points
  unsigned int n_q_points;

  /// Number of face quadrature points
  unsigned int n_face_q_points;

  /// Number of degrees of freedom of the current interior face (both cells)
  unsigned int n_interface_dofs;

  /// JxW values at the cell quadrature points
  std::vector<double> JxW;

  /// Cell quadrature points, in real space
  std::vector<Point<dim>> quadrature_points;

  /// Shape functions at the cell quadrature points: phi[q][i]
  Table<2, double> phi;

  /// Gradients of the shape functions at the cell quadrature points:
  /// grad_phi[q][i]
  Table<2, Tensor<1, dim>> grad_phi;

  /// Previous solutions at the cell quadrature points:
  /// previous_values[component][p][q] is the solution p + 1 time steps ago
  std::vector<std::vector<std::vector<double>>> previous_values;

  /// Source term at the cell quadrature points: source_values[q][component]
  std::vector<Vector<double>> source_values;

  /// Fluid velocity at the cell quadrature points
  std::vector<Tensor<1, dim>> velocity_values;

  /// Shape functions of the "here" cell at the face quadrature points:
  /// face_phi_here[q][i], with i an interface degree of freedom
  Table<2, double> face_phi_here;

  /// Shape functions of the "there" cell at the face quadrature points:
  /// face_phi_there[q][i], with i an interface degree of freedom
  Table<2, double> face_phi_there;

  /// Jump of the shape functions, here minus there, at the face quadrature
  /// points: face_jump_phi[q][i], with i an interface degree of freedom
  Table<2, double> face_jump_phi;

  /// Fluid velocity at the face quadrature points
  std::vector<Tensor<1, dim>> face_velocity_values;
};

#endif
