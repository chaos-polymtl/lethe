// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file qmom_scratch_data.h
 * @brief Scratch and copy data used to evaluate the discontinuous Galerkin
 * transport operator of the QMOM physics.
 *
 * The scratch data gathers, at the cell and face quadrature points, the QMOM
 * fields (ethanol concentration, free lipids concentration and moments) and
 * the fluid velocity that advects them. The velocity is interpolated linearly
 * in time between the previous and the present fluid solutions, so that the
 * transport can be sub-cycled within a time step of the fluid dynamics.
 */

#ifndef lethe_qmom_scratch_data_h
#define lethe_qmom_scratch_data_h

#include <core/vector.h>

#include <solvers/physics_scratch_data.h>

#include <deal.II/base/quadrature.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/vector.h>

#include <memory>
#include <vector>

using namespace dealii;

/**
 * @brief Local contributions of one cell to the transport residual of the
 * QMOM physics. The contributions of the interior faces are stored
 * separately, since they couple the degrees of freedom of two cells.
 */
struct QMOMCopyData
{
  /**
   * @brief Allocate the cell contributions.
   *
   * @param[in] n_dofs Number of degrees of freedom of a cell.
   */
  explicit QMOMCopyData(const unsigned int n_dofs)
    : local_rhs(n_dofs)
    , local_dof_indices(n_dofs)
  {}

  /// Contributions of the cell and of its boundary faces
  Vector<double> local_rhs;

  /// Global indices of the degrees of freedom of the cell
  std::vector<types::global_dof_index> local_dof_indices;

  /// Whether the cell is locally owned, i.e. whether local_rhs is valid
  bool cell_is_local = false;

  /// Contributions of an interior face to the degrees of freedom of the two
  /// cells that share it
  struct FaceData
  {
    /// Face contributions, ordered as the interface degrees of freedom
    Vector<double> face_rhs;

    /// Global indices of the interface degrees of freedom
    std::vector<types::global_dof_index> joint_dof_indices;
  };

  /// Contributions of the interior faces assembled with the cell
  std::vector<FaceData> face_data;
};

/**
 * @brief Scratch data of the transport operator of the QMOM physics.
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
   * @param[in] fe_qmom FESystem of the QMOM fields.
   * @param[in] quadrature Cell quadrature.
   * @param[in] face_quadrature Face quadrature.
   * @param[in] mapping Mapping of the QMOM physics, also used to evaluate the
   * fluid velocity at the same quadrature points.
   * @param[in] fe_fluid FESystem of the fluid dynamics (velocity and
   * pressure).
   */
  QMOMScratchData(const FiniteElement<dim>  &fe_qmom,
                  const Quadrature<dim>     &quadrature,
                  const Quadrature<dim - 1> &face_quadrature,
                  const Mapping<dim>        &mapping,
                  const FiniteElement<dim>  &fe_fluid);

  /**
   * @brief Copy constructor, required by MeshWorker::mesh_loop. The FEValues
   * objects cannot be copied and are rebuilt.
   *
   * @param[in] sd Scratch data to copy.
   */
  QMOMScratchData(const QMOMScratchData<dim> &sd);

  /**
   * @brief Allocate the vectors of values at the quadrature points.
   */
  void
  allocate() override;

  /**
   * @brief Gather the QMOM fields at the quadrature points of a cell.
   *
   * @param[in] cell Cell of the QMOM DoFHandler.
   * @param[in] solution Ghosted QMOM solution.
   */
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const GlobalVectorType                               &solution);

  /**
   * @brief Gather the fluid velocity at the quadrature points of a cell.
   *
   * @tparam VectorType Vector type of the fluid dynamics solution.
   * @param[in] fluid_cell Cell of the fluid dynamics DoFHandler.
   * @param[in] previous_fluid_solution Fluid solution at the beginning of the
   * time step.
   * @param[in] present_fluid_solution Fluid solution at the end of the time
   * step.
   * @param[in] interpolation_factor Position in the time step, between 0
   * (previous solution) and 1 (present solution).
   */
  template <typename VectorType>
  void
  reinit_velocity(
    const typename DoFHandler<dim>::active_cell_iterator &fluid_cell,
    const VectorType &previous_fluid_solution,
    const VectorType &present_fluid_solution,
    const double      interpolation_factor);

  /**
   * @brief Gather the QMOM fields on both sides of an interior face.
   *
   * @param[in] cell Cell on the first side of the face.
   * @param[in] face_no Face number in cell.
   * @param[in] sub_face_no Subface number in cell, or
   * numbers::invalid_unsigned_int if the face is not refined from this side.
   * @param[in] neighbor_cell Cell on the second side of the face.
   * @param[in] neighbor_face_no Face number in neighbor_cell.
   * @param[in] neighbor_sub_face_no Subface number in neighbor_cell.
   * @param[in] solution Ghosted QMOM solution.
   */
  void
  reinit_internal_face(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const unsigned int                                    face_no,
    const unsigned int                                    sub_face_no,
    const typename DoFHandler<dim>::active_cell_iterator &neighbor_cell,
    const unsigned int                                    neighbor_face_no,
    const unsigned int                                    neighbor_sub_face_no,
    const GlobalVectorType                               &solution);

  /**
   * @brief Gather the QMOM fields on a boundary face.
   *
   * @param[in] cell Cell of the QMOM DoFHandler.
   * @param[in] face_no Face number in cell.
   * @param[in] solution Ghosted QMOM solution.
   */
  void
  reinit_boundary_face(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const unsigned int                                    face_no,
    const GlobalVectorType                               &solution);

  /**
   * @brief Gather the fluid velocity at the quadrature points of a face, as
   * seen from the first side of the face. These quadrature points are the
   * ones of the last call to reinit_internal_face or reinit_boundary_face.
   *
   * @tparam VectorType Vector type of the fluid dynamics solution.
   * @param[in] fluid_cell Cell of the fluid dynamics DoFHandler on the first
   * side of the face.
   * @param[in] face_no Face number in fluid_cell.
   * @param[in] sub_face_no Subface number in fluid_cell, or
   * numbers::invalid_unsigned_int if the face is not refined from this side.
   * @param[in] previous_fluid_solution Fluid solution at the beginning of the
   * time step.
   * @param[in] present_fluid_solution Fluid solution at the end of the time
   * step.
   * @param[in] interpolation_factor Position in the time step, between 0
   * (previous solution) and 1 (present solution).
   */
  template <typename VectorType>
  void
  reinit_face_velocity(
    const typename DoFHandler<dim>::active_cell_iterator &fluid_cell,
    const unsigned int                                    face_no,
    const unsigned int                                    sub_face_no,
    const VectorType &previous_fluid_solution,
    const VectorType &present_fluid_solution,
    const double      interpolation_factor);

  /// FEValues of the QMOM fields on cells
  FEValues<dim> fe_values;

  /// FEInterfaceValues of the QMOM fields on interior and boundary faces
  FEInterfaceValues<dim> fe_interface_values;

  /// FEValues of the fluid dynamics on cells
  FEValues<dim> fe_values_fluid;

  /// FEFaceValues of the fluid dynamics on faces
  FEFaceValues<dim> fe_face_values_fluid;

  /// FESubfaceValues of the fluid dynamics on refined faces. Only allocated
  /// for hypercube meshes, since simplex meshes have no hanging nodes.
  std::unique_ptr<FESubfaceValues<dim>> fe_subface_values_fluid;

  /// Extractor of the velocity in the fluid dynamics FESystem
  const FEValuesExtractors::Vector velocity_extractor;

  /// Number of components of the QMOM FESystem
  unsigned int n_components;

  /// Number of degrees of freedom of a cell
  unsigned int n_dofs;

  /// Number of cell quadrature points
  unsigned int n_q_points;

  /// Number of face quadrature points
  unsigned int n_face_q_points;

  /// Component of each degree of freedom of a cell
  std::vector<unsigned int> component_of_dof;

  /// QMOM fields at the cell quadrature points, indexed by [q](component)
  std::vector<Vector<double>> values;

  /// Fluid velocity at the cell quadrature points
  std::vector<Tensor<1, dim>> velocity_values;

  /// Number of interface degrees of freedom of the current face
  unsigned int n_interface_dofs;

  /// Component of each interface degree of freedom of the current face
  std::vector<unsigned int> component_of_interface_dof;

  /// QMOM fields on the first and second sides of the current face, indexed
  /// by [q](component). values_there is not used on boundary faces.
  std::vector<Vector<double>> values_here;
  std::vector<Vector<double>> values_there;

  /// Fluid velocity at the face quadrature points
  std::vector<Tensor<1, dim>> face_velocity_values;

  /// Upwind values of the QMOM fields at a boundary quadrature point, used
  /// as a buffer when assembling the boundary fluxes
  Vector<double> upwind_values;

private:
  /// Fluid velocity at the beginning of the time step, used as a temporary
  /// when interpolating the velocity in time
  std::vector<Tensor<1, dim>> previous_velocity_values;

  /// Fluid velocity at the face quadrature points at the beginning of the
  /// time step, used as a temporary when interpolating the velocity in time
  std::vector<Tensor<1, dim>> previous_face_velocity_values;
};

#endif
