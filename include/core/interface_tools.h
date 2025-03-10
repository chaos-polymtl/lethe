// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_interface_tools_h
#define lethe_interface_tools_h

#include <core/utilities.h>
#include <core/vector.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

namespace InterfaceTools
{
  /**
   * @brief Scalar function defined by the DOF values of a single cell. Based on
   * the CellWiseFunction and RefSpaceFEFieldFunction of dealii.
   *
   * @tparam dim An integer that denotes the dimension of the space in which
   * the problem is solved.
   *
   * @tparam VectorType The vector type of the solution vector.
   *
   * @tparam FEType The finite element type used to discretize the problem.
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
     * @param[in] p_fe_degree Finite element degree discretizing the field we
     * want to convert to a CellWiseFunction.
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
     * @brief Return the function gradient of the specified component at the
     * given point in the reference cell.
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

    /// Value of the level-set field at the DoFs of the cell
    VectorType cell_dof_values;
  };

  template <int dim, typename VectorType, typename FEType>
  CellWiseFunction<dim, VectorType, FEType>::CellWiseFunction(
    const unsigned int p_fe_degree)
    : fe(p_fe_degree)
  {
    n_cell_wise_dofs = fe.dofs_per_cell;
  }

  template <int dim, typename VectorType, typename FEType>
  inline void
  CellWiseFunction<dim, VectorType, FEType>::set_active_cell(
    const VectorType &in_local_dof_values)
  {
    cell_dof_values = in_local_dof_values;
  }

  template <int dim, typename VectorType, typename FEType>
  inline double
  CellWiseFunction<dim, VectorType, FEType>::value(
    const Point<dim>  &point,
    const unsigned int component) const
  {
    double value = 0;
    for (unsigned int i = 0; i < n_cell_wise_dofs; ++i)
      value +=
        cell_dof_values[i] * fe.shape_value_component(i, point, component);

    return value;
  }

  template <int dim, typename VectorType, typename FEType>
  inline Tensor<1, dim>
  CellWiseFunction<dim, VectorType, FEType>::gradient(
    const Point<dim>  &point,
    const unsigned int component) const
  {
    Tensor<1, dim> gradient;
    for (unsigned int i = 0; i < n_cell_wise_dofs; ++i)
      gradient +=
        cell_dof_values[i] * fe.shape_grad_component(i, point, component);

    return gradient;
  }

  template <int dim, typename VectorType, typename FEType>
  inline SymmetricTensor<2, dim>
  CellWiseFunction<dim, VectorType, FEType>::hessian(
    const Point<dim>  &point,
    const unsigned int component) const
  {
    Tensor<2, dim> hessian;
    for (unsigned int i = 0; i < n_cell_wise_dofs; ++i)
      hessian +=
        cell_dof_values[i] * fe.shape_grad_grad_component(i, point, component);

    return symmetrize(hessian);
  }

  /**
   * @brief
   * Compute the volume enclosed by the 0 level of a level-set field
   * inside a cell. The inside volume is computed, defined by negative values of
   * the level-set field.
   *
   * @tparam dim An integer that denotes the dimension of the space in which
   * the problem is solved.
   *
   * @param[in] fe_point_evaluation FePointEvaluation
   *
   * @param[in] cell Cell for which the volume is computed
   *
   * @param[in] cell_dof_values Cell DOFs value of the level set field
   *
   * @param[in] corr Correction to apply to the DOF values (constant for all
   * DOFs). It can be used if we want the volume enclosed by the iso-level equal
   * to the calue of the correction.
   *
   * @param[in] n_quad_points Number of quadrature points for the volume
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
   * Compute the volume enclosed by a given level of a level set field
   * in the domain. The inside volume is computed, defined by negative values of
   * the level-set field.
   *
   * @tparam dim An integer that denotes the dimension of the space in which
   * the problem is solved.
   *
   * @tparam VectorType The vector type of the solution vector.
   *
   * @param[in] mapping Mapping of the domain
   *
   * @param[in] dof_handler DofHandler associated to the triangulation on which
   * the volume is computed
   *
   * @param[in] level_set_vector Level-set vector
   *
   * @param[in] iso_level Given level of the level-set field enclosing the
   * volume of interest
   *
   * @param[in] mpi_communicator MPI communicator
   *
   * @return Volume enclosed by the specified level
   */
  template <int dim, typename VectorType>
  double
  compute_volume(const Mapping<dim>       &mapping,
                 const DoFHandler<dim>    &dof_handler,
                 const FiniteElement<dim> &fe,
                 const VectorType         &level_set_vector,
                 const double              iso_level,
                 const MPI_Comm           &mpi_communicator);

  /**
   * @brief
   * Reconstruct the interface defined by a given level of a level set field
   * in the domain.
   *
   * @tparam dim An integer that denotes the dimension of the space in which
   * the problem is solved.
   *
   * @tparam VectorType The vector type of the solution vector.
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
   * @param[in] iso_level Given level of the level-set field defining the
   * interface of interest
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
    const double              iso_level,
    std::map<types::global_cell_index, std::vector<Point<dim>>>
      &interface_reconstruction_vertices,
    std::map<types::global_cell_index, std::vector<CellData<dim - 1>>>
                                      &interface_reconstruction_cells,
    std::set<types::global_dof_index> &intersected_dofs);


  template <int dim, typename VectorType>
  class SignedDistanceSolver
  {
  public:
    /**
     * @brief Base constructor.
     *
     * @param[in] background_triangulation Shared pointer to the triangulation
     * of the domain
     *
     * @param[in] background_fe Shared pointer to the dinite element
     * discretizing the domain
     *
     * @param[in] p_max_distance Maximum reinitialization distance value
     *
     */
    SignedDistanceSolver(
      std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                          background_triangulation,
      std::shared_ptr<FiniteElement<dim>> background_fe,
      const double                        p_max_distance)
      : dof_handler(*background_triangulation)
      , fe(background_fe)
      , max_distance(p_max_distance)
      , pcout(std::cout,
              (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
    {
      mapping = std::make_shared<MappingFE<dim>>(*fe);
    }

    /**
     * @brief setup_dofs
     *
     * Initialize the degree of freedom and the memory
     * associated with them for fluid dynamics and enabled auxiliary physics.
     */
    void
    setup_dofs(const MPI_Comm &mpi_communicator);

    /**
     * @brief set_level_set_from_background_mesh
     *
     * Set the level-set field from the main solver. For example, when using the
     * current SignedDistanceSolver for the geometric redistanciation of the VOF
     * phase fraction, the level-set field comes from the VOF solver and is
     * described by the corresponding DoFHandler.
     *
     * @param[in] background_dof_handler DoFHandler corresponding to the
     * level-set field solver
     *
     * @param[in] background_level_set_vector level-set solution vector
     *
     * @param[in] mpi_communicator MPI communicator
     */
    void
    set_level_set_from_background_mesh(
      const DoFHandler<dim> &background_dof_handler,
      const VectorType      &background_level_set_vector,
      const MPI_Comm        &mpi_communicator);

    /**
     * @brief solve
     *
     * Solve for the signed distance from the 0-level of the level-set vector.
     */
    void
    solve(const MPI_Comm &mpi_communicator);

    /**
     * @brief Getter methods to get the private attribute signed_distance
     *
     * @param[in] mpi_communicator MPI communicator
     *
     * @return vector storing the computed signed distance
     */
    VectorType &
    get_signed_distance(const MPI_Comm &mpi_communicator);

    /**
     * @brief Output the interface reconstruction used for the signed distance
     * computations
     *
     * @param[in] output_name name of the output file
     *
     * @param[in] output_path path to the output file
     *
     * @param[in] it iteration number
     *
     */
    void
    output_interface_reconstruction(const std::string  output_name,
                                    const std::string  output_path,
                                    const unsigned int it) const;

    /**
     * @brief Attach the solution vector to the DataOut provided. This function
     * enable the auxiliary physics to output their solution via the core
     * solver.
     *
     * @param[in,out] data_out DataOut reponsible for solution output
     */
    void
    attach_solution_to_output(DataOut<dim> &data_out);

    /// DoFHandler describing the problem
    DoFHandler<dim> dof_handler;

  private:
    /**
     * @brief Zero the ghost dofs entries of the solution vectors to gain write
     * access to the ghost elements
     */
    void
    zero_out_ghost_values();

    /**
     * @brief Update the ghost dofs entries of the solution vectors to gain read
     * access to the ghost elements
     */
    void
    update_ghost_values();

    /**
     * @brief Exchange the distance and set it to the minimum value across the
     * processors by using the compress(VectorOperation::min) function
     */
    void
    exchange_distance();

    /**
     * @brief Initialize the distance solution vectors to the given max_distance
     */
    void
    initialize_distance();

    /**
     * @brief Compute the geometric distance (brute force) between the interface
     * reconstruction and the dofs of the intersected cells (first neighbors)
     */
    void
    compute_first_neighbors_distance();

    /**
     * @brief Compute the geometric distance (marching method) between the
     * interface reconstruction and the rest of the dofs (second neighbors)
     *
     * @param[in] mpi_communicator MPI communicator
     */
    void
    compute_second_neighbors_distance(const MPI_Comm &mpi_communicator);

    /**
     * @brief Compute the signed_distance from the distance and the sign of the
     * signed_distance entries
     */
    void
    compute_signed_distance_from_distance();

    /**
     * @brief Compute the cell-wise volume correction to match the volume englobed
     * by the level 0 of the level_set field in each element
     */
    void
    compute_cell_wise_volume_correction();

    /**
     * @brief Correct the global volume to match the volume englobed by the level
     * 0 of the level_set field
     *
     * @param[in] mpi_communicator MPI communicator
     */
    void
    conserve_global_volume(const MPI_Comm &mpi_communicator);

    /**
     * @brief Return the local id of the opposite faces to the given local DOF
     * (works for quad only).
     *
     * @param[in] local_dof_id Local id of the DOF
     *
     * @param[out] local_opposite_faces The vector containing the id of the
     * opposite faces
     */
    inline void
    get_dof_opposite_faces(unsigned int               local_dof_id,
                           std::vector<unsigned int> &local_opposite_faces);

    /**
     * @brief Return the face transformation jacobian (dim-1 x dim-1). This is
     * required because the distance minimization problem is resolved in the
     * reference face space (dim-1).
     *
     * @param[in] cell_transformation_jac transformation jacobian of the cell
     * (dim x dim)
     *
     * @param[in] local_face_id local id of the face
     *
     * @param[out] face_transformation_jac face transformation jacobian (dim-1 x
     * dim-1) faces
     */
    inline void
    get_face_transformation_jacobian(
      const DerivativeForm<1, dim, dim> &cell_transformation_jac,
      const unsigned int                 local_face_id,
      DerivativeForm<1, dim - 1, dim>   &face_transformation_jac);

    /**
     * @brief
     * Transform a point dim-1 in a reference face to a point dim in the
     * reference cell. This is required because the distance minimization
     * problem is resolved in the reference face space (dim-1).
     *
     * @param[in] x_ref_face point dim-1 in the reference face
     *
     * @param[in] local_face_id local id of the face
     *
     * @return Point dim in the reference cell
     */
    inline Point<dim>
    transform_ref_face_point_to_ref_cell(const Point<dim - 1> &x_ref_face,
                                         const unsigned int    local_face_id);

    /**
     * @brief Compute the residual at the point x_n of the distance minimization
     * problem for the DoF x_I in the reference face space (dim - 1). In the
     * real space the residual correspond to: R = grad(d) - (x_I - x_n)/||x_I -
     * x_n||
     *
     *
     * @param[in] x_n_to_x_I_real vector from x_n to x_I in the real space
     *
     * @param[in] distance_gradient gradient of the distance at the point x_n
     *
     * @param[in] transformation_jac transformation jacobian of the face
     *
     * @param[out] residual_ref residual in the reference face space
     */
    inline void
    compute_residual(const Tensor<1, dim>                  &x_n_to_x_I_real,
                     const Tensor<1, dim>                  &distance_gradient,
                     const DerivativeForm<1, dim - 1, dim> &transformation_jac,
                     Tensor<1, dim - 1>                    &residual_ref);

    /**
     * @brief Compute the stencil for the numerical jacobian of the distance
     * minimization problem in the face where the problem is solved
     *
     * @param[in] x_ref point in the reference cell space where the jacobian is
     * evaluated
     *
     * @param[in] local_face_id local id of the face in which the numerical
     * jacobian is required
     *
     * @param[in] perturbation value of the perturbation use for the
     * numerical jacobian computation
     *
     * @return vector containing the 2*dim - 1 stencil point. The entries of the
     * vector are the following:
     *                                    4
     *
     *                               1    0    2
     *
     *                                    3
     * The entry 0 is the current evaluation point x_ref
     */
    inline std::vector<Point<dim>>
    compute_numerical_jacobian_stencil(const Point<dim>   x_ref,
                                       const unsigned int local_face_id,
                                       const double       perturbation);

    /**
     * @brief
     * Transform the Newton correction in a reference face to a tensor (dim) in
     * the reference cell. This is required because the distance minimization
     * problem is resolved in the reference face space (dim-1).
     *
     * @param[in] x_ref_face point dim-1 in the reference face
     *
     * @param[in] local_face_id local id of the face
     *
     * @return tensor dim in the reference cell
     */
    inline Tensor<1, dim>
    transform_ref_face_correction_to_ref_cell(
      const Vector<double> &correction_ref_face,
      const unsigned int    local_face_id);

    /**
     * @brief
     * Compute the numerical jacobian at the point x_n of the distance
     * minimization problem for the DoF x_I in the reference face space (dim -
     * 1).
     *
     * @param[in] stencil_real stencil in the real space. The evaluation point
     * x_n correspon to the first entry.
     *
     * @param[in] x_I_real coordinate of the DoF x_I in the real space
     *
     * @param[in] distance_gradients vector storing the distance gradients at
     * each stencil points
     *
     * @param[in] transformation_jacobians vector storing the face
     * transformation jacobians at each stencil points
     *
     * @param[in] perturbation value of the perturbation use for the
     * numerical jacobian computation
     *
     * @param[out] jacobian_matrix jacobian matrix of the minimization problem
     */
    inline void
    compute_numerical_jacobian(
      const std::vector<Point<dim>>     &stencil_real,
      const Point<dim>                  &x_I_real,
      const std::vector<Tensor<1, dim>> &distance_gradients,
      const std::vector<DerivativeForm<1, dim - 1, dim>>
                               &transformation_jacobians,
      const double              perturbation,
      LAPACKFullMatrix<double> &jacobian_matrix);

    /**
     * @brief
     * Compute the distance accordinf to: d(x_I) = d(x_n) + ||x_I - x_n||
     *
     * @param[in] x_n_to_x_I_real vector from x_n to x_I in the real space
     *
     * @param[in] distance value of the distance at the point x_n
     *
     * @return distance between x_I and the interface
     */
    inline double
    compute_distance(const Tensor<1, dim> &x_n_to_x_I_real,
                     const double          distance);

    /// Finite element discretizing the problem
    std::shared_ptr<FiniteElement<dim>> fe;

    /// Mapping between the real and reference space
    std::shared_ptr<Mapping<dim>> mapping;

    /// Maximum redistanciation distance
    const double max_distance;

    /// Parallel output stream
    ConditionalOStream pcout;

    /// Set of locally owned DoFs
    IndexSet locally_owned_dofs;

    /// Set of locally relevant DoFs
    IndexSet locally_relevant_dofs;

    /// Set of locally active DoFs
    IndexSet locally_active_dofs;

    /// Level-set field coming from the main solver
    VectorType level_set;

    /// Solution vector of the signed distance (write only)
    LinearAlgebra::distributed::Vector<double> signed_distance;

    /// Solution vector of the signed distance with ghost values (read-only
    /// vesion of signed_distance)
    LinearAlgebra::distributed::Vector<double> signed_distance_with_ghost;

    /// Solution vector of the distance (write only)
    LinearAlgebra::distributed::Vector<double> distance;

    /// Solution vector of the distance with ghost values (read-only vesion of
    /// distance)
    LinearAlgebra::distributed::Vector<double> distance_with_ghost;

    /// Value of the orrection to apply to the signed_distance to match the
    /// cell-wise volume encompassed by the level 0 of level_set
    LinearAlgebra::distributed::Vector<double> volume_correction;

    /// Hanging node contraint
    AffineConstraints<double> constraints;

    /// Surface vertices of the interface recontruction stored in a cell-wise
    /// map (volume cell)
    std::map<types::global_cell_index, std::vector<Point<dim>>>
      interface_reconstruction_vertices;
    /// Surface cells of the interface recontruction stored in a cell-wise map
    /// (volume cell)
    std::map<types::global_cell_index, std::vector<CellData<dim - 1>>>
      interface_reconstruction_cells;

    /// Set of DoFs belonging to intersected cells
    std::set<types::global_dof_index> intersected_dofs;
  };

  // template <int dim, typename VectorType>
  // void
  // SignedDistanceSolver<dim, VectorType>::setup_dofs(const MPI_Comm
  // &mpi_communicator)
  // {
  //   dof_handler.distribute_dofs(fe);
  //
  //   locally_owned_dofs = dof_handler.locally_owned_dofs();
  //   locally_relevant_dofs =
  //     DoFTools::extract_locally_relevant_dofs(dof_handler);
  //   locally_active_dofs = DoFTools::extract_locally_active_dofs(dof_handler);
  //
  //   level_set.reinit(locally_owned_dofs,
  //                    locally_relevant_dofs,
  //                    mpi_communicator);
  //
  //   signed_distance.reinit(locally_owned_dofs,
  //                          locally_active_dofs,
  //                          mpi_communicator);
  //   signed_distance_with_ghost.reinit(locally_owned_dofs,
  //                                     locally_active_dofs,
  //                                     mpi_communicator);
  //
  //   distance.reinit(locally_owned_dofs, locally_active_dofs,
  //   mpi_communicator); distance_with_ghost.reinit(locally_owned_dofs,
  //                              locally_active_dofs,
  //                              mpi_communicator);
  //
  //   volume_correction.reinit(locally_owned_dofs,
  //                            locally_active_dofs,
  //                            mpi_communicator);
  //
  //   constraints.clear();
  //   constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  //   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  //   constraints.close();
  // }
} // namespace InterfaceTools

#endif
