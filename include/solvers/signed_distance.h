// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_signed_distance_h
#  define lethe_signed_distance_h

#  include <deal.II/distributed/tria.h>

#  include <deal.II/dofs/dof_handler.h>

#  include <deal.II/fe/fe_q.h>
#  include <deal.II/fe/mapping_q.h>

using namespace dealii;

template <int dim, typename VectorType>
class SignedDistanceSolver
{
public:
  /**
   * @brief Base constructor.
   *
   * @param[in] background_triangulation Triangulation of the domain
   *
   * @param[in] background_fe Finite element discretizing the domain
   *
   * @param[in] p_max_distance Maximum reinitialization distance value
   *
   */
  SignedDistanceSolver(
    const parallel::DistributedTriangulationBase<dim> &background_triangulation,
    const FE_Q<dim>                                   &background_fe,
    const double                                       p_max_distance)
    : dof_handler(background_triangulation)
    , fe(background_fe)
    , mapping(fe.degree)
    , max_distance(p_max_distance)
    , pcout(std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
  {}

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
   * @param[in] background_dof_handler DoFHandler corresponding to the level-set
   * field solver
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
  VectorType
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
   * enable the auxiliary physics to output their solution via the core solver.
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
  initialize_local_distance();

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
   * @param[in] cell_transformation_jac transformation jacobian of the cell (dim
   * x dim)
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
   * reference cell. This is required because the distance minimization problem
   * is resolved in the reference face space (dim-1).
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
   * problem for the DoF x_I in the reference face space (dim - 1). In the real
   * space the residual correspond to:
   *   R = grad(d) - (x_I - x_n)/||x_I - x_n||
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
   * minimization problem for the DoF x_I in the reference face space (dim - 1).
   *
   * @param[in] stencil_real stencil in the real space. The evaluation point x_n
   * correspon to the first entry.
   *
   * @param[in] x_I_real coordinate of the DoF x_I in the real space
   *
   * @param[in] distance_gradients vector storing the distance gradients at each
   * stencil points
   *
   * @param[in] transformation_jacobians vector storing the face transformation
   * jacobians at each stencil points
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
  FE_Q<dim> fe;

  /// Mapping between the real and reference space
  MappingQ<dim> mapping;

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

  /// Solution vector of the signed distance with ghost values (read-only vesion
  /// of signed_distance)
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

  /// Surface vertices of the interface recontruction stored in a cell-wise map
  /// (volume cell)
  std::map<types::global_cell_index, std::vector<Point<dim>>>
    interface_reconstruction_vertices;
  /// Surface cells of the interface recontruction stored in a cell-wise map
  /// (volume cell)
  std::map<types::global_cell_index, std::vector<CellData<dim - 1>>>
    interface_reconstruction_cells;

  /// Set of DoFs belonging to intersected cells
  std::set<types::global_dof_index> intersected_dofs;
}
