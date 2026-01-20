// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_interface_tools_h
#define lethe_interface_tools_h

#include <core/lethe_grid_tools.h>
#include <core/output_struct.h>
#include <core/solutions_output.h>
#include <core/utilities.h>
#include <core/vector.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/lapack_full_matrix.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/non_matching/fe_immersed_values.h>

#include <deal.II/numerics/vector_tools.h>

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
    CellWiseFunction(const unsigned int p_fe_degree)
      : fe(p_fe_degree)
      , n_cell_wise_dofs(fe.dofs_per_cell)
    {}

    /**
     * @brief Set the cell that the function should be evaluated on.
     *
     * @param[in] in_local_dof_values Cell's DOF values
     *
     */
    void
    set_active_cell(const VectorType &in_local_dof_values)
    {
      cell_dof_values = in_local_dof_values;
    }

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

    /// Number of DoFs per element
    unsigned int n_cell_wise_dofs;

    /// Value of the level-set field at the DoFs of the cell
    VectorType cell_dof_values;
  };


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
   * @param[in] cell_dof_level_set_values Cell DoF values of the level-set field
   *
   * @param[in] corr Correction to apply to the DOF values (constant for all
   * DoFs). It can be used if we want the volume enclosed by the iso-level equal
   * to the value of the correction.
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
   * Compute the surface area of a given level of a level set field and the
   * volume enclosed by it. The inside volume is computed, defined by negative
   * values of the level-set field.
   *
   * @tparam dim An integer that denotes the dimension of the space in which
   * the problem is solved.
   *
   * @tparam VectorType The vector type of the solution vector.
   *
   * @param[in] dof_handler DofHandler associated to the triangulation on which
   * the volume is computed
   *
   * @param[in] fe Finite element
   *
   * @param[in] level_set_vector Level-set vector
   *
   * @param[in] iso_level Given level of the level-set field enclosing the
   * volume of interest
   *
   * @param[in] mpi_communicator MPI communicator
   *
   * @return Surface area of the specified level and the volume enclosed by it
   */
  template <int dim, typename VectorType>
  std::pair<double, double>
  compute_surface_and_volume(const DoFHandler<dim>    &dof_handler,
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
   * @param[in,out] intersected_dofs Set of DoFs that belong to intersected
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


  /**
   * @brief Interface to build the patches of the interface reconstruction
   * vertices for visualization. It reproduce the same structure as particle
   * visualization.
   *
   * @tparam dim An integer that denotes the number of spatial dimensions.
   *
   */
  template <int dim>
  class InterfaceReconstructionDataOut : public dealii::DataOutInterface<0, dim>
  {
  public:
    /**
     * @brief Default constructor.
     */
    InterfaceReconstructionDataOut()
    {}

    /**
     * @brief Build the patches of the interface reconstruction vertices for
     * visualization.
     *
     * @param[in,out] interface_reconstruction_vertices Cell-wise map of the
     * reconstructed surface vertices. The map contains vectors storing the
     * vertices of the reconstructed surface for each intersected volume cell
     * (dim).
     */
    void
    build_patches(
      const std::map<types::global_cell_index, std::vector<Point<dim>>>
        &interface_reconstruction_vertices);

  private:
    /**
     * @brief Implementation of the corresponding function of the base class.
     */
    const std::vector<DataOutBase::Patch<0, dim>> &
    get_patches() const override
    {
      return patches;
    }

    /**
     * @brief Implementation of the corresponding function of the base class.
     */
    std::vector<std::string>
    get_dataset_names() const override
    {
      return dataset_names;
    }

    /// Output information that is filled by build_patches() and
    /// written by the write function of the base class.
    std::vector<DataOutBase::Patch<0, dim>> patches;

    /// A list of field names for all data components stored in patches.
    std::vector<std::string> dataset_names;
  };

  template <int dim>
  void
  InterfaceReconstructionDataOut<dim>::build_patches(
    const std::map<types::global_cell_index, std::vector<Point<dim>>>
      &interface_reconstruction_vertices)
  {
    for (auto const &cell : interface_reconstruction_vertices)
      {
        const std::vector<Point<dim>> &vertices = cell.second;
        for (const Point<dim> &vertex : vertices)
          {
            DataOutBase::Patch<0, dim> temp;
            temp.vertices[0] = vertex;
            patches.push_back(temp);
          }
      }
  }

  /**
   * @brief Solver to compute the signed distance from a given level of a level-set field. It is based on the method proposed in the following article: Ausas, R.F., Dari, E.A. and Buscaglia, G.C. (2011), A geometric mass-preserving redistancing scheme for the level set function. Int. J. Numer. Meth. Fluids, 65: 989-1010. https://doi.org/10.1002/fld.2227.
   *
   * @tparam dim An integer that denotes the dimension of the space in which
   * the problem is solved.
   *
   * @tparam VectorType The vector type of the level-set vector.
   */
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
     * @param[in] background_fe Shared pointer to the finite element
     * discretizing the domain
     *
     * @param[in] p_max_distance Maximum reinitialization distance value
     *
     * @param[in] p_iso_level Iso-level (before scaling) from which the signed
     * distance is computed
     *
     * @param[in] p_scaling Scaling factor to apply to the input level-set
     * field. It is used to switch the inside and outside subdomains, if needed.
     *
     * @param[in] p_verbosity Verbosity level
     *
     */
    SignedDistanceSolver(
      std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                          background_triangulation,
      std::shared_ptr<FiniteElement<dim>> background_fe,
      const double                        p_max_distance,
      const double                        p_iso_level,
      const double                        p_scaling,
      const Parameters::Verbosity         p_verbosity)
      : dof_handler(*background_triangulation)
      , fe(background_fe)
      , max_distance(p_max_distance)
      , iso_level(p_iso_level * p_scaling)
      , scaling(p_scaling)
      , verbosity(p_verbosity)
      , pcout(std::cout,
              (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
    {
      // MappingQ is required to reduce significantly the computational time of
      // the FEPointEvaluation.reinit(). Using MappingQ and FE_Q allow to take
      // the "fast path" in FEPointEvaluation calls.
      mapping = std::make_shared<MappingQ<dim>>(fe->degree);
      set_face_opposite_dofs_map();
      set_face_dofs_map();
    }

    /**
     * @brief Initialize the degrees of freedom and associated memory.
     *
     * This function sets up the DoF handler, initializes the solution vectors,
     * and configures the constraints for the signed distance solver.
     */
    void
    setup_dofs();

    /**
     * @brief Set the level-set field from the background mesh solver.
     *
     * This function transfers the level-set field from the main solver to the
     * signed distance solver. For example, when using geometric redistanciation
     * of the VOF phase fraction, the level-set field comes from the VOF solver.
     *
     * @param[in] background_dof_handler DoFHandler corresponding to the
     * level-set field solver.
     * @param[in] background_level_set_vector Level-set solution vector from
     * the background solver.
     */
    void
    set_level_set_from_background_mesh(
      const DoFHandler<dim> &background_dof_handler,
      const VectorType      &background_level_set_vector);

    /**
     * @brief Solve for the signed distance from the given level of the level-set vector.
     *
     * This function computes the signed distance field using the method
     * described in Ausas et al. (2011). The algorithm preserves the interface
     * location while computing accurate distances throughout the domain.
     */
    void
    solve();

    /**
     * @brief Get the computed signed distance field.
     *
     * @return Reference to the vector storing the computed signed distance values.
     */
    VectorType &
    get_signed_distance();

    /**
     * @brief Output the interface reconstruction used for signed distance computations.
     *
     * This function writes the reconstructed interface vertices and cells to
     * VTU files for visualization purposes.
     *
     * @param[in] output_name Name of the output file.
     * @param[in] output_path Path to the output directory.
     * @param[in] time Current simulation time.
     * @param[in] it Iteration number for file naming.
     */
    void
    output_interface_reconstruction(const std::string &output_name,
                                    const std::string &output_path,
                                    const double       time,
                                    const unsigned int it);

    /**
     * @brief Output the signed distance field for visualization.
     *
     * This function writes the computed signed distance field to VTU files
     * for post-processing and visualization.
     *
     * @param[in] output_name Name of the output file.
     * @param[in] output_path Path to the output directory.
     * @param[in] time Current simulation time.
     * @param[in] it Iteration number for file naming.
     */
    void
    output_signed_distance(const std::string &output_name,
                           const std::string &output_path,
                           const double       time,
                           const unsigned int it);

    /**
     * @brief Gather and return vector of output structs that are particular to some applications.
     *
     * @return Vector of OutputStructs that will be used to write the output results as VTU files.
     */
    std::vector<OutputStruct<dim, GlobalVectorType>>
    gather_output_hook();

    /// DoFHandler describing the problem
    DoFHandler<dim> dof_handler;

  private:
    /**
     * @brief Zero the ghost DoFs entries of the solution vectors to gain write
     * access to the ghost elements
     */
    void
    zero_out_ghost_values() const;

    /**
     * @brief Update the ghost DoFs entries of the solution vectors to gain read
     * access to the ghost elements
     */
    void
    update_ghost_values() const;

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
     * reconstruction and the DoFs of the intersected cells (first neighbors)
     */
    void
    compute_first_neighbors_distance();

    /**
     * @brief Compute the geometric distance (marching method) between the
     * interface reconstruction and the rest of the DoFs (second neighbors)
     */
    void
    compute_second_neighbors_distance();

    /**
     * @brief Compute the signed_distance from the distance vector and the sign of the
     * signed_distance vector entries
     */
    void
    compute_signed_distance_from_distance();

    /**
     * @brief Compute the cell-wise volume correction to match the volume englobed
     * by the given level of the level_set field in each element
     */
    void
    compute_cell_wise_volume_correction();

    /**
     * @brief Correct the global volume to match the volume englobed by the given level
     * of the level_set field
     */
    void
    conserve_global_volume();

    /**
     * @brief Return the local id of the opposite faces to the given local Dof
     * (works for quad only).
     *
     * @param[in] local_dof_id Local id of the DoF in the cell
     *
     * @param[out] local_opposite_faces The vector containing the id of the
     * opposite faces
     */
    inline void
    get_dof_opposite_faces(unsigned int               local_dof_id,
                           std::vector<unsigned int> &local_opposite_faces)
    {
      unsigned int local_dof_id_2d = local_dof_id % 4;

      local_opposite_faces[0] = (local_dof_id_2d + 1) % 2;
      local_opposite_faces[1] = 3 - local_dof_id_2d / 2;

      if constexpr (dim == 3)
        local_opposite_faces[2] = 5 - local_dof_id / 4;
    };

    /**
     * @brief Set the map of local id of the opposite DoFs to the given face
     * (works for quad only).
     */
    inline void
    set_face_opposite_dofs_map()
    {
      if constexpr (dim == 2)
        {
          face_opposite_dofs_map[0] = {1, 3};
          face_opposite_dofs_map[1] = {0, 2};
          face_opposite_dofs_map[2] = {2, 3};
          face_opposite_dofs_map[3] = {0, 1};
        }

      if constexpr (dim == 3)
        {
          face_opposite_dofs_map[0] = {1, 3, 5, 7};
          face_opposite_dofs_map[1] = {0, 2, 4, 6};
          face_opposite_dofs_map[2] = {2, 3, 6, 7};
          face_opposite_dofs_map[3] = {0, 1, 4, 5};
          face_opposite_dofs_map[4] = {4, 5, 6, 7};
          face_opposite_dofs_map[5] = {0, 1, 2, 3};
        }
    }

    /**
     * @brief Set the map of local id of the opposite DoFs to the given face
     * (works for quad only).
     */
    inline void
    set_face_dofs_map()
    {
      if constexpr (dim == 2)
        {
          face_dofs_map[0] = {0, 2};
          face_dofs_map[1] = {1, 3};
          face_dofs_map[2] = {0, 1};
          face_dofs_map[3] = {2, 3};
        }

      if constexpr (dim == 3)
        {
          face_dofs_map[0] = {0, 2, 4, 6};
          face_dofs_map[1] = {1, 3, 5, 7};
          face_dofs_map[2] = {0, 1, 4, 5};
          face_dofs_map[3] = {2, 3, 6, 7};
          face_dofs_map[4] = {0, 1, 2, 3};
          face_dofs_map[5] = {4, 5, 6, 7};
        }
    }

    /**
     * @brief Return the local id of the opposite DoFs to the given face
     * (works for quad only).
     *
     * @param[in] local_face_id Local id of the face in the cell
     *
     * @param[out] local_opposite_dofs The vector containing the id of the
     * opposite faces
     */
    inline void
    get_face_opposite_dofs(unsigned int               local_face_id,
                           std::vector<unsigned int> &local_opposite_dofs)
    {
      local_opposite_dofs = face_opposite_dofs_map.at(local_face_id);
    };

    /**
     * @brief Return the local id of the opposite DoFs to the given face
     * (works for quad only).
     *
     * @param[in] local_face_id Local id of the face in the cell
     *
     * @param[out] local_opposite_dofs The vector containing the id of the
     * opposite faces
     */
    inline void
    get_face_local_dofs(unsigned int               local_face_id,
                        std::vector<unsigned int> &local_dofs)
    {
      local_dofs = face_dofs_map.at(local_face_id);
    };

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
      DerivativeForm<1, dim - 1, dim>   &face_transformation_jac)
    {
      for (unsigned int i = 0; i < dim; ++i)
        {
          int k = 0;
          for (unsigned int j = 0; j < dim; ++j)
            {
              if (local_face_id / 2 == j)
                continue;

              if (k < dim - 1)
                face_transformation_jac[i][k] = cell_transformation_jac[i][j];
              k += 1;
            }
        }
    };

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
                                         const unsigned int    local_face_id)
    {
      Point<dim> x_ref_cell;

      unsigned int j = 0;
      for (unsigned int i = 0; i < dim; ++i)
        {
          if (local_face_id / 2 == i)
            {
              x_ref_cell[i] = double(local_face_id % 2);
              continue;
            }
          x_ref_cell[i] = x_ref_face[j];
          j += 1;
        }

      return x_ref_cell;
    };

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
                     Tensor<1, dim - 1>                    &residual_ref)
    {
      Tensor<1, dim> residual_real =
        distance_gradient - (1.0 / x_n_to_x_I_real.norm()) * x_n_to_x_I_real;

      DerivativeForm<1, dim, dim - 1> transformation_jac_transpose =
        transformation_jac.transpose();

      for (unsigned int i = 0; i < dim - 1; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          residual_ref[i] +=
            transformation_jac_transpose[i][j] * residual_real[j];
    };

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
    inline void
    compute_numerical_jacobian_stencil(const Point<dim>         x_ref,
                                       const unsigned int       local_face_id,
                                       const double             perturbation,
                                       std::vector<Point<dim>> &stencil)
    {
      for (unsigned int i = 0; i < 2 * dim - 1; ++i)
        {
          stencil[i] = x_ref;
        }

      unsigned int              skip_index = local_face_id / 2;
      std::vector<unsigned int> j_index(dim - 1);

      // Set the coordinates (x,y,z) to be skipped (the coordinates not
      // perturbed)
      unsigned int j = 0;
      for (unsigned int i = 0; i < dim; ++i)
        {
          if (i == skip_index)
            continue;
          j_index[j] = i;
          j += 1;
        }

      for (unsigned int i = 1; i < dim; ++i)
        {
          j = j_index[i - 1];
          stencil[2 * i - 1][j] -= perturbation;
          stencil[2 * i][j] += perturbation;
        }
    };

    /**
     * @brief
     * Transform the Newton correction in a reference face to a tensor (dim) in
     * the reference cell. This is required because the distance minimization
     * problem is resolved in the reference face space (dim-1).
     *
     * @param[in] correction_ref_face Newton correction in the reference face
     *
     * @param[in] local_face_id local id of the face
     *
     * @return tensor dim in the reference cell
     */
    inline Tensor<1, dim>
    transform_ref_face_correction_to_ref_cell(
      const Vector<double> &correction_ref_face,
      const unsigned int    local_face_id)
    {
      Tensor<1, dim> correction_ref_cell;

      unsigned int j = 0;
      for (unsigned int i = 0; i < dim; ++i)
        {
          if (local_face_id / 2 == i)
            {
              correction_ref_cell[i] = 0.0;
              continue;
            }
          correction_ref_cell[i] = correction_ref_face[j];
          j += 1;
        }

      return correction_ref_cell;
    };

    /**
     * @brief
     * Compute the numerical jacobian at the point x_n of the distance
     * minimization problem for the DoF x_I in the reference face space (dim -
     * 1).
     *
     * @param[in] stencil_real stencil in the real space. The evaluation point
     * x_n corresponds to the first entry.
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
      LAPACKFullMatrix<double> &jacobian_matrix)
    {
      for (unsigned int i = 0; i < dim - 1; ++i)
        {
          const Tensor<1, dim> x_n_to_x_I_real_m1 =
            x_I_real - stencil_real[2 * i + 1];

          Tensor<1, dim - 1> residual_ref_m1;
          compute_residual(x_n_to_x_I_real_m1,
                           distance_gradients[2 * i + 1],
                           transformation_jacobians[2 * i + 1],
                           residual_ref_m1);

          const Tensor<1, dim> x_n_to_x_I_real_p1 =
            x_I_real - stencil_real[2 * i + 2];

          Tensor<1, dim - 1> residual_ref_p1;
          compute_residual(x_n_to_x_I_real_p1,
                           distance_gradients[2 * i + 2],
                           transformation_jacobians[2 * i + 2],
                           residual_ref_p1);

          for (unsigned int j = 0; j < dim - 1; ++j)
            {
              jacobian_matrix.set(j,
                                  i,
                                  (residual_ref_p1[j] - residual_ref_m1[j]) /
                                    (2.0 * perturbation));
            }
        }
    };

    /**
     * @brief
     * Compute the analytical jacobian at the point x_n of the distance
     * minimization problem for the DoF x_I in the reference face space (dim -
     * 1): J_R = J^T*H*J, where J is the face transformation jacobian, and
     * H is the Hessian matrix in the real space (H = H(d) + H(||x_I - x_n||)).
     *
     * @param[in] x_real evaluation point x_n in the real spac.
     *
     * @param[in] x_I_real coordinate of the DoF x_I in the real space
     *
     * @param[in] transformation_jacobian face transformation jacobian
     *
     * @param[in] face_local_dof_values values of the DoFs of the face
     *
     * @param[out] jacobian_matrix jacobian matrix of the minimization problem
     */
    inline void
    compute_analytical_jacobian(
      const Point<dim>                      &x_real,
      const Point<dim>                      &x_I_real,
      const DerivativeForm<1, dim - 1, dim> &transformation_jacobian,
      const std::vector<double>             &face_local_dof_values,
      LAPACKFullMatrix<double>              &jacobian_matrix)
    {
      const Tensor<1, dim> x_n_to_x_I_real_p1 = x_I_real - x_real;

      const double x_n_to_x_I_real_p1_norm = x_n_to_x_I_real_p1.norm();

      const double x_n_to_x_I_real_p1_norm_inv = 1.0 / x_n_to_x_I_real_p1_norm;

      const double x_n_to_x_I_real_p1_norm_cubic_inv =
        x_n_to_x_I_real_p1_norm_inv * x_n_to_x_I_real_p1_norm_inv *
        x_n_to_x_I_real_p1_norm_inv;

      LAPACKFullMatrix<double> hessian_matrix(dim, dim);

      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int j = 0; j < dim; ++j)
            {
              double h_ij = -(x_n_to_x_I_real_p1[i] * x_n_to_x_I_real_p1[j]) *
                            x_n_to_x_I_real_p1_norm_cubic_inv;
              if (i == j)
                h_ij += x_n_to_x_I_real_p1_norm_inv;

              hessian_matrix.set(i, j, h_ij);
            }
        }

      const DerivativeForm<1, dim, dim - 1> transformation_jacobian_tanspose =
        transformation_jacobian.transpose();

      LAPACKFullMatrix<double> H_x_transformation_jacobian(dim, dim - 1);
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int j = 0; j < dim - 1; ++j)
            {
              double matrix_ij = 0;
              for (unsigned int k = 0; k < dim; ++k)
                {
                  matrix_ij +=
                    hessian_matrix(i, k) * transformation_jacobian[k][j];
                }
              H_x_transformation_jacobian.set(i, j, matrix_ij);
            }
        }

      for (unsigned int i = 0; i < dim - 1; ++i)
        {
          for (unsigned int j = 0; j < dim - 1; ++j)
            {
              double matrix_ij = 0;
              for (unsigned int k = 0; k < dim; ++k)
                {
                  matrix_ij += transformation_jacobian_tanspose[i][k] *
                               H_x_transformation_jacobian(k, j);
                }
              jacobian_matrix.set(i, j, matrix_ij);
            }
        }

      // In Q1, off-diagonal terms are coming from the term H(d) of the
      // real space Hessian. FEPointEvaluation doesn't include the
      // implementation of the Hessian, hence it is computed by hand.
      double off_diag_H = 0;
      if constexpr (dim == 3)
        {
          off_diag_H = jacobian_matrix(0, 1) + face_local_dof_values[0] -
                       face_local_dof_values[1] - face_local_dof_values[2] +
                       face_local_dof_values[3];
          jacobian_matrix.set(0, 1, off_diag_H);
          jacobian_matrix.set(1, 0, off_diag_H);
        }
    };

    /**
     * @brief
     * Compute the distance according to: d(x_I) = d(x_n) + ||x_I - x_n||
     *
     * @param[in] x_n_to_x_I_real vector from x_n to x_I in the real space
     *
     * @param[in] distance value of the distance at the point x_n
     *
     * @return distance between x_I and the interface
     */
    inline double
    compute_distance(const Tensor<1, dim> &x_n_to_x_I_real,
                     const double          distance)
    {
      return distance + x_n_to_x_I_real.norm();
    };

    /// Finite element discretizing the problem
    std::shared_ptr<FiniteElement<dim>> fe;

    /// Mapping between the real and reference space
    std::shared_ptr<MappingQ<dim>> mapping;

    /// Maximum redistanciation distance
    const double max_distance;

    /// Iso-level describing the interface from which the signed distance is
    /// computed (after scaling!!)
    const double iso_level;

    /// Scaling factor to apply to the input level-set field
    const double scaling;

    /// Verbosity level
    const Parameters::Verbosity verbosity;

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

    /// Solution vector of the signed distance as used for output only
    GlobalVectorType signed_distance_output;

    /// Solution vector of the signed distance with ghost values (read-only
    /// version of signed_distance)
    LinearAlgebra::distributed::Vector<double> signed_distance_with_ghost;

    /// Solution vector of the distance (write only)
    LinearAlgebra::distributed::Vector<double> distance;

    /// Solution vector of the distance with ghost values (read-only version of
    /// distance)
    LinearAlgebra::distributed::Vector<double> distance_with_ghost;

    /// Value of the correction to apply to the signed_distance to match the
    /// cell-wise volume encompassed by the level 0 of level_set
    LinearAlgebra::distributed::Vector<double> volume_correction;

    /// Hanging node constraints
    AffineConstraints<double> constraints;

    /// Surface vertices of the interface reconstruction stored in a cell-wise
    /// map (volume cell)
    std::map<types::global_cell_index, std::vector<Point<dim>>>
      interface_reconstruction_vertices;

    /// Surface cells of the interface reconstruction stored in a cell-wise map
    /// (volume cell)
    std::map<types::global_cell_index, std::vector<CellData<dim - 1>>>
      interface_reconstruction_cells;

    /// Set of DoFs belonging to intersected cells
    std::set<types::global_dof_index> intersected_dofs;

    /// PVDHandler for interface reconstruction
    PVDHandler pvd_handler_reconstruction;

    /// PVDHandler for signed distance
    PVDHandler pvd_handler_signed_distance;

    std::map<unsigned int, std::vector<unsigned int>> face_opposite_dofs_map;

    std::map<unsigned int, std::vector<unsigned int>> face_dofs_map;
  };
} // namespace InterfaceTools

#endif
