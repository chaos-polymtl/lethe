// SPDX-FileCopyrightText: Copyright (c) 2021, 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_copy_data_h
#define lethe_copy_data_h

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <vector>

using namespace dealii;

/**
 * @brief Class responsible for storing the information calculated using the assembly of regular (meaning
 * non-stabilized) equations. It is also used to
 * initialize, zero (reset) and store the cell_matrix, the cell_rhs and the
 * dof indices associated with the dofs of the cell.
 **/

class CopyData
{
public:
  /**
   * @brief Constructor. Allocates the memory for the cell_matrix, cell_rhs and the local_dof_indices
   *
   * @param n_dofs Number of degrees of freedom per cell in the problem
   */
  CopyData(const unsigned int n_dofs)
    : local_matrix(n_dofs, n_dofs)
    , local_rhs(n_dofs)
    , local_dof_indices(n_dofs){};

  /**
   * @brief Resets the cell_matrix and the cell_rhs to zero
   */
  void
  reset()
  {
    local_matrix = 0;
    local_rhs    = 0;
  }

  FullMatrix<double>                   local_matrix;
  Vector<double>                       local_rhs;
  std::vector<types::global_dof_index> local_dof_indices;

  // Boolean used to indicate if the cell being assembled is local or not
  // This information is used to indicate to the copy_local_to_global function
  // if it should indeed copy or not.
  bool cell_is_local;
};


/**
 * @brief Class responsible for storing the information calculated using the assembly of stabilized
 * scalar equations. Like the CopyData class, this class is used to initialize,
 * zero (reset) and store the cell_matrix and the cell_rhs.
 * Contrary to the regular CopyData class, this class
 * also stores the strong_residual and the strong_jacobian of the equation being
 * assembled. This is useful for equations that implement residual-based
 * stabilization such as SUPG. This class is specialized for single component
 * equations because the strong jacobian is stored using a Vector<double>
 **/
class StabilizedMethodsCopyData
{
public:
  /**
   * @brief Constructor. Allocates the memory for the cell_matrix, cell_rhs
   * and dof-indices using the number of dofs and the strong_residual using the
   * number of quadrature points and, the strong_jacobian using both
   *
   * @param n_dofs Number of degrees of freedom per cell in the problem
   *
   * @param n_q_points Number of quadrature points
   */
  StabilizedMethodsCopyData(const unsigned int n_dofs,
                            const unsigned int n_q_points)
    : local_matrix(n_dofs, n_dofs)
    , local_rhs(n_dofs)
    , local_dof_indices(n_dofs)
    , strong_residual(n_q_points)
    , strong_jacobian(n_q_points, Vector<double>(n_dofs)){};


  /**
   * @brief Resets the cell_matrix, cell_rhs, strong_residual
   * and strong_jacobian to zero
   */
  void
  reset()
  {
    local_matrix = 0;
    local_rhs    = 0;

    strong_residual = 0;
    for (unsigned int i = 0; i < strong_jacobian.size(); ++i)
      {
        strong_jacobian[i] = 0;
      }
  }

  FullMatrix<double>                   local_matrix;
  Vector<double>                       local_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  Vector<double>                       strong_residual;
  std::vector<Vector<double>>          strong_jacobian;

  // Boolean used to indicate if the cell being assembled is local or not
  // This information is used to indicate to the copy_local_to_global function
  // if it should indeed copy or not.
  bool cell_is_local;
};



/**
 * @brief Responsible for storing the information calculated using the discontinuous
 * galerkin (DG) assembly of stabilized scalar equations. Like the
 * StabilizedCopyData class, this class is used to initialize, zero (reset) and
 * store the cell_matrix and the cell_rhs while having support for stabilized
 * terms. However, it also has to support the storage of information at internal
 * faces. This is required for the DG methods.
 **/
class StabilizedDGMethodsCopyData : public StabilizedMethodsCopyData
{
public:
  /**
   * @brief Constructor. Allocates the memory using the base class Constructor.
   * @param[in] n_dofs Number of degrees of freedom per cell in the problem
   *
   * @param[in] n_q_points Number of quadrature points
   */
  StabilizedDGMethodsCopyData(const unsigned int n_dofs,
                              const unsigned int n_q_points)
    : StabilizedMethodsCopyData(n_dofs, n_q_points){};


  // Data structure for internal face contributions
  struct CopyDataFace
  {
    FullMatrix<double> face_matrix;
    Vector<double>     face_rhs;

    std::vector<types::global_dof_index> joint_dof_indices;
  };

  std::vector<CopyDataFace> face_data;
};

/**
 * @brief Class responsible for storing the information calculated using the assembly of stabilized
 * Tensor<1,dim> equations. Like the CopyData class, this class is used to
 * initialize, zero (reset) and store the cell_matrix and the cell_rhs. Contrary
 * to the regular CopyData class, this class also stores the strong_residual and
 * the strong_jacobian of the equation being assembled. This is useful for
 * equations that implement residual-based stabilization such as SUPG. This
 * class is specialized for Tensor<1,dim> equations because the strong
 * jacobian is stored using a Tensor<1,dim>
 **/


template <int dim>
class StabilizedMethodsTensorCopyData
{
public:
  /**
   * @brief Constructor. Allocates the memory for the cell_matrix and cell_rhs
   * using the number of dofs and the strong_residual using the number of
   * quadrature points and, the strong_jacobian using both
   *
   * @param n_dofs Number of degrees of freedom per cell in the problem
   *
   * @param n_q_points Number of quadrature points
   */
  StabilizedMethodsTensorCopyData(const unsigned int n_dofs,
                                  const unsigned int n_q_points)
    : local_matrix(n_dofs, n_dofs)
    , local_rhs(n_dofs)
    , local_dof_indices(n_dofs)
    , strong_residual(n_q_points)
    , strong_jacobian(n_q_points, std::vector<Tensor<1, dim>>(n_dofs)){};

  /**
   * @brief Resets the cell_matrix, cell_rhs, strong_residual
   * and strong_jacobian to zero
   */
  void
  reset()
  {
    local_matrix = 0;
    local_rhs    = 0;

    for (unsigned int q = 0; q < strong_jacobian.size(); ++q)
      {
        strong_residual[q] = 0;

        for (unsigned int i = 0; i < strong_jacobian[q].size(); ++i)
          strong_jacobian[q][i] = 0;
      }
  }

  FullMatrix<double>                       local_matrix;
  Vector<double>                           local_rhs;
  std::vector<types::global_dof_index>     local_dof_indices;
  std::vector<Tensor<1, dim>>              strong_residual;
  std::vector<std::vector<Tensor<1, dim>>> strong_jacobian;


  // Boolean used to indicate if the cell being assembled is local or not
  // This information is used to indicate to the copy_local_to_global function
  // if it should indeed copy or not.
  bool cell_is_local;
  bool cell_is_cut;
};


/**
 * @brief Class responsible for storing the local system of an ultraweak Discontinuous
 * Petrov-Galerkin (DPG) discretization. The assemblers fill the uncondensed
 * local system: the Gram matrix of the test space \f$G\f$, the matrices of the
 * bilinear form in the cell interior \f$B\f$ and on the skeleton \f$\hat{B}\f$,
 * and the load vector \f$l\f$. The solver then condenses this system on the
 * skeleton unknowns using the operators
 * \f$M_1 = B^\dagger G^{-1} B\f$, \f$M_2 = B^\dagger G^{-1} \hat{B}\f$,
 * \f$M_3 = \hat{B}^\dagger G^{-1} \hat{B}\f$, \f$M_4 l = B^\dagger G^{-1} l\f$
 * and \f$M_5 l = \hat{B}^\dagger G^{-1} l\f$. Since \f$G\f$ and \f$M_1\f$ are
 * symmetric positive definite, they are factorized with a Cholesky
 * decomposition and their inverses are only applied through solves. The copy
 * data then stores either the condensed skeleton system (assembly) or the
 * reconstructed interior solution and the DPG residual (interior
 * reconstruction).
 **/
class DPGCopyData
{
public:
  /**
   * @brief Constructor. Allocates the memory of all the local matrices and
   * vectors from the number of degrees of freedom per cell of each finite
   * element space.
   *
   * @param[in] n_dofs_test Number of degrees of freedom per cell of the test
   * space.
   *
   * @param[in] n_dofs_trial_interior Number of degrees of freedom per cell of
   * the interior trial space.
   *
   * @param[in] n_dofs_trial_skeleton Number of degrees of freedom per cell of
   * the skeleton trial space.
   */
  DPGCopyData(const unsigned int n_dofs_test,
              const unsigned int n_dofs_trial_interior,
              const unsigned int n_dofs_trial_skeleton)
    : G_matrix(n_dofs_test, n_dofs_test)
    , B_matrix(n_dofs_test, n_dofs_trial_interior)
    , B_hat_matrix(n_dofs_test, n_dofs_trial_skeleton)
    , l_vector(n_dofs_test)
    , G_inverse_B(n_dofs_test, n_dofs_trial_interior)
    , G_inverse_B_hat(n_dofs_test, n_dofs_trial_skeleton)
    , G_inverse_l(n_dofs_test)
    , M1_matrix(n_dofs_trial_interior, n_dofs_trial_interior)
    , M2_matrix(n_dofs_trial_interior, n_dofs_trial_skeleton)
    , M3_matrix(n_dofs_trial_skeleton, n_dofs_trial_skeleton)
    , M4_l(n_dofs_trial_interior)
    , M5_l(n_dofs_trial_skeleton)
    , M1_inverse_M2(n_dofs_trial_interior, n_dofs_trial_skeleton)
    , M1_inverse_M4_l(n_dofs_trial_interior)
    , tmp_matrix_M2M1M2(n_dofs_trial_skeleton, n_dofs_trial_skeleton)
    , local_matrix(n_dofs_trial_skeleton, n_dofs_trial_skeleton)
    , local_rhs(n_dofs_trial_skeleton)
    , local_dof_indices(n_dofs_trial_skeleton)
    , local_skeleton_solution(n_dofs_trial_skeleton)
    , local_interior_solution(n_dofs_trial_interior)
    , tmp_vector_interior(n_dofs_trial_interior)
    , tmp_vector_error_indicator(n_dofs_test)
    , local_residual(n_dofs_test)
    , local_dof_indices_trial_interior(n_dofs_trial_interior)
    , local_residual_norm_squared(0.)
    , active_cell_index(0)
    , cell_is_local(false){};

  /**
   * @brief Resets the uncondensed local system to zero. The \f$M_1\f$ matrix
   * is also reset since LAPACKFullMatrix keeps track of its factorization
   * status and forbids to factorize it again if it has already been
   * factorized.
   */
  void
  reset()
  {
    G_matrix     = 0;
    B_matrix     = 0;
    B_hat_matrix = 0;
    l_vector     = 0;
    M1_matrix    = 0;
  }

  // Uncondensed local DPG system filled by the assemblers
  LAPACKFullMatrix<double> G_matrix;
  LAPACKFullMatrix<double> B_matrix;
  LAPACKFullMatrix<double> B_hat_matrix;
  Vector<double>           l_vector;

  // Solves with the Gram matrix: G^{-1} B, G^{-1} B_hat and G^{-1} l
  LAPACKFullMatrix<double> G_inverse_B;
  LAPACKFullMatrix<double> G_inverse_B_hat;
  Vector<double>           G_inverse_l;

  // Operators of the static condensation, solves with M_1 and temporary
  // matrix tmp_matrix_M2M1M2 = M_2^\dagger M_1^{-1} M_2
  LAPACKFullMatrix<double> M1_matrix;
  LAPACKFullMatrix<double> M2_matrix;
  LAPACKFullMatrix<double> M3_matrix;
  Vector<double>           M4_l;
  Vector<double>           M5_l;
  LAPACKFullMatrix<double> M1_inverse_M2;
  Vector<double>           M1_inverse_M4_l;
  LAPACKFullMatrix<double> tmp_matrix_M2M1M2;

  // Condensed skeleton system distributed in the global system
  FullMatrix<double>                   local_matrix;
  Vector<double>                       local_rhs;
  std::vector<types::global_dof_index> local_dof_indices;

  // Interior reconstruction and DPG residual. The temporary vectors are
  // tmp_vector_interior = M_1^{-1} M_2 * x_skeleton, used when reconstructing
  // the interior solution, and tmp_vector_error_indicator = B * x_interior +
  // B_hat * x_skeleton, used when computing the error indicator.
  Vector<double>                       local_skeleton_solution;
  Vector<double>                       local_interior_solution;
  Vector<double>                       tmp_vector_interior;
  Vector<double>                       tmp_vector_error_indicator;
  Vector<double>                       local_residual;
  std::vector<types::global_dof_index> local_dof_indices_trial_interior;
  double                               local_residual_norm_squared;
  unsigned int                         active_cell_index;

  // Boolean used to indicate if the cell being assembled is local or not
  // This information is used to indicate to the copy_local_to_global function
  // if it should indeed copy or not. It replace the call to the function
  // cell->is_locally_owned().
  bool cell_is_local;
};


#endif
