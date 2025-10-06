// SPDX-FileCopyrightText: Copyright (c) 2021, 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_copy_data_h
#define lethe_copy_data_h

#include <deal.II/lac/full_matrix.h>
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


#endif
