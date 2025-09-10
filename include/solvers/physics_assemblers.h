// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_physics_assemblers_h
#define lethe_physics_assemblers_h

/**
 * @brief Base class for assembler objects of physics and subequations.
 *
 * @tparam ScratchDataType Type of scratch data object used for linear system
 * assembly.
 *
 * @tparam CopyDataType Type of copy data object used for transferring from
 * local system assembly to global system.
 */
template <typename ScratchDataType, typename CopyDataType>
class PhysicsAssemblerBase
{
public:
  /**
   * @brief Assemble the matrix.
   *
   * @param[in] scratch_data Scratch data containing the information required
   * for system assembly.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   *
   * @param[in,out] copy_data Destination where the local_matrix is copied to.
   */
  virtual void
  assemble_matrix(const ScratchDataType &scratch_data,
                  CopyDataType          &copy_data) = 0;

  /**
   * @brief Assemble the right-hand side (rhs).
   *
   * @param[in] scratch_data Scratch data containing the information required
   * for system assembly.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for rhs assembly.
   *
   * @param[in,out] copy_data Destination where the local_rhs is copied to.
   */
  virtual void
  assemble_rhs(const ScratchDataType &scratch_data,
               CopyDataType          &copy_data) = 0;

  /**
   * @brief Default destructor
   */

  virtual ~PhysicsAssemblerBase() = default;
};


/**
 * @brief Base class for assembler objects of physics and subequations that require a face assembly.
 * This type of assembler is generally used for DG equations.
 *
 * @tparam ScratchDataType Type of scratch data object used for face
 * assembly.
 *
 * @tparam CopyDataType Type of copy data object used for transferring from
 * local system assembly to global system.
 */
template <typename ScratchDataType, typename CopyDataType>
class PhysicsFaceAssemblerBase
{
public:
  /**
   * @brief Interface for the call to matrix assembly
   * @param[in]  scratch_data Scratch data containing the physical solver
   * information. It is important to note that the scratch data has to have been
   * re-inited before calling for matrix assembly.
   * @param[in,out]  copy_data Destination where the local_rhs and local_matrix
   * should be copied
   */

  virtual void
  assemble_matrix(const ScratchDataType &scratch_data,
                  CopyDataType          &copy_data) = 0;


  /**
   * @brief Interface for the call to rhs
   * @param[in]  scratch_data Scratch data containing the physical solver
   * information. It is important to note that the scratch data has to have been
   * re-inited before calling for matrix assembly.
   * @param[in,out] copy_data Destination where the local_rhs and local_matrix
   * should be copied
   */

  virtual void
  assemble_rhs(const ScratchDataType &scratch_data,
               CopyDataType          &copy_data) = 0;

  /**
   * @brief Default destructor
   */

  virtual ~PhysicsFaceAssemblerBase() = default;
};

#endif
