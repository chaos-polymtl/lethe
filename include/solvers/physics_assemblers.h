// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_physics_assemblers_h
#define lethe_physics_assemblers_h

/**
 * @brief Base class for assembler objects of physics' subequations.
 *
 * @note It allows to call subequations' assembler through the
 * VOFSubequationsInterface.
 */
template <typename ScratchDataType, typename CopyDataType>
class PhysicsAssemblerBase
{
public:
  /**
   * @brief Assemble the matrix.
   *
   * @param[in] scratch_data Scratch data containing the VOF phase gradient
   * projection information.
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
   * @param[in] scratch_data Scratch data containing the VOF phase gradient
   * projection information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for rhs assembly.
   *
   * @param[in,out] copy_data Destination where the local_rhs is copied to.
   */
  virtual void
  assemble_rhs(const ScratchDataType &scratch_data,
               CopyDataType          &copy_data) = 0;
};

#endif
