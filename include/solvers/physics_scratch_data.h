// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_physics_scratch_data_h
#define lethe_physics_scratch_data_h

/**
 * @brief Base class for scratch data objects of physics and subequations.
 */
class PhysicsScratchDataBase
{
public:
  /**
   * @brief Default virtual destructor
   */
  virtual ~PhysicsScratchDataBase() = default;

  /**
   * @brief Allocate the memory necessary memory for all members of the
   * scratch.
   */
  virtual void
  allocate() = 0;
};

#endif
