// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_mesh_controller_h
#define lethe_mesh_controller_h

/**
 * @Class Controller that target a given number of elements in the mesh. This
 * controller is used to define the coarsening factor.
 */

class MeshController
{
public:
  /**
   * @brief Constructor of the controller.
   * @param target_number_of_elements Number of element the controller tries to obtain.
   */
  MeshController(const unsigned int target_number_of_elements)
    : target_number_of_elements(target_number_of_elements)
    , previous_number_of_elements(0)
    , previous_mesh_control_error(0.0)
  {}
  virtual ~MeshController()
  {}
  /**
   * @brief Returns the coarsening factor.
   */
  double
  calculate_coarsening_factor(unsigned int current_number_of_elements);

private:
  const unsigned int target_number_of_elements;
  unsigned int       previous_number_of_elements;
  double             previous_mesh_control_error;
};

#endif
