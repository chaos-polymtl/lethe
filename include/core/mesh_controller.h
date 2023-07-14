/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------*/


#ifndef lethe_mesh_controler_h
#define lethe_mesh_controler_h

#include <iostream>


/**
 * @Class Controller that target a given number of elements in the mesh. This
 * controller is used to define the coarsening factor.
 */

class MeshController
{
public:
  /**
   * @brief Constructor of the controller.
   * @param target_number_of_elements, the number of element the controller try to obtained.
   */
  MeshController(const unsigned int target_number_of_elements)
    : target_number_of_elements(target_number_of_elements)
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
