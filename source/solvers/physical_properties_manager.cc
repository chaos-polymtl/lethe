/*---------------------------------------------------------------------
 *
 * Copyright (C) 2021 - by the Lethe authors
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
 * ---------------------------------------------------------------------
 */

#include <solvers/physical_properties_manager.h>


PhysicalPropertiesManager::PhysicalPropertiesManager(
  Parameters::PhysicalProperties physical_properties)
  : number_of_fluids(physical_properties.number_of_fluids)
{
  // For each fluid, declare the physical properties
  for (unsigned int f = 0; f < number_of_fluids; ++f)
    {
    }
  // For
}
