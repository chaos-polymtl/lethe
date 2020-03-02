/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

 *
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include <deal.II/base/data_out_base.h>

#include <dem/dem_solver_parameters.h>
#include <dem/visualization.h>
#include <sys/stat.h>

#include <fstream>
#include <iostream>

#ifndef DEM_WRITEVTU_H_
#  define DEM_WRITEVTU_H_

/**
 * Writing the simulation output as .pvtu and .vtu formats
 *
 * @note This function is taken from Aspect and dealii and implemented here
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class WriteVTU
{
public:
  WriteVTU<dim>();

  /**
   * Carries out writing .pvtu file which contains the general information of
   * simulation outputs and creating the output folder
   *
   * @param visulization_object An object of visualization class
   * @param dem_parameters DEM parameters declared in the .prm file
   */

  void
  write_master_files(Visualization<dim> &            visulization_object,
                     const DEMSolverParameters<dim> &dem_parameters);

  /**
   * Carries out writing .vtu files which contain the information of all the
   * particles in the system at the writing time-step
   *
   * @param visulization_object An object of visualization class
   * @param DEM_step Current DEM step number
   * @param DEM_time Current DEM time
   * @param dem_parameters DEM parameters declared in the .prm file
   */

  void
  writeVTUFiles(Visualization<dim> &            visulization_object,
                const int &                     DEM_step,
                const double &                  DEM_time,
                const DEMSolverParameters<dim> &dem_parameters);
};

#endif /* DEM_WRITEVTU_H_ */
