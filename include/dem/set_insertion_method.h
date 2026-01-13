// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_set_insertion_method_h
#define lethe_set_insertion_method_h

#include <core/parameters_lagrangian.h>

#include <dem/insertion_file.h>
#include <dem/insertion_list.h>
#include <dem/insertion_plane.h>
#include <dem/insertion_volume.h>

/**
 * @brief Set up the particle insertion method for the simulation.
 *
 * Creates and configures the appropriate insertion object based on the
 * simulation parameters (e.g., volume insertion, plane insertion, etc.).
 *
 * @param size_distribution_object_container Contains all distribution for each
 * particle type.
 * @param triangulation Triangulation to access the cells in which the
 * particles will be inserted
 * @param dem_parameters DEM parameters declared in the .prm file
 * @param maximum_particle_diameter Maximum particle diameter based on values
 * defined in the parameter handler
 *
 * @return Shared pointer to the configured insertion object
 */
template <int dim, typename PropertiesIndex>
std::shared_ptr<Insertion<dim, PropertiesIndex>>
set_insertion_type(std::vector<std::shared_ptr<Distribution>>
                     &size_distribution_object_container,
                   parallel::distributed::Triangulation<dim> &triangulation,
                   DEMSolverParameters<dim>                  &dem_parameters,
                   const double &maximum_particle_diameter)
{
  using namespace Parameters::Lagrangian;
  typename InsertionInfo<dim>::InsertionMethod insertion_method =
    dem_parameters.insertion_info.insertion_method;

  switch (insertion_method)
    {
      case InsertionInfo<dim>::InsertionMethod::file:
        return std::make_shared<InsertionFile<dim, PropertiesIndex>>(
          size_distribution_object_container, triangulation, dem_parameters);
      case InsertionInfo<dim>::InsertionMethod::list:
        return std::make_shared<InsertionList<dim, PropertiesIndex>>(
          size_distribution_object_container, triangulation, dem_parameters);
      case InsertionInfo<dim>::InsertionMethod::plane:
        return std::make_shared<InsertionPlane<dim, PropertiesIndex>>(
          size_distribution_object_container, triangulation, dem_parameters);
      case InsertionInfo<dim>::InsertionMethod::volume:
        return std::make_shared<InsertionVolume<dim, PropertiesIndex>>(
          size_distribution_object_container,
          triangulation,
          dem_parameters,
          maximum_particle_diameter);
      default:
        throw(std::runtime_error("Invalid insertion method."));
    }
}
#endif
