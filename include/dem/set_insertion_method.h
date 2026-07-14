// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_set_insertion_method_h
#define lethe_set_insertion_method_h

#include <core/parameters_lagrangian.h>

#include <dem/insertion_file.h>
#include <dem/insertion_list.h>
#include <dem/insertion_packed.h>
#include <dem/insertion_plane.h>
#include <dem/insertion_volume.h>

/**
 * @brief Set up the particle insertion method for the simulation.
 * Creates and configures the appropriate insertion object based on the
 * simulation parameters (e.g., volume insertion, plane insertion, etc.).
 *
 * @param size_distribution_object_container Contains all distribution for each
 * particle type.
 * @param triangulation Triangulation to access the cells in which the
 * particles will be inserted.
 * @param dem_parameters DEM parameters declared in the .prm file.
 * @param maximum_particle_diameter Maximum particle diameter based on values
 * defined in the parameter handler.
 * @param packing_method Set to true if the packed insertion method is selected.
 *
 * @return Shared pointer to the configured insertion object
 */
template <int dim, typename PropertiesIndex>
std::shared_ptr<Insertion<dim, PropertiesIndex>>
set_insertion_type(std::vector<std::shared_ptr<Distribution>>
                     &size_distribution_object_container,
                   parallel::distributed::Triangulation<dim> &triangulation,
                   DEMSolverParameters<dim>                  &dem_parameters,
                   const double &maximum_particle_diameter,
                   bool         &packing_method)
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
      case InsertionInfo<dim>::InsertionMethod::packed:
        {
          if constexpr (std::is_same_v<PropertiesIndex,
                                       DEM::DEMProperties::PropertiesIndex>)
            {
              packing_method = true;
              return std::make_shared<InsertionPacked<dim, PropertiesIndex>>(
                size_distribution_object_container,
                triangulation,
                dem_parameters);
            }
          else
            {
              packing_method = false;
              AssertThrow(false,
                          ExcMessage("InsertionPacked is not valid for "
                                     "the current solver type. It only works "
                                     "for the standard DEM solver."));
            }
        }
      default:
        AssertThrow(false, ExcMessage("Invalid insertion method."));
    }
}

/**
 * @brief Overload of set_insertion_type for callers that do not need to know
 * whether the packed insertion method was selected.
 *
 * @param size_distribution_object_container Contains all distribution for each
 * particle type.
 * @param triangulation Triangulation to access the cells in which the
 * particles will be inserted.
 * @param dem_parameters DEM parameters declared in the .prm file.
 * @param maximum_particle_diameter Maximum particle diameter based on values
 * defined in the parameter handler.
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
  bool dummy = false;
  return set_insertion_type<dim, PropertiesIndex>(
    size_distribution_object_container,
    triangulation,
    dem_parameters,
    maximum_particle_diameter,
    dummy);
}

#endif
