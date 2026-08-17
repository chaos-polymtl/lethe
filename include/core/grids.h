// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_grids_h
#define lethe_grids_h

#include <core/boundary_conditions.h>
#include <core/manifolds.h>
#include <core/parameters.h>

#include <deal.II/distributed/tria_base.h>

using namespace dealii;

/**
 * @brief Attaches a grid to a triangulation using mesh parameters
 *
 * @param[in,out] triangulation The triangulation to which a grid is attached
 *
 * @param[in] mesh_parameters The mesh parameters used to decide what type of
 * mesh or primitive is used
 */
template <int dim, int spacedim = dim>
void
attach_grid_to_triangulation(Triangulation<dim, spacedim> &triangulation,
                             const Parameters::Mesh       &mesh_parameters);

/**
 * @brief Modifies the triangulation to set up its periodic boundary conditions
 *
 * The periodic boundary information is the part of the boundary conditions that
 * is shared by the CFD and the DEM solvers, so this function is used by both.
 *
 * @param[in,out] triangulation The triangulation to which a grid is attached
 *
 * @param[in] periodic_boundary_information The periodic boundary pairs and
 * their direction of periodicity. This is used to set up the periodicity of the
 * domain
 */
template <int dim, int spacedim = dim>
void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<dim, spacedim> &triangulation,
  const Parameters::PeriodicBoundaryInformation &periodic_boundary_information);

/**
 * @brief Completely sets up a mesh and its manifolds.
 *
 * This function reads and constructs a mesh according to the given mesh
 * parameters, sets up manifolds on boundaries, and applies boundary conditions.
 * It is intended for distributed triangulations.
 *
 * @tparam dim        The topological dimension of the mesh.
 * @tparam spacedim   The embedding space dimension (default = dim).
 * @param[in,out] triangulation         The triangulation to which the grid will
 * be attached.
 * @param[in]     mesh_parameters      Parameters defining the mesh or primitive
 * to use.
 * @param[in]     manifolds_parameters Information about the manifolds attached
 * to boundaries.
 * @param[in]     restart              Flag indicating whether this is a
 * restart.
 * @param[in]     periodic_boundary_information The periodic boundary pairs and
 * their direction of periodicity, used to set up the periodicity of the domain.
 */
template <int dim, int spacedim = dim>
void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<dim, spacedim> &triangulation,
  const Parameters::Mesh                                &mesh_parameters,
  const Parameters::Manifolds                           &manifolds_parameters,
  bool                                                   restart,
  const Parameters::PeriodicBoundaryInformation &periodic_boundary_information);

/**
 * @brief Completely sets up a mesh and its manifolds for rotor-stator domains.
 *
 * The `triangulation` initially contains the stator
 * triangulation, which is later replaced by the merged rotor-stator
 * triangulation.
 *
 * @tparam dim        The topological dimension of the mesh.
 * @tparam spacedim   The embedding space dimension (default = dim).
 * @param[in,out] triangulation        The triangulation to which the grid will
 * be attached.
 * @param[in]     mesh_parameters      Parameters defining the mesh or primitive
 * to use.
 * @param[in]     manifolds_parameters Information about the manifolds attached
 * to boundaries.
 * @param[in]     restart              Flag indicating whether this is a
 * restart.
 * @param[in]     periodic_boundary_information The periodic boundary pairs and
 * their direction of periodicity, used to set up the periodicity of the domain.
 * @param[in]     mortar_parameters    Parameters controlling the mortar method,
 * including rotor mesh info.
 */
template <int dim, int spacedim = dim>
void
read_mesh_and_manifolds_for_stator_and_rotor(
  parallel::DistributedTriangulationBase<dim, spacedim> &triangulation,
  const Parameters::Mesh                                &mesh_parameters,
  const Parameters::Manifolds                           &manifolds_parameters,
  bool                                                   restart,
  const Parameters::PeriodicBoundaryInformation &periodic_boundary_information,
  const Parameters::Mortar<dim>                 &mortar_parameters);

/**
 * @brief Refine a mesh around specific boundary ids
 *
 * @param[in] n_refinement The number of times the boundary should be refined
 *
 * @param[in] boundary_ids The boundary ids for which the cells should be
 * refined
 *
 * @param[in, out] triangulation The triangulation on which refinement must be
 * carried out
 *
 */
template <int dim, int spacedim = dim>
void
refine_triangulation_at_boundaries(
  const std::vector<int>                                 boundary_ids,
  const unsigned int                                     n_refinement,
  parallel::DistributedTriangulationBase<dim, spacedim> &triangulation)
{
  // For the amount of refinements required
  for (unsigned int r = 0; r < n_refinement; ++r)
    {
      // Loop over the cells and flag all cells which are within the list of
      // boundary ids
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              for (const auto f : cell->face_indices())

                {
                  if (cell->face(f)->at_boundary())
                    {
                      if (std::find(boundary_ids.begin(),
                                    boundary_ids.end(),
                                    cell->face(f)->boundary_id()) !=
                          boundary_ids.end())
                        {
                          cell->set_refine_flag();
                        }
                    }
                }
            }
        }
      triangulation.execute_coarsening_and_refinement();
    }
}

/**
 * @brief Build the triangulation which delimits a mesh refinement box
 *
 * The box is described by a full set of mesh parameters, so any mesh type
 * supported by attach_grid_to_triangulation() can be used to delimit a
 * refinement region. The initial refinement of the box mesh is applied here, so
 * the resulting triangulation is ready to be used to flag cells.
 *
 * @param[in] box_mesh_parameters The mesh parameters describing the box
 *
 * @param[in,out] box_triangulation The triangulation which will hold the box
 * mesh
 */
template <int dim, int spacedim = dim>
void
build_refinement_box_triangulation(
  const Parameters::Mesh       &box_mesh_parameters,
  Triangulation<dim, spacedim> &box_triangulation);

/**
 * @brief Apply scaling, rotation and translation to the mesh in the listed
 * order.
 *
 * @tparam dim An integer that denotes the dimensionality of the geometry.
 * @tparam spacedim An integer that denotes the dimension of the space occupied
 * by the embedded geometry.
 *
 * @param[in] mesh_parameters Parameters of the mesh.
 * @param[in, out] triangulation Triangulation object on which the mesh
 * transformations are applied.
 */
template <int dim, int spacedim = dim>
void
apply_mesh_transformation(const Parameters::Mesh       &mesh_parameters,
                          Triangulation<dim, spacedim> &triangulation);

#endif
