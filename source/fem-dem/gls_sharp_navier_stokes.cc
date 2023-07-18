/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
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

#include <core/bdf.h>
#include <core/grids.h>
#include <core/lethe_grid_tools.h>
#include <core/sdirk.h>
#include <core/solutions_output.h>
#include <core/tensors_and_points_dimension_manipulation.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/navier_stokes_vof_assemblers.h>
#include <solvers/postprocessing_cfd.h>

#include <fem-dem/gls_sharp_navier_stokes.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/full_matrix.h>

// Constructor for class GLSNavierStokesSolver
template <int dim>
GLSSharpNavierStokesSolver<dim>::GLSSharpNavierStokesSolver(
  CFDDEMSimulationParameters<dim> &p_nsparam)
  : GLSNavierStokesSolver<dim>(p_nsparam.cfd_parameters)
  , cfd_dem_parameters(p_nsparam)
  , all_spheres(true)

{}

template <int dim>
GLSSharpNavierStokesSolver<dim>::~GLSSharpNavierStokesSolver()
{}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::vertices_cell_mapping()
{
  // Find all the cells around each vertices
  TimerOutput::Scope t(this->computing_timer, "vertices_to_cell_map");

  LetheGridTools::vertices_cell_mapping(this->dof_handler, vertices_to_cell);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::check_whether_all_particles_are_sphere()
{
  all_spheres = false;

  // WIP The optimized cut-cell mapping seems to lead to instability
  /*
  for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
    {
      if (particles[p_i].shape->get_shape_name().second !=
          Shape<dim>::ShapeType::sphere)
        {
          all_spheres = false;
          std::cout
            << "A non-spherical particle was found: using regular
  cut_cells_mapping."
            << std::endl;
          break;
        }
    }
    */
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::generate_cut_cells_map()
{
  // check all the cells if they are cut or not. Put the information in a map
  // with the key being the cell.
  TimerOutput::Scope t(this->computing_timer, "cut_cells_mapping");
  std::map<types::global_dof_index, Point<dim>> support_points;

  // A vector of the unordered map. Each map stores if a point is inside or
  // outside of a given particle.
  std::vector<std::unordered_map<types::global_dof_index, bool>>
    inside_outside_support_point_vector(particles.size());
  DoFTools::map_dofs_to_support_points(*this->mapping,
                                       this->dof_handler,
                                       support_points);

  // When the finite element order > 1, overconstrained cells are impossible
  // since there is always at least a DOF inside the element that is not
  // overconstrained. We therefore only have to check when the velocity order
  // == 1.
  const bool mapping_overconstrained_cells =
    this->simulation_parameters.fem_parameters.velocity_order == 1;

  cut_cells_map.clear();
  cells_inside_map.clear();
  overconstrained_fluid_cell_map.clear();
  if (mapping_overconstrained_cells)
    {
      local_dof_overconstrained.reinit(this->locally_owned_dofs,
                                       this->mpi_communicator);
      dof_overconstrained.reinit(this->locally_owned_dofs,
                                 this->locally_relevant_dofs,
                                 this->mpi_communicator);
    }

  const auto &       cell_iterator = this->dof_handler.active_cell_iterators();
  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int dofs_per_face = this->fe->dofs_per_face;

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> local_face_dof_indices(dofs_per_face);

  auto &             v_x_fe                  = this->fe->get_sub_fe(0, 1);
  const unsigned int dofs_per_cell_local_v_x = v_x_fe.dofs_per_cell;
  // // Loop on all the cells and check if they are cut.
  for (const auto &cell : cell_iterator)
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          bool         cell_is_cut                                = false;
          bool         cell_is_inside                             = false;
          unsigned int particle_id_which_cuts_this_cell           = 0;
          unsigned int particle_id_in_which_this_cell_is_embedded = 0;
          unsigned int number_of_particles_cutting_this_cell      = 0;
          // is the particle index that cuts the cell if it's cut. If the cell
          // is not cut the default value is stored (0). If the cell is not cut
          // this value will never be used.
          // number_of_particles_cutting_this_cell count the number of particles
          // that cut a cell if multiple particles cut this cell. This is used
          // to treat cells that are cut by multiple particles differently.
          cell->get_dof_indices(local_dof_indices);

          for (unsigned int p = 0; p < particles.size(); ++p)
            {
              // Particles defined from a stl or an iges don't have a signed
              // distance function. To identify cells that are cut a special
              // algorithm is required. For this type of particle, no cell are
              // identified as being inside the particle. This means that the
              // fluid is always solved inside the particle.
              if (particles[p].shape->additional_info_on_shape == "iges")
                {
                  cell_is_cut =
                    cell_cut_by_p_absolute_distance(cell, support_points, p);
                  particle_id_which_cuts_this_cell = p;
                  if (cell_is_cut)
                    {
                      number_of_particles_cutting_this_cell += 1;
                    }
                }
              else
                {
                  unsigned int nb_dof_inside = 0;
                  for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
                    {
                      // Count the number of DOFs that are inside
                      // of the particles. If all the DOfs are on one side
                      // the cell is not cut by the boundary.
                      if (0 == this->fe->system_to_component_index(j).first)
                        {
                          // Check if we already have checked this point in a
                          // previously evaluated cell. If we didn't find it in
                          // the previous evaluation, we assess whether the DOF
                          // is inside or outside the shape.
                          auto iterator =
                            inside_outside_support_point_vector[p].find(
                              local_dof_indices[j]);
                          if (iterator ==
                              inside_outside_support_point_vector[p].end())
                            {
                              if (particles[p].get_levelset(
                                    support_points[local_dof_indices[j]],
                                    cell) <= 0)
                                {
                                  ++nb_dof_inside;
                                  inside_outside_support_point_vector
                                    [p][local_dof_indices[j]] = true;
                                }
                              else
                                inside_outside_support_point_vector
                                  [p][local_dof_indices[j]] = false;
                            }
                          else
                            {
                              if (inside_outside_support_point_vector
                                    [p][local_dof_indices[j]] == true)
                                {
                                  ++nb_dof_inside;
                                }
                            }
                        }
                    }

                  // If some DOFs are inside the boundary, the cell is inside
                  // the particle or cut by the particle.
                  if (nb_dof_inside != 0)
                    {
                      // If all the DOFs are inside the boundary this cell is
                      // inside the particle. Otherwise the particle is cut.
                      if (nb_dof_inside == dofs_per_cell_local_v_x)
                        {
                          //  We register only the particle with the lowest id
                          //  as the particle in which this cell is embedded if
                          //  the cell is in multiple particles.
                          if (number_of_particles_cutting_this_cell == 0)
                            {
                              cell_is_cut                      = false;
                              particle_id_which_cuts_this_cell = 0;
                              cell_is_inside                   = true;
                              particle_id_in_which_this_cell_is_embedded = p;
                            }
                          break;
                        }
                      else
                        {
                          // We only register the particle with the lowest id as
                          // the particle by this cell is cut if the cell is cut
                          // in multiple particles.
                          if (number_of_particles_cutting_this_cell == 0)
                            {
                              cell_is_cut                      = true;
                              particle_id_which_cuts_this_cell = p;
                              cell_is_inside                   = false;
                              particle_id_in_which_this_cell_is_embedded = 0;
                            }
                          number_of_particles_cutting_this_cell += 1;
                        }
                    }
                  else
                    {
                      if (number_of_particles_cutting_this_cell == 0)
                        {
                          cell_is_cut                                = false;
                          particle_id_which_cuts_this_cell           = 0;
                          cell_is_inside                             = false;
                          particle_id_in_which_this_cell_is_embedded = 0;
                        }
                    }
                }
            }

          if (mapping_overconstrained_cells)
            {
              // If a cell is cut, we register its DOFs as "cut" also. This
              // will allow us to detect fluid cells that have all their DOFs
              // cut, and consider them as cut cells.
              if (cell_is_cut)
                {
                  size_t id;
                  for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
                    {
                      id                            = local_dof_indices[j];
                      local_dof_overconstrained(id) = 1;
                    }
                }

              // We loop on every face of the cell, and if the face is at a
              // boundary we count it as constrained.
              size_t nb_faces = cell->n_faces();
              for (unsigned int f = 0; f < nb_faces; ++f)
                {
                  const auto face = cell->face(f);
                  if (face->at_boundary())
                    {
                      face->get_dof_indices(local_face_dof_indices);
                      size_t id;
                      for (unsigned int v = 0;
                           v < local_face_dof_indices.size();
                           ++v)
                        {
                          id = local_face_dof_indices[v];
                          local_dof_overconstrained(id) = 1;
                        }
                    }
                }
            }

          cut_cells_map[cell]    = {cell_is_cut,
                                 particle_id_which_cuts_this_cell,
                                 number_of_particles_cutting_this_cell};
          cells_inside_map[cell] = {cell_is_inside,
                                    particle_id_in_which_this_cell_is_embedded};
        }
    }

  if (mapping_overconstrained_cells)
    {
      dof_overconstrained = local_dof_overconstrained;
      dof_overconstrained.compress(VectorOperation::insert);

      for (const auto &cell : cell_iterator)
        {
          if (cell->is_locally_owned() || cell->is_ghost())
            {
              bool cell_is_cut;
              bool cell_is_inside;
              cell->get_dof_indices(local_dof_indices);
              std::tie(cell_is_cut, std::ignore, std::ignore) =
                cut_cells_map[cell];
              std::tie(cell_is_inside, std::ignore) = cells_inside_map[cell];
              if (!cell_is_cut && !cell_is_inside)
                {
                  unsigned int number_of_vertices_in_cell =
                    local_dof_indices.size();
                  unsigned int number_of_vertices_cut = 0;
                  unsigned int particle_current_id =
                    std::numeric_limits<int>::max();

                  double particle_current_distance(
                    std::numeric_limits<double>::max());
                  double particle_candidate_distance(
                    std::numeric_limits<double>::max());

                  unsigned int global_dof_id;
                  // We count the number of vertices that are constrained.
                  for (unsigned int j = 0; j < number_of_vertices_in_cell; ++j)
                    {
                      global_dof_id = local_dof_indices[j];
                      if (dof_overconstrained(global_dof_id) > 0)
                        {
                          number_of_vertices_cut += 1;
                        }
                    }
                  // We check which particle is the closest to the cell, so
                  // that we know which one to use for applying the immersed
                  // boundary constraints.
                  if (number_of_vertices_in_cell == number_of_vertices_cut)
                    {
                      for (unsigned int p = 0; p < particles.size(); p++)
                        {
                          particle_candidate_distance =
                            particles[p].get_levelset(cell->barycenter(), cell);
                          if (abs(particle_candidate_distance) <
                              abs(particle_current_distance))
                            {
                              particle_current_id = p;
                              particle_current_distance =
                                particle_candidate_distance;
                            }
                        }
                      overconstrained_fluid_cell_map[cell] = {
                        true, particle_current_id, particle_current_distance};
                    }
                }
            }
        }
    }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::refinement_control(
  const bool initial_refinement)
{
  //  This function applies the various refinement steps depending on the
  //  parameters and the state.
  if (initial_refinement)
    {
      // Apply the initial box refinement
      this->box_refine_mesh();
      update_precalculations_for_ib();
    }
  if (this->simulation_parameters.particlesParameters
        ->time_extrapolation_of_refinement_zone ||
      initial_refinement)
    {
      // Stores variable for refinement around the particle.
      double temp_refine =
        this->simulation_parameters.mesh_adaptation.variables.begin()
          ->second.refinement_fraction;
      double temp_coarse =
        this->simulation_parameters.mesh_adaptation.variables.begin()
          ->second.coarsening_fraction;
      this->simulation_parameters.mesh_adaptation.variables.begin()
        ->second.refinement_fraction = 0;
      this->simulation_parameters.mesh_adaptation.variables.begin()
        ->second.coarsening_fraction = 0;

      for (unsigned int i = 0;
           i <
           this->simulation_parameters.particlesParameters->initial_refinement;
           ++i)
        {
          this->pcout << "Initial refinement around IB particles - Step : "
                      << i + 1 << " of "
                      << this->simulation_parameters.particlesParameters
                           ->initial_refinement
                      << std::endl;
          refine_ib();
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          update_precalculations_for_ib();
        }
      this->simulation_parameters.mesh_adaptation.variables.begin()
        ->second.refinement_fraction = temp_refine;
      this->simulation_parameters.mesh_adaptation.variables.begin()
        ->second.coarsening_fraction = temp_coarse;
    }
  if (initial_refinement == false)
    {
      refine_ib();
      NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
        refine_mesh();
      update_precalculations_for_ib();
    }
}



template <int dim>
bool
GLSSharpNavierStokesSolver<dim>::cell_cut_by_p_absolute_distance(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  std::map<types::global_dof_index, Point<dim>> &       support_points,
  unsigned int                                          p)
{
  // This function aims at defining if a cell is cut when the level set used to
  // define the particle is not signed. This function works by analyzing the
  // location of the projection of the support point and the cell centroid on
  // the particle's surface. If the projected point at the surface of the
  // particle is inside the cell then the cell is cut.

  // First step: define useful variables.
  const unsigned int                   dofs_per_cell = this->fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  cell->get_dof_indices(local_dof_indices);
  bool       cell_is_cut      = false;
  Point<dim> centroid_of_cell = cell->barycenter();
  Point<dim> projected_point;
  particles[p].closest_surface_point(centroid_of_cell, projected_point, cell);

  // Check the centroid of the cell first.
  if (cell->point_inside(projected_point))
    {
      cell_is_cut = true;
      return cell_is_cut;
    }
  else
    {
      // If the projection of the centroid is not inside the cell we also check
      // for each dof support point of the cell.
      for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
        {
          // Only check for the support point of the velocity in x.
          if (0 == this->fe->system_to_component_index(j).first)
            {
              particles[p].closest_surface_point(
                support_points[local_dof_indices[j]], projected_point, cell);

              if (cell->point_inside(projected_point))
                {
                  cell_is_cut = true;
                  break;
                }
              if ((projected_point - support_points[local_dof_indices[j]])
                    .norm() < cell->diameter() * 0.1)
                {
                  cell_is_cut = true;
                  break;
                }
            }
        }
    }
  return cell_is_cut;
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::optimized_generate_cut_cells_map()
{
  TimerOutput::Scope t(this->computing_timer, "optimized_cut_cells_mapping");
  MappingQ1<dim>     mapping;
  unsigned int       max_children = GeometryInfo<dim>::max_children_per_cell;

  std::set<typename DoFHandler<dim>::cell_iterator> all_candidate_cells;
  std::set<typename DoFHandler<dim>::cell_iterator>
    previous_all_candidate_cells;

  cut_cells_map.clear();
  cells_inside_map.clear();
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          cut_cells_map[cell]    = {false, 0, 0};
          cells_inside_map[cell] = {false, 0};
        }
    }

  // Loop over particles in reverse.
  // This is done in reverse to guarantee that the lowest particle ID is
  // associated to the cut cell or the inside particle cell. For instance, if
  // you have two particles cutting the same cell, it guarantees that the lowest
  // ID particle is associated with the cell in the map.
  for (int p = particles.size() - 1; p >= 0; --p)
    {
      bool         empty    = true;
      unsigned int lvl_iter = 0;

      // Fix max level search
      unsigned int max_lvl_search =
        this->dof_handler.get_triangulation().n_levels() - 2;
      if (max_lvl_search < 4)
        max_lvl_search = this->dof_handler.get_triangulation().n_levels();

      // Search for candidates until at least one is found or until the
      // max_lvl_search is reached
      while (empty || lvl_iter < max_lvl_search)
        {
          const auto &cell_iterator =
            this->dof_handler.cell_iterators_on_level(lvl_iter);

          // Loop over the cells on level lvl_iter of the mesh
          for (const auto &cell : cell_iterator)
            {
              auto is_candidate = generate_cut_cell_candidates(cell, p);

              // check whether the cell is candidate for inside of cut
              if (is_candidate.first)
                {
                  all_candidate_cells.insert(cell);
                }
              if (is_candidate.second)
                {
                  all_candidate_cells.insert(cell);
                }
            }
          if (all_candidate_cells.size() > 0)
            empty = false;

          lvl_iter += 1;
        }

      // Store list with all candidate cells
      previous_all_candidate_cells = all_candidate_cells;

      // If there are candidate cells
      if (all_candidate_cells.size() != 0)
        {
          bool all_cells_are_active = false;

          // Loop until candidate cells are active
          while (!all_cells_are_active)
            {
              all_cells_are_active = true;
              all_candidate_cells.clear();

              // loop over the last set of candidates to check if it is still a
              // candidate
              for (auto cell : previous_all_candidate_cells)
                {
                  auto is_candidate = generate_cut_cell_candidates(cell, p);

                  // Check if cell is inside particles
                  if (is_candidate.first)
                    {
                      // If it is an active cell, add it to cells_inside_map
                      if (cell->is_active())
                        {
                          cells_inside_map[cell] = {true, p};
                        }
                      else
                        {
                          // If we are here, the cell has children.
                          all_cells_are_active = false;
                          for (unsigned int j = 0; j < max_children; ++j)
                            {
                              // Store the children of the cell for
                              // further checks.
                              all_candidate_cells.insert(cell->child(j));
                            }
                        }
                    }

                  if (is_candidate.second)
                    {
                      // If it is an active cell, add it to cells_cut_map
                      if (cell->is_active())
                        {
                          cut_cells_map[cell] = {true, p, 0};
                        }
                      else
                        {
                          // If we are here, the cell has children.
                          all_cells_are_active = false;
                          for (unsigned int j = 0; j < max_children; ++j)
                            {
                              // Store the children of the cell for
                              // further checks.
                              all_candidate_cells.insert(cell->child(j));
                            }
                        }
                    }
                }
              previous_all_candidate_cells.clear();
              previous_all_candidate_cells = all_candidate_cells;
            }
        }
    }
}


template <int dim>
std::pair<bool, bool>
GLSSharpNavierStokesSolver<dim>::generate_cut_cell_candidates(
  const typename DoFHandler<dim>::cell_iterator &cell,
  const unsigned int                             p_id)
{
  bool cell_is_inside = false;
  bool cell_is_cut    = false;

  double search_radius = particles[p_id].radius * (1.0 - 1e-07);
  auto   position      = particles[p_id].position;

  bool point_inside_cell = cell->point_inside(position);

  // Check how many vertices are inside the particle
  unsigned int nb_vertices_inside = 0;
  for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
    {
      if ((cell->vertex(i) - position).norm() < search_radius)
        ++nb_vertices_inside;
    }

  // If vertices are found inside the particle
  if (nb_vertices_inside > 0)
    {
      // If the number of vertices inside the cell is equal to the number of
      // vertices per cell, the cell is inside the particle
      if (nb_vertices_inside == GeometryInfo<dim>::vertices_per_cell)
        {
          cell_is_inside = true;
          cell_is_cut    = false;
        }
      // Otherwise, the cell is cut
      else
        {
          cell_is_inside = false;
          cell_is_cut    = true;
        }
      return {cell_is_inside, cell_is_cut};
    }

  // If the particles is inside the cell and it is not known whether all
  // vertices are inside the particles or not, set all true by default
  if (point_inside_cell)
    {
      cell_is_inside = true;
      cell_is_cut    = true;
      return {cell_is_inside, cell_is_cut};
    }

  // The last check consists of projecting the particle's position on the cells'
  // face and checking whether the projected points are within the particle. If
  // one of the projected points is inside the particle, either the cell is
  // inside or is cut by the particles. If this is true, the next
  // level will be tested again by this function.

  // Initialize superpoint of manifold
  std::vector<Point<dim>> manifold_points(
    GeometryInfo<dim - 1>::vertices_per_cell);

  // Loop through faces
  for (const auto face : cell->face_indices())
    {
      auto  face_iter           = cell->face(face);
      auto &local_face_manifold = face_iter->get_manifold();

      // Loop through face vertices
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face; ++i)
        {
          // Assign vertex to manifold superpoint
          manifold_points[i] = face_iter->vertex(i);
        }

      // Create array of points on face
      auto surrounding_face_points =
        make_array_view(manifold_points.begin(), manifold_points.end());

      // Loop through face vertices
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face; ++i)
        {
          // Project points to face
          Point<dim> projected_point =
            local_face_manifold.project_to_manifold(surrounding_face_points,
                                                    position);

          bool           projected_point_over_face = true;
          Point<dim>     unit_cell_projected_point;
          Tensor<1, dim> projected_point_tensor(projected_point);

          // Check whether the projected point is within the cell
          try
            {
              projected_point_over_face = cell->point_inside(projected_point);

              if ((position - projected_point).norm() < search_radius &&
                  projected_point_over_face)
                {
                  cell_is_inside = true;
                  cell_is_cut    = true;
                  return {cell_is_inside, cell_is_cut};
                }
            }
          // If the method crashes, change default for false
          catch (...)
            {
              projected_point_over_face = false;
            }

          for (unsigned int d = 0; d != dim; d++)
            {
              if (projected_point_tensor[d] < 0 ||
                  projected_point_tensor[d] > 1)
                {
                  projected_point_over_face = false;
                }
            }
        }
    }
  return {cell_is_inside, cell_is_cut};
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::define_particles()
{
  some_particles_are_coupled = false;
  // initialized the particles
  if (this->simulation_parameters.particlesParameters
        ->load_particles_from_file == false)
    {
      particles.resize(this->simulation_parameters.particlesParameters->nb);
      for (unsigned int i = 0;
           i < this->simulation_parameters.particlesParameters->nb;
           ++i)
        {
          particles[i] =
            this->simulation_parameters.particlesParameters->particles[i];
          if (particles[i].integrate_motion == true)
            {
              if (typeid(*particles[i].shape) != typeid(Sphere<dim>))
                throw std::runtime_error(
                  "Shapes other than sphere cannot have their motion integrated through fluid-structure interaction");
              some_particles_are_coupled = true;
            }
        }
    }
  else
    {
      load_particles_from_file();
      for (unsigned int i = 0;
           i < this->simulation_parameters.particlesParameters->nb;
           ++i)
        {
          if (particles[i].integrate_motion == true)
            {
              if (typeid(*particles[i].shape) != typeid(Sphere<dim>))
                throw std::runtime_error(
                  "Shapes other than sphere cannot have their motion integrated through fluid-structure interaction");
              some_particles_are_coupled = true;
            }
        }
    }

  table_p.resize(particles.size());
  ib_dem.initialize(
    this->simulation_parameters.particlesParameters,
    std::make_shared<Parameters::Lagrangian::FloatingWalls<dim>>(
      cfd_dem_parameters.dem_parameters.floating_walls),
    this->mpi_communicator,
    particles);

  check_whether_all_particles_are_sphere();
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::refine_ib()
{
  TimerOutput::Scope t(this->computing_timer, "refine_around_ib");
  Point<dim>         center_immersed;
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(*this->mapping,
                                       this->dof_handler,
                                       support_points);
  double dt = this->simulation_control->get_time_steps_vector()[0];

  const unsigned int                   dofs_per_cell = this->fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  bool extrapolate_particle_position =
    this->simulation_parameters.particlesParameters
      ->time_extrapolation_of_refinement_zone;

  const auto &cell_iterator = this->dof_handler.active_cell_iterators();
  for (const auto &cell : cell_iterator)
    {
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(local_dof_indices);
          for (unsigned int p = 0; p < particles.size(); ++p)
            {
              unsigned int count_small = 0;
              center_immersed          = particles[p].position;
              Tensor<1, dim> r;
              r[0] = particles[p].radius;

              Point<dim>   particle_position;
              Tensor<1, 3> particle_orientation;
              if (extrapolate_particle_position)
                {
                  particle_position    = particles[p].position;
                  particle_orientation = particles[p].orientation;
                  particles[p].orientation =
                    particles[p].orientation + particles[p].omega * dt;
                  particles[p].position[0] =
                    particles[p].position[0] + particles[p].velocity[0] * dt;
                  particles[p].position[1] =
                    particles[p].position[1] + particles[p].velocity[1] * dt;

                  if constexpr (dim == 3)
                    {
                      particles[p].position[2] = particles[p].position[2] +
                                                 particles[p].velocity[2] * dt;
                    }
                  particles[p].set_position(particles[p].position);
                  particles[p].set_orientation(particles[p].orientation);
                }
              // Check if a point on the random point on the IB is contained in
              // that cell. If the particle is much smaller than the cell, all
              // its vertices may be outside of the particle. In that case the
              // cell won't be refined. To prevent that, we check if a random
              // point on the boundary is contained in the cell.
              bool cell_as_ib_inside =
                cell->point_inside(particles[p].position + r);
              for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
                {
                  // Only check the dof of velocity in x.
                  if (this->fe->system_to_component_index(j).first == 0)
                    {
                      // Count the number of DOFs that fall in the refinement
                      // zone around the particle. To fall in the zone, the
                      // radius of the DOFs to the center of the particle must
                      // be bigger than the inside radius of the refinement zone
                      // and smaller than the outside radius of the refinement
                      // zone.
                      if (particles[p].is_inside_crown(
                            support_points[local_dof_indices[j]],
                            this->simulation_parameters.particlesParameters
                              ->outside_radius,
                            this->simulation_parameters.particlesParameters
                              ->inside_radius,
                            cell))
                        {
                          ++count_small;
                        }
                    }
                }

              if (extrapolate_particle_position)
                {
                  particles[p].position    = particle_position;
                  particles[p].orientation = particle_orientation;
                  particles[p].set_position(particles[p].position);
                  particles[p].set_orientation(particles[p].orientation);
                }
              if (count_small > 0 || cell_as_ib_inside)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }
    }
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::force_on_ib()
{
  // This function defines the force and torque applied on an Immersed Boundary
  // based on the sharp edge method on a hyper_sphere of dim=2 or dim=3

  TimerOutput::Scope t(this->computing_timer, "new force_eval");

  const FEValuesExtractors::Scalar pressure(dim);
  const FEValuesExtractors::Vector velocities(0);

  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  // Define a map to all dofs and their support points
  std::map<types::global_dof_index, Point<dim>> support_points;

  DoFTools::map_dofs_to_support_points(*this->mapping,
                                       this->dof_handler,
                                       support_points);

  // Initalize fe value objects in order to do calculation with it later
  QGauss<dim>            q_formula(this->number_quadrature_points);
  MappingQ<dim - 1, dim> local_face_map(
    (2 > this->simulation_parameters.fem_parameters.velocity_order) ?
      2 :
      this->simulation_parameters.fem_parameters.velocity_order);

  FESystem<dim - 1, dim> local_face_fe(
    FE_Q<dim - 1, dim>(
      this->simulation_parameters.fem_parameters.velocity_order),
    dim,
    FE_Q<dim - 1, dim>(
      this->simulation_parameters.fem_parameters.pressure_order),
    1);
  FEValues<dim - 1, dim> fe_face_projection_values(local_face_map,
                                                   local_face_fe,
                                                   *this->face_quadrature,
                                                   update_values |
                                                     update_quadrature_points |
                                                     update_gradients |
                                                     update_JxW_values |
                                                     update_normal_vectors);



  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int dofs_per_face = this->fe->dofs_per_face;

  Assert(this->simulation_parameters.physical_properties_manager
           .density_is_constant(),
         RequiresConstantDensity(
           "GLSSharpNavierStokesSolver<dim>::force_on_ib"));

  int  order = this->simulation_parameters.particlesParameters->order;
  auto density_model =
    this->simulation_parameters.physical_properties_manager.get_density(0);
  std::map<field, double> field_values;
  field_values[field::temperature] = 0;
  double fluid_density             = density_model->value(field_values);
  double length_ratio =
    this->simulation_parameters.particlesParameters->length_ratio;
  IBStencil<dim>      stencil;
  std::vector<double> ib_coef = stencil.coefficients(order, length_ratio);

  // Rheological model for viscosity properties
  double     viscosity;
  const auto rheological_model =
    this->simulation_parameters.physical_properties_manager.get_rheology();

  const unsigned int vertices_per_face = GeometryInfo<dim>::vertices_per_face;
  const unsigned int n_q_points_face   = this->face_quadrature->size();

  std::vector<double>                      pressure_values(ib_coef.size());
  Tensor<1, dim>                           normal_vector;
  std::vector<Tensor<2, dim>>              velocity_gradients(ib_coef.size());
  std::vector<std::vector<Tensor<1, dim>>> velocity_gradients_component(dim +
                                                                        1);
  for (unsigned int i = 0; i < dim + 1; ++i)
    velocity_gradients_component[i].resize(ib_coef.size());
  Tensor<2, dim>              fluid_viscous_stress;
  Tensor<2, dim>              fluid_pressure;
  Tensor<2, dim>              fluid_viscous_stress_at_ib;
  Tensor<2, dim>              fluid_pressure_stress_at_ib;
  Tensor<2, dim>              shear_rate;
  DoFHandler<dim - 1, dim>    local_face_dof_handler;
  Triangulation<dim - 1, dim> local_face_projection_triangulation;

  std::vector<Point<dim>>        vertices_of_face_projection(vertices_per_face);
  std::vector<CellData<dim - 1>> local_face_cell_data(1);

  typename Mapping<dim>::InternalDataBase mapping_data;
  // Define multiple local_dof_indices one for the cell iterator one for the
  // cell with the second point for the sharp edge stencil and one for
  // manipulation on the neighbour’s cell.

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_2(dofs_per_cell);

  std::vector<types::global_dof_index> local_face_dof_indices(dofs_per_face);

  std::vector<Point<dim>> unite_cell_interpolation_points(ib_coef.size());
  std::vector<Point<dim>> cell_interpolation_points(ib_coef.size());
  std::vector<double>     local_interp_sol(ib_coef.size());

  std::unordered_map<unsigned int, std::pair<Tensor<2, dim>, Tensor<2, dim>>>
    force_eval;

  // Define cell iterator
  const auto &cell_iterator = this->dof_handler.active_cell_iterators();

  // Clear particle force and torque
  for (unsigned int i = 0; i < particles.size(); ++i)
    {
      particles[i].fluid_forces          = 0;
      particles[i].fluid_viscous_forces  = 0;
      particles[i].fluid_pressure_forces = 0;
      particles[i].fluid_torque          = 0;
    }

  double       total_area    = 0;
  unsigned int nb_evaluation = 0;
  // Loop over all the cell
  for (const auto &cell : cell_iterator)
    {
      if (cell->is_locally_owned())
        {
          // Particle id that cut the cell.
          unsigned int p;
          unsigned int p_count;
          bool         cell_is_cut;
          std::tie(cell_is_cut, p, p_count) = cut_cells_map[cell];
          // If the cell is cut
          if (cell_is_cut)
            {
              // Loop over all the face of the cell that is cut.
              for (const auto face : cell->face_indices())
                {
                  auto local_face = cell->face(face);
                  cell->face(face)->get_dof_indices(local_face_dof_indices);

                  // Check if the face is cut
                  unsigned int nb_dof_inside = 0;
                  for (unsigned int j = 0; j < local_face_dof_indices.size();
                       ++j)
                    {
                      // Count the number of DOFs that are inside
                      // of the particles. If all the DOfs are on one side
                      // the cell is not cut by the boundary.
                      if (particles[p].get_levelset(
                            support_points[local_face_dof_indices[j]], cell) <=
                          0)
                        ++nb_dof_inside;
                    }

                  // If the face is not cut and the face is outside of the IB,
                  // the face of the cell is on the boundary of the
                  // computational domain.
                  if (nb_dof_inside == 0)
                    {
                      // Create a vector to approximate the normal of the cell
                      // to orient the projection of the face.
                      Tensor<1, dim> approximate_surface_cell_normal;
                      // Projects the face on the surface of the IB. This
                      // creates a surface cell where we can evaluate the
                      // solution. Define the triangulation of the surface cell.
                      for (unsigned int i = 0; i < vertices_per_face; ++i)
                        {
                          Point<dim> vertex_projection;
                          particles[p].closest_surface_point(
                            local_face->vertex(i), vertex_projection, cell);
                          approximate_surface_cell_normal +=
                            (local_face->vertex(i) - vertex_projection) /
                            ((local_face->vertex(i) - vertex_projection)
                               .norm() +
                             DBL_MIN);

                          // Define the vertices of the surface cell.
                          // Create the list of vertices
                          for (unsigned int j = 0; j < dim; ++j)
                            {
                              vertices_of_face_projection[i][j] =
                                vertex_projection[j];
                            }
                          // Create the connectivity of the vertices of the cell
                          local_face_cell_data[0].vertices[i] = i;
                        }
                      approximate_surface_cell_normal =
                        approximate_surface_cell_normal / vertices_per_face;

                      local_face_cell_data[0].material_id = 0;

                      local_face_projection_triangulation =
                        Triangulation<dim - 1, dim>();
                      // Create a dof handler that contains the triangulation of
                      // the projection of the face on the IB. This create the
                      // surface cell on the IB
                      local_face_projection_triangulation.create_triangulation(
                        vertices_of_face_projection,
                        local_face_cell_data,
                        SubCellData());
                      local_face_projection_triangulation.set_all_manifold_ids(
                        0);
                      local_face_projection_triangulation.set_manifold(
                        0, *particles[p].shape->get_shape_manifold());

                      local_face_dof_handler.reinit(
                        local_face_projection_triangulation);
                      local_face_dof_handler.distribute_dofs(local_face_fe);

                      // Defined the solution on  IB surface cell by using the
                      // IB stencil to extrapolate the fluid stress tensor.

                      std::vector<Tensor<2, dim>>
                        local_face_viscous_stress_tensor(dofs_per_face);
                      std::vector<Tensor<2, dim>> local_face_pressure_tensor(
                        dofs_per_face);
                      for (unsigned int i = 0;
                           i < local_face_dof_indices.size();
                           ++i)
                        {
                          const unsigned int component_i =
                            this->fe->face_system_to_component_index(i).first;
                          // Check if that dof already have been used to
                          // extrapolate the fluid stress tensor on the IB
                          // surface.
                          if (force_eval.find(local_face_dof_indices[i]) ==
                              force_eval.end())
                            {
                              if (component_i == 0)
                                {
                                  // Only need one extrapolation by dof
                                  // location;

                                  // Count the number of evaluation

                                  auto [point, interpolation_points] =
                                    stencil.support_points_for_interpolation(
                                      order,
                                      length_ratio,
                                      particles[p],
                                      support_points[local_face_dof_indices[i]],
                                      cell);

                                  auto cell_2 =
                                    ib_done[local_face_dof_indices[i]].second;
                                  // Check if we already have the cell used to
                                  // defined the IB constraint of that dof. We
                                  // always have that information except if the
                                  // dof is not owned.
                                  if (ib_done[local_face_dof_indices[i]]
                                        .first == false)
                                    {
                                      // Get the cell use for the extrapolation
                                      auto point_to_find_cell =
                                        stencil.point_for_cell_detection(
                                          particles[p],
                                          support_points
                                            [local_face_dof_indices[i]],
                                          cell);

                                      try
                                        {
                                          cell_2 = LetheGridTools::
                                            find_cell_around_point_with_neighbors<
                                              dim>(this->dof_handler,
                                                   vertices_to_cell,
                                                   cell,
                                                   point_to_find_cell);
                                        }
                                      catch (...)
                                        {
                                          cell_2 = cell;
                                        }
                                    }

                                  cell_2->get_dof_indices(local_dof_indices_2);

                                  unite_cell_interpolation_points[0] =
                                    this->mapping->transform_real_to_unit_cell(
                                      cell_2, point);
                                  for (unsigned int j = 1; j < ib_coef.size();
                                       ++j)
                                    {
                                      unite_cell_interpolation_points[j] =
                                        this->mapping
                                          ->transform_real_to_unit_cell(
                                            cell_2,
                                            interpolation_points[j - 1]);
                                    }

                                  fluid_viscous_stress_at_ib  = 0;
                                  fluid_pressure_stress_at_ib = 0;

                                  // Create a quadrature that is based on the IB
                                  // stencil
                                  Quadrature<dim> q_local(
                                    unite_cell_interpolation_points, ib_coef);
                                  FEValues<dim> fe_values_cell2(
                                    *this->fe,
                                    q_local,
                                    update_quadrature_points |
                                      update_gradients | update_values);

                                  // Evaluate the relevant information at the
                                  // quadrature points to do the extrapolation.
                                  fe_values_cell2.reinit(cell_2);
                                  fe_values_cell2[velocities]
                                    .get_function_gradients(
                                      this->evaluation_point,
                                      velocity_gradients);
                                  fe_values_cell2[pressure].get_function_values(
                                    this->evaluation_point, pressure_values);

                                  // Extrapolate the fluid stress tensor on the
                                  // surface of the IB.

                                  for (unsigned int k = 0; k < ib_coef.size();
                                       ++k)
                                    {
                                      fluid_pressure = 0;
                                      for (int d = 0; d < dim; ++d)
                                        {
                                          fluid_pressure[d][d] =
                                            pressure_values[k];
                                        }
                                      shear_rate =
                                        velocity_gradients[k] +
                                        transpose(velocity_gradients[k]);


                                      const double shear_rate_magnitude =
                                        calculate_shear_rate_magnitude(
                                          shear_rate);

                                      std::map<field, double> field_values;
                                      field_values[field::shear_rate] =
                                        shear_rate_magnitude;

                                      viscosity =
                                        rheological_model->value(field_values);

                                      fluid_viscous_stress =
                                        -viscosity * shear_rate;

                                      fluid_viscous_stress_at_ib -=
                                        fluid_viscous_stress * ib_coef[k];

                                      fluid_pressure_stress_at_ib -=
                                        fluid_pressure * ib_coef[k];
                                    }
                                  // Store the stress tensor that results from
                                  // the extrapolation in the local evaluation
                                  // vector of the IB surface cell and in a map
                                  // that is used if the same extrapolation is
                                  // needed in another face.
                                  local_face_viscous_stress_tensor[i] =
                                    fluid_viscous_stress_at_ib;
                                  local_face_pressure_tensor[i] =
                                    fluid_pressure_stress_at_ib;

                                  force_eval[local_face_dof_indices[i]] =
                                    std::make_pair(fluid_viscous_stress_at_ib,
                                                   fluid_pressure_stress_at_ib);
                                }
                            }
                          else
                            {
                              // Use the results from a previously evaluated
                              // extrapolation. This step comes with an error
                              // due to the curvature of the surface in Q2 and
                              // higher order elements.
                              local_face_viscous_stress_tensor[i] =
                                force_eval[local_face_dof_indices[i]].first;
                              local_face_pressure_tensor[i] =
                                force_eval[local_face_dof_indices[i]].second;
                            }
                        }
                      // Use the extrapolation of fluid stress tensor at the
                      // dof location of the IB surface cell to integrate the
                      // stress tensor on the surface of the IB
                      auto local_face_viscous_stress_tensor_old =
                        local_face_viscous_stress_tensor;
                      auto local_face_pressure_tensor_old =
                        local_face_pressure_tensor;

                      local_face_viscous_stress_tensor.clear();
                      local_face_pressure_tensor.clear();

                      local_face_viscous_stress_tensor.resize(
                        local_face_dof_indices.size());
                      local_face_pressure_tensor.resize(
                        local_face_dof_indices.size());

                      for (const auto &projection_cell_face :
                           local_face_dof_handler.active_cell_iterators())
                        {
                          fe_face_projection_values.reinit(
                            projection_cell_face);
                          std::vector<Point<dim>> q_points =
                            fe_face_projection_values.get_quadrature_points();
                          try
                            {
                              if (this->simulation_parameters.fem_parameters
                                    .velocity_order > 1)
                                {
                                  FullMatrix<double> interpolation_matrix(
                                    local_face_dof_indices.size(),
                                    local_face_dof_indices.size());
                                  FullMatrix<double> inv_interpolation_matrix(
                                    local_face_dof_indices.size(),
                                    local_face_dof_indices.size());

                                  // Define the interpolation matrix of the
                                  // surface cell
                                  for (unsigned int i = 0;
                                       i < local_face_dof_indices.size();
                                       ++i)
                                    {
                                      Point<dim> point_projection;
                                      particles[p].closest_surface_point(
                                        support_points
                                          [local_face_dof_indices[i]],
                                        point_projection,
                                        cell);

                                      auto projected_point_unit =
                                        local_face_map
                                          .transform_real_to_unit_cell(
                                            projection_cell_face,
                                            point_projection);

                                      for (unsigned int j = 0;
                                           j < local_face_dof_indices.size();
                                           ++j)
                                        {
                                          interpolation_matrix[i][j] = 0;
                                          if (
                                            this->fe
                                              ->face_system_to_component_index(
                                                j)
                                              .first ==
                                            this->fe
                                              ->face_system_to_component_index(
                                                i)
                                              .first)
                                            interpolation_matrix[i][j] +=
                                              local_face_fe.shape_value(
                                                j, projected_point_unit);
                                        }
                                    }
                                  inv_interpolation_matrix.invert(
                                    interpolation_matrix);
                                  // Define the value of the fluid stress tensor
                                  // on the surface cell at the DOF support
                                  // points location.
                                  for (unsigned int i = 0;
                                       i < local_face_dof_indices.size();
                                       ++i)
                                    {
                                      for (unsigned int j = 0;
                                           j < local_face_dof_indices.size();
                                           ++j)
                                        {
                                          local_face_viscous_stress_tensor[i] +=
                                            inv_interpolation_matrix[i][j] *
                                            local_face_viscous_stress_tensor_old
                                              [j];
                                          local_face_pressure_tensor[i] +=
                                            inv_interpolation_matrix[i][j] *
                                            local_face_pressure_tensor_old[j];
                                        }
                                    }
                                }
                              else
                                {
                                  local_face_viscous_stress_tensor =
                                    local_face_viscous_stress_tensor_old;
                                  local_face_pressure_tensor =
                                    local_face_pressure_tensor_old;
                                }
                            }
                          catch (...)
                            {
                              local_face_viscous_stress_tensor =
                                local_face_viscous_stress_tensor_old;
                              local_face_pressure_tensor =
                                local_face_pressure_tensor_old;
                            }
                          for (unsigned int q = 0; q < n_q_points_face; q++)
                            {
                              // Evaluate the total surface
                              // Redefined the normal at the quadrature point
                              // since we don't control the orientation of the
                              // cell.
                              normal_vector =
                                fe_face_projection_values.normal_vector(q);
                              if (scalar_product(
                                    normal_vector,
                                    approximate_surface_cell_normal) < 0)
                                normal_vector = -normal_vector;

                              fluid_viscous_stress = 0;
                              fluid_pressure       = 0;

                              double local_weight = 0;
                              // Integrate
                              for (unsigned int i = 0;
                                   i < local_face_dof_indices.size();
                                   ++i)
                                {
                                  const unsigned int component_i =
                                    local_face_fe.system_to_component_index(i)
                                      .first;
                                  if (component_i == 0)
                                    {
                                      fluid_viscous_stress +=
                                        fe_face_projection_values.shape_value(
                                          i, q) *
                                        local_face_viscous_stress_tensor[i];

                                      fluid_pressure +=
                                        fe_face_projection_values.shape_value(
                                          i, q) *
                                        local_face_pressure_tensor[i];

                                      total_area +=
                                        fe_face_projection_values.JxW(q) *
                                        fe_face_projection_values.shape_value(
                                          i, q);
                                      local_weight +=
                                        fe_face_projection_values.shape_value(
                                          i, q);
                                    }
                                }

                              auto viscous_force =
                                fluid_viscous_stress * normal_vector *
                                fe_face_projection_values.JxW(q);

                              auto pressure_force =
                                fluid_pressure * normal_vector *
                                fe_face_projection_values.JxW(q);

                              auto force = viscous_force + pressure_force;

                              if (force != force)
                                {
                                  // The force is nan; this happens when the
                                  // face projection is a line. This generally
                                  // happens when the surface on which the force
                                  // is evaluated is perfectly flat and aligned
                                  // with the mesh. Since the area associated
                                  // with this face projection is zero, we set
                                  // the local force contribution to zero.
                                  force          = 0;
                                  viscous_force  = 0;
                                  pressure_force = 0;
                                }
                              if (force.norm() > 0)
                                {
                                  nb_evaluation += local_weight;
                                }

                              if (force.norm() > 0)
                                {
                                  nb_evaluation += local_weight;
                                }


                              // Add the local contribution of this surface
                              // cell.

                              particles[p].fluid_viscous_forces +=
                                tensor_nd_to_3d(viscous_force);

                              particles[p].fluid_pressure_forces +=
                                tensor_nd_to_3d(pressure_force);

                              auto distance =
                                q_points[q] - particles[p].position;
                              if (dim == 2)
                                {
                                  particles[p].fluid_torque[0] += 0.;
                                  particles[p].fluid_torque[1] += 0.;
                                  particles[p].fluid_torque[2] +=
                                    distance[0] * force[1] -
                                    distance[1] * force[0];
                                }
                              else if (dim == 3)
                                {
                                  particles[p].fluid_torque[0] +=
                                    distance[1] * force[2] -
                                    distance[2] * force[1];
                                  particles[p].fluid_torque[1] +=
                                    distance[2] * force[0] -
                                    distance[0] * force[2];
                                  particles[p].fluid_torque[2] +=
                                    distance[0] * force[1] -
                                    distance[1] * force[0];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  // Sums the force evaluation on each of the processor.
  for (unsigned int i = 0; i < particles.size(); ++i)
    {
      particles[i].fluid_viscous_forces =
        Utilities::MPI::sum(particles[i].fluid_viscous_forces,
                            this->mpi_communicator) *
        fluid_density;
      particles[i].fluid_pressure_forces =
        Utilities::MPI::sum(particles[i].fluid_pressure_forces,
                            this->mpi_communicator) *
        fluid_density;

      particles[i].fluid_forces =
        particles[i].fluid_viscous_forces + particles[i].fluid_pressure_forces;

      particles[i].fluid_torque =
        Utilities::MPI::sum(particles[i].fluid_torque, this->mpi_communicator) *
        fluid_density;
    }

  total_area = Utilities::MPI::sum(total_area, this->mpi_communicator);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::write_force_ib()
{
  TimerOutput::Scope t(this->computing_timer, "output_forces_ib");
  for (unsigned int p = 0; p < particles.size(); ++p)
    {
      {
        if (this->this_mpi_process == 0)
          {
            std::string filename =
              this->simulation_parameters.simulation_control.output_folder +
              this->simulation_parameters.particlesParameters
                ->ib_force_output_file +
              "." + Utilities::int_to_string(p, 2) + ".dat";
            std::ofstream output(filename.c_str());

            table_p[p].write_text(output);
            std::string filename_residual =
              this->simulation_parameters.simulation_control.output_folder +
              "residual" + ".dat";
            std::ofstream output_residual(filename_residual.c_str());
          }
        MPI_Barrier(this->mpi_communicator);
      }
    }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::output_field_hook(DataOut<dim> &data_out)
{
  std::vector<std::shared_ptr<Shape<dim>>> all_shapes;
  for (const IBParticle<dim> &particle : particles)
    {
      all_shapes.push_back(particle.shape);
    }
  std::shared_ptr<Shape<dim>> combined_shapes =
    std::make_shared<CompositeShape<dim>>(all_shapes, Point<dim>(), Point<3>());

  levelset_postprocessor =
    std::make_shared<LevelsetPostprocessor<dim>>(combined_shapes);
  data_out.add_data_vector(this->present_solution, *levelset_postprocessor);
  Vector<float> cell_cuts(this->triangulation->n_active_cells());
  Vector<float> cell_overconstrained(this->triangulation->n_active_cells());
  const unsigned int                   dofs_per_cell = this->fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  // If the enable extra verbose output is activated we add an output field to
  // the results where we identify which particle cuts each cell of the domain.
  if (this->simulation_parameters.particlesParameters
        ->enable_extra_sharp_interface_vtu_output_field)
    {
      // Define cell iterator
      const auto & cell_iterator = this->dof_handler.active_cell_iterators();
      unsigned int i             = 0;

      for (const auto &cell : cell_iterator)
        {
          if (cell->is_locally_owned())
            {
              // Here we identify all cut cells, and if a cell is cut we output
              // the ID of the cell that cuts it. If the cell is not cut, a
              // default value of -1 is output.
              bool cell_is_cut;
              int  particle_id;
              std::tie(cell_is_cut, particle_id, std::ignore) =
                cut_cells_map[cell];
              if (cell_is_cut)
                cell_cuts(i) = particle_id;
              else
                cell_cuts(i) = -1;

              // Here we identify all overconstrained cells. If a cell is
              // overconstrained, we output the ID of the closest particle (the
              // one used for extrapolation of the velocity immersed boundary
              // condition). If the cell is not overconstrained, a default value
              // of -1 is output.
              bool cell_is_overconstrained;
              std::tie(cell_is_overconstrained, particle_id, std::ignore) =
                overconstrained_fluid_cell_map[cell];
              if (cell_is_overconstrained)
                cell_overconstrained(i) = particle_id;
              else
                cell_overconstrained(i) = -1;
            }
          i += 1;
        }
      data_out.add_data_vector(cell_cuts, "cell_cut");
      data_out.add_data_vector(cell_overconstrained, "cell_overconstrained");

      levelset_gradient_postprocessor =
        std::make_shared<LevelsetGradientPostprocessor<dim>>(combined_shapes);
      data_out.add_data_vector(this->present_solution,
                               *levelset_gradient_postprocessor);
    }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::postprocess_fd(bool firstIter)
{
  if (this->simulation_control->is_output_iteration())
    {
      this->write_output_results(this->present_solution);
    }

  bool enable =
    this->simulation_parameters.analytical_solution->calculate_error();
  this->simulation_parameters.analytical_solution->set_enable(false);
  NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
    postprocess_fd(firstIter);
  this->simulation_parameters.analytical_solution->set_enable(enable);
  // Calculate the error with respect to the analytical solution
  if (!firstIter &&
      this->simulation_parameters.analytical_solution->calculate_error())
    {
      // Update the time of the exact solution to the actual time
      this->exact_solution->set_time(
        this->simulation_control->get_current_time());
      const std::pair<double, double> error =
        this->calculate_L2_error_particles();

      if (this->simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        {
          this->error_table.add_value(
            "cells", this->triangulation->n_global_active_cells());
          this->error_table.add_value("error_velocity", error.first);
          this->error_table.add_value("error_pressure", error.second);

          auto summary = this->computing_timer.get_summary_data(
            this->computing_timer.total_wall_time);
          double total_time = 0;
          for (auto it = summary.begin(); it != summary.end(); ++it)
            {
              total_time += summary[it->first];
            }
          this->error_table.add_value("total_time", total_time);
        }
      else
        {
          this->error_table.add_value(
            "time", this->simulation_control->get_current_time());
          this->error_table.add_value("error_velocity", error.first);
          this->error_table.add_value("error_pressure", error.second);

          if (this->simulation_parameters.timer.write_time_in_error_table)
            {
              auto summary = this->computing_timer.get_summary_data(
                this->computing_timer.total_wall_time);
              double total_time = 0;
              for (auto it = summary.begin(); it != summary.end(); ++it)
                {
                  total_time += summary[it->first];
                }
              this->error_table.add_value("total_time", total_time);
            }
        }
      if (this->simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error velocity : " << error.first
                      << " L2 error pressure: " << error.second << std::endl;
        }
    }
}


template <int dim>
std::pair<double, double>
GLSSharpNavierStokesSolver<dim>::calculate_L2_error_particles()
{
  TimerOutput::Scope t(this->computing_timer, "error");
  QGauss<dim>        quadrature_formula(this->number_quadrature_points + 1);
  FEValues<dim>      fe_values(*this->mapping,
                          *this->fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
  FEFaceValues<dim>  fe_face_values(*this->mapping,
                                   *this->fe,
                                   *this->face_quadrature,
                                   update_values | update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  const unsigned int dofs_per_cell =
    this->fe->dofs_per_cell; // This gives you dofs per cell
  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int n_q_points      = quadrature_formula.size();
  const unsigned int n_q_points_face = this->face_quadrature->size();

  std::vector<Vector<double>> q_exactSol(n_q_points, Vector<double>(dim + 1));

  std::vector<Tensor<1, dim>> local_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> local_face_velocity_values(n_q_points_face);
  std::vector<double>         local_pressure_values(n_q_points);
  std::vector<double>         div_phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);

  Function<dim> *l_exact_solution = this->exact_solution;

  Point<dim>                                    center_immersed;
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(*this->mapping,
                                       this->dof_handler,
                                       support_points);

  double l2errorU                  = 0.;
  double l2errorU_boundary         = 0.;
  double l2errorP                  = 0.;
  double total_velocity_divergence = 0.;
  double pressure_integral         = 0;
  double exact_pressure_integral   = 0;
  double volume                    = 0;

  // loop over elements
  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();

  // loop over elements to calculate average pressure
  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          bool cell_is_cut;
          // std::ignore is used because we don't care about what particle cut
          // the cell or the number of particles that cut the cell.
          std::tie(cell_is_cut, std::ignore, std::ignore) = cut_cells_map[cell];
          bool cell_is_inside;
          std::tie(cell_is_inside, std::ignore) = cells_inside_map[cell];
          bool cell_is_overconstrained;
          std::tie(cell_is_overconstrained, std::ignore, std::ignore) =
            overconstrained_fluid_cell_map[cell];

          if ((!cell_is_cut && !cell_is_inside) && !cell_is_overconstrained)
            {
              auto &evaluation_point = this->evaluation_point;
              fe_values.reinit(cell);

              fe_values[pressure].get_function_values(evaluation_point,
                                                      local_pressure_values);
              // Get the exact solution at all gauss points
              l_exact_solution->vector_value_list(
                fe_values.get_quadrature_points(), q_exactSol);


              // Retrieve the effective "connectivity matrix" for this element
              cell->get_dof_indices(local_dof_indices);

              for (unsigned int q = 0; q < n_q_points; q++)
                {
                  pressure_integral +=
                    local_pressure_values[q] * fe_values.JxW(q);
                  exact_pressure_integral +=
                    q_exactSol[q][dim] * fe_values.JxW(q);
                  volume += fe_values.JxW(q);
                }
            }
        }
    }

  pressure_integral =
    Utilities::MPI::sum(pressure_integral, this->mpi_communicator);
  exact_pressure_integral =
    Utilities::MPI::sum(exact_pressure_integral, this->mpi_communicator);
  volume = Utilities::MPI::sum(volume, this->mpi_communicator);


  double average_pressure       = pressure_integral / volume;
  double average_exact_pressure = exact_pressure_integral / volume;
  cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();

  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(local_dof_indices);

          bool cell_is_cut;
          // std::ignore is used because we don't care about what particle cut
          // the cell or the number of particles that cut the cell.
          std::tie(cell_is_cut, std::ignore, std::ignore) = cut_cells_map[cell];

          bool cell_is_overconstrained;
          std::tie(cell_is_overconstrained, std::ignore, std::ignore) =
            overconstrained_fluid_cell_map[cell];

          bool cell_is_inside;
          std::tie(cell_is_inside, std::ignore) = cells_inside_map[cell];
          if (cell->at_boundary() &&
              this->check_existance_of_bc(
                BoundaryConditions::BoundaryType::function_weak))
            {
              for (unsigned int i_bc = 0;
                   i_bc < this->simulation_parameters.boundary_conditions.size;
                   ++i_bc)
                {
                  if (this->simulation_parameters.boundary_conditions
                        .type[i_bc] ==
                      BoundaryConditions::BoundaryType::function_weak)
                    {
                      for (const auto face : cell->face_indices())
                        {
                          if (cell->face(face)->at_boundary())
                            {
                              unsigned int boundary_id =
                                cell->face(face)->boundary_id();
                              if (boundary_id ==
                                  this->simulation_parameters
                                    .boundary_conditions.id[i_bc])
                                {
                                  NavierStokesFunctionDefined<dim> function_v(
                                    &this->simulation_parameters
                                       .boundary_conditions
                                       .bcFunctions[boundary_id]
                                       .u,
                                    &this->simulation_parameters
                                       .boundary_conditions
                                       .bcFunctions[boundary_id]
                                       .v,
                                    &this->simulation_parameters
                                       .boundary_conditions
                                       .bcFunctions[boundary_id]
                                       .w);

                                  fe_face_values.reinit(cell, face);
                                  fe_face_values[velocities]
                                    .get_function_values(
                                      this->present_solution,
                                      local_face_velocity_values);
                                  for (unsigned int q = 0; q < n_q_points_face;
                                       q++)
                                    {
                                      double u_x =
                                        local_face_velocity_values[q][0];
                                      double u_y =
                                        local_face_velocity_values[q][1];

                                      double u_x_a = function_v.value(
                                        fe_face_values.quadrature_point(q), 0);
                                      double u_y_a = function_v.value(
                                        fe_face_values.quadrature_point(q), 1);

                                      l2errorU_boundary +=
                                        ((u_x - u_x_a) * (u_x - u_x_a) +
                                         (u_y - u_y_a) * (u_y - u_y_a)) *
                                        fe_face_values.JxW(q);
                                      if (dim == 3)
                                        {
                                          double u_z =
                                            local_face_velocity_values[q][2];
                                          double u_z_a = function_v.value(
                                            fe_face_values.quadrature_point(q),
                                            2);
                                          l2errorU_boundary +=
                                            (u_z - u_z_a) * (u_z - u_z_a) *
                                            fe_face_values.JxW(q);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

          if (!cell_is_cut && !cell_is_overconstrained)
            {
              auto &evaluation_point = this->evaluation_point;
              auto &present_solution = this->present_solution;
              fe_values.reinit(cell);
              fe_values[velocities].get_function_values(present_solution,
                                                        local_velocity_values);
              fe_values[pressure].get_function_values(present_solution,
                                                      local_pressure_values);
              fe_values[velocities].get_function_gradients(
                evaluation_point, present_velocity_gradients);


              // Retrieve the effective "connectivity matrix" for this element
              cell->get_dof_indices(local_dof_indices);

              // Get the exact solution at all gauss points
              l_exact_solution->vector_value_list(
                fe_values.get_quadrature_points(), q_exactSol);
              for (unsigned int q = 0; q < n_q_points; q++)
                {
                  double present_velocity_divergence =
                    trace(present_velocity_gradients[q]);
                  double mass_source =
                    this->simulation_parameters.source_term
                      ->navier_stokes_source.value(
                        fe_values.get_quadrature_points()[q], dim);


                  // Find the values of x and u_h (the finite element solution)
                  // at the quadrature points
                  double ux_sim   = local_velocity_values[q][0];
                  double ux_exact = q_exactSol[q][0];

                  double uy_sim   = local_velocity_values[q][1];
                  double uy_exact = q_exactSol[q][1];
                  l2errorU += (ux_sim - ux_exact) * (ux_sim - ux_exact) *
                              fe_values.JxW(q);
                  l2errorU += (uy_sim - uy_exact) * (uy_sim - uy_exact) *
                              fe_values.JxW(q);

                  if (dim == 3)
                    {
                      double uz_sim   = local_velocity_values[q][2];
                      double uz_exact = q_exactSol[q][2];
                      l2errorU += (uz_sim - uz_exact) * (uz_sim - uz_exact) *
                                  fe_values.JxW(q);
                    }
                  if (!cell_is_inside)
                    {
                      total_velocity_divergence +=
                        (present_velocity_divergence - mass_source) *
                        fe_values.JxW(q);
                      double p_sim =
                        local_pressure_values[q] - average_pressure;
                      double p_exact =
                        q_exactSol[q][dim] - average_exact_pressure;
                      l2errorP += (p_sim - p_exact) * (p_sim - p_exact) *
                                  fe_values.JxW(q);
                    }
                }
            }
        }
    }
  l2errorU = Utilities::MPI::sum(l2errorU, this->mpi_communicator);
  l2errorU_boundary =
    Utilities::MPI::sum(l2errorU_boundary, this->mpi_communicator);
  l2errorP = Utilities::MPI::sum(l2errorP, this->mpi_communicator);
  total_velocity_divergence =
    Utilities::MPI::sum(total_velocity_divergence, this->mpi_communicator);

  this->pcout << "div u : " << total_velocity_divergence << std::endl;
  if (this->check_existance_of_bc(
        BoundaryConditions::BoundaryType::function_weak))
    {
      this->pcout << "error at the weak BC : " << std::sqrt(l2errorU_boundary)
                  << std::endl;
    }

  return std::make_pair(std::sqrt(l2errorU), std::sqrt(l2errorP));
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::integrate_particles()
{
  particle_residual = 0;


  TimerOutput::Scope t(this->computing_timer, "integrate particles");
  // Integrate the velocity of the particle. If integrate motion is set to
  // true in the parameter this function will also integrate the force to update
  // the velocity. Otherwise the velocity is kept constant

  // To integrate the forces and update the velocity, this function uses the
  // implicit Euler algorithm. To find the force at t+dt the function use the
  // fix point algorithm in parallel to the newton iteration used for the fluid
  // resolution.
  using numbers::PI;
  double dt    = this->simulation_control->get_time_steps_vector()[0];
  double time  = this->simulation_control->get_current_time();
  double alpha = this->simulation_parameters.particlesParameters->alpha;
  this->simulation_parameters.particlesParameters->f_gravity->set_time(time);
  double dr = GridTools::minimal_cell_diameter(*this->triangulation) / sqrt(2);

  const auto rheological_model =
    this->simulation_parameters.physical_properties_manager.get_rheology();
  ib_dem.update_particles(particles, time - dt);
  std::map<field, double> field_values;
  field_values[field::shear_rate]  = 1;
  field_values[field::temperature] = 1;

  // Check if the parameters' combination is compatible. This is temporary and
  // will be moved to a new class that tests the parameter combination in the
  // parameter initialization.

  Assert(!(this->simulation_parameters.physical_properties_manager
             .is_non_newtonian() and
           this->simulation_parameters.particlesParameters
             ->enable_lubrication_force and
           some_particles_are_coupled),
         RequiresConstantViscosity(
           "GLSSharpNavierStokesSolver<dim>::integrate_particles"));

  double h_min =
    dr * this->simulation_parameters.particlesParameters->lubrication_range_min;
  double h_max =
    dr * this->simulation_parameters.particlesParameters->lubrication_range_max;
  // Deactivated the lubrication force if the flui is Non-newtonian.
  if (this->simulation_parameters.particlesParameters
        ->enable_lubrication_force == false)
    {
      h_max = 0;
    }

  particle_residual = 0;
  if (some_particles_are_coupled && time > 0)
    {
      Assert(this->simulation_parameters.physical_properties_manager
               .density_is_constant(),
             RequiresConstantDensity(
               "GLSSharpNavierStokesSolver<dim>::integrate_particles"));
      this->simulation_parameters.physical_properties_manager.get_density(0);
      auto density_model =
        this->simulation_parameters.physical_properties_manager.get_density(0);
      double fluid_density = density_model->value(field_values);
      double viscosity = rheological_model->value(field_values) * fluid_density;

      Vector<double> particles_residual_vect;
      particles_residual_vect.reinit(particles.size());
      ib_dem.integrate_particles_motion(
        dt, h_max, h_min, fluid_density, viscosity);
      unsigned int worst_residual_particle_id = UINT_MAX;

      for (unsigned int p = 0; p < particles.size(); ++p)
        {
          if (particles[p].integrate_motion)
            {
              // calculate the volume of fluid displaced by the particle.
              double volume =
                particles[p].shape->displaced_volume(fluid_density);

              // Transfers the impulsion evaluated in the sub-time-stepping
              // scheme to the particle at the CFD time scale.
              particles[p].impulsion = ib_dem.dem_particles[p].impulsion;
              particles[p].omega_impulsion =
                ib_dem.dem_particles[p].omega_impulsion;
              particles[p].contact_impulsion =
                ib_dem.dem_particles[p].contact_impulsion;
              particles[p].omega_contact_impulsion =
                ib_dem.dem_particles[p].omega_contact_impulsion;
              // Time stepping information
              const auto method =
                this->simulation_control->get_assembly_method();
              std::vector<double> time_steps_vector =
                this->simulation_control->get_time_steps_vector();

              // Vector for the BDF coefficients
              Vector<double> bdf_coefs =
                bdf_coefficients(method, time_steps_vector);

              // Define the residual of the particle dynamics.
              Tensor<1, 3> residual_velocity =
                -(bdf_coefs[0] * particles[p].velocity);

              for (unsigned int i = 1;
                   i < number_of_previous_solutions(method) + 1;
                   ++i)
                {
                  residual_velocity +=
                    -(bdf_coefs[i] * particles[p].previous_velocity[i - 1]);
                }
              residual_velocity +=
                (particles[p].impulsion) / particles[p].mass / dt;
              // Approximate a diagonal Jacobian with a secant method.

              double inverse_of_relaxation_coefficient_velocity =
                -bdf_coefs[0] -
                0.5 * volume * fluid_density / particles[p].mass / dt + DBL_MIN;

              // Evaluate the relaxation parameter using a generalization of the
              // secant method.
              if ((particles[p].velocity - particles[p].velocity_iter).norm() !=
                  0)
                {
                  auto vector_of_velocity_variation =
                    (particles[p].velocity - particles[p].velocity_iter);
                  auto vector_of_residual_variation =
                    (bdf_coefs[0] *
                       (-particles[p].velocity + particles[p].velocity_iter) +
                     (particles[p].impulsion - particles[p].impulsion_iter) /
                       particles[p].mass / dt);
                  double dot_product_of_the_variation_vectors =
                    scalar_product(vector_of_velocity_variation,
                                   vector_of_residual_variation);
                  inverse_of_relaxation_coefficient_velocity =
                    1 / (vector_of_velocity_variation.norm() /
                           vector_of_residual_variation.norm() *
                           dot_product_of_the_variation_vectors /
                           (vector_of_velocity_variation.norm() *
                            vector_of_residual_variation.norm()) +
                         DBL_MIN);
                }
              // Relaxation parameter for the particle dynamics.
              double local_alpha = 1.;
              // We keep in memory the state of the particle before the update
              // in case something went wrong and the particle is now outside of
              // the domain. If this is the case, the particle won't be updated.
              IBParticle<dim> save_particle_state         = particles[p];
              bool            save_particle_state_is_used = false;
              // Define the correction vector.
              Tensor<1, 3> velocity_correction_vector =
                residual_velocity * 1. /
                inverse_of_relaxation_coefficient_velocity;

              // Update the particle state and keep in memory the last iteration
              // information.

              particles[p].velocity_iter = particles[p].velocity;
              particles[p].velocity =
                particles[p].velocity_iter -
                velocity_correction_vector * alpha * local_alpha;
              particles[p].impulsion_iter      = particles[p].impulsion;
              particles[p].previous_d_velocity = velocity_correction_vector;
              particles[p].previous_local_alpha_velocity = local_alpha;


              // If the particles have impacted a wall or another particle, we
              // want to use the sub-time step position. Otherwise, we solve the
              // new position directly with the new velocity found.
              if (particles[p].contact_impulsion.norm() < 1e-12)
                {
                  particles[p].position.clear();
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      for (unsigned int i = 1;
                           i < number_of_previous_solutions(method) + 1;
                           ++i)
                        {
                          double position_update =
                            -bdf_coefs[i] *
                            particles[p].previous_positions[i - 1][d] /
                            bdf_coefs[0];

                          particles[p].set_position(particles[p].position[d] +
                                                      position_update,
                                                    d);
                        }

                      double position_update =
                        particles[p].velocity[d] / bdf_coefs[0];
                      particles[p].set_position(particles[p].position[d] +
                                                  position_update,
                                                d);
                    }
                }
              else
                {
                  Tensor<1, dim> position_update =
                    local_alpha *
                    (ib_dem.dem_particles[p].position - particles[p].position);
                  particles[p].set_position(particles[p].position +
                                            position_update);
                }
              // Check if the particle is in the domain. Throw an error if it's
              // the case.
              try
                {
                  const auto &cell =
                    LetheGridTools::find_cell_around_point_with_tree(
                      this->dof_handler, particles[p].position);
                  (void)cell;
                }
              catch (...)
                {
                  this->pcout
                    << "particle " << p
                    << " is now outside the domain we do not update this particle. Position: "
                    << particles[p].position[0] << ", "
                    << particles[p].position[1] << ", "
                    << particles[p].position[2] << std::endl;
                  save_particle_state_is_used = true;
                }

              // For the rotation velocity : same logic as the velocity.
              auto         inv_inertia = invert(particles[p].inertia);
              Tensor<1, 3> residual_omega =
                -(bdf_coefs[0] * particles[p].omega);
              for (unsigned int i = 1;
                   i < number_of_previous_solutions(method) + 1;
                   ++i)
                {
                  residual_omega +=
                    -(bdf_coefs[i] * particles[p].previous_omega[i - 1]);
                }
              residual_omega +=
                inv_inertia * (particles[p].omega_impulsion) / dt;

              double inverse_of_relaxation_coefficient_omega =
                -bdf_coefs[0] - 0.5 * 2. / 5 * volume * fluid_density *
                                  particles[p].radius * particles[p].radius *
                                  inv_inertia.norm() / dt;
              // Evaluate the relaxation parameter using a generalization of the
              // secant method.
              if ((particles[p].omega - particles[p].omega_iter).norm() != 0)
                {
                  auto vector_of_omega_variation =
                    (particles[p].omega - particles[p].omega_iter);
                  auto vector_of_residual_variation =
                    (bdf_coefs[0] *
                       (-particles[p].omega + particles[p].omega_iter) +
                     (particles[p].omega_impulsion -
                      particles[p].omega_impulsion_iter) *
                       inv_inertia.norm() / dt);
                  double dot_product_of_the_variation_vectors =
                    scalar_product(vector_of_omega_variation,
                                   vector_of_residual_variation);
                  inverse_of_relaxation_coefficient_omega =
                    1 / (vector_of_omega_variation.norm() /
                           vector_of_residual_variation.norm() *
                           dot_product_of_the_variation_vectors /
                           (vector_of_omega_variation.norm() *
                            vector_of_residual_variation.norm()) +
                         DBL_MIN);
                }
              // Define the correction vector.
              Tensor<1, 3> omega_correction_vector =
                residual_omega * 1 / inverse_of_relaxation_coefficient_omega;

              double local_alpha_omega = 1;

              particles[p].omega_iter = particles[p].omega;
              particles[p].omega =
                particles[p].omega_iter -
                omega_correction_vector * alpha * local_alpha_omega;
              particles[p].omega_impulsion_iter = particles[p].omega_impulsion;
              particles[p].previous_d_omega     = omega_correction_vector;
              particles[p].previous_local_alpha_omega = local_alpha_omega;


              // If something went wrong during the update, the particle state
              // would be reversed to its original state here.
              if (save_particle_state_is_used)
                particles[p] = save_particle_state;

              // Evaluate global residual of the particle dynamics.
              double this_particle_residual =
                sqrt(std::pow(residual_velocity.norm(), 2) +
                     std::pow(residual_omega.norm(), 2)) *
                dt;
              // Keep in memory the residual.
              particles[p].residual_velocity = residual_velocity.norm();
              particles[p].residual_omega =
                (particles[p].omega - particles[p].omega_iter).norm();
              particles_residual_vect[p] = this_particle_residual;
              // L_inf of all the particles' residual.
              if (this_particle_residual > particle_residual)
                {
                  particle_residual          = this_particle_residual;
                  worst_residual_particle_id = p;
                }
            }
          else
            {
              if (this->simulation_parameters.particlesParameters
                    ->load_particles_from_file == false)
                {
                  particles[p].f_position->set_time(time);
                  particles[p].f_velocity->set_time(time);
                  particles[p].f_omega->set_time(time);
                  particles[p].f_orientation->set_time(time);

                  particles[p].position[0] =
                    particles[p].f_position->value(particles[p].position, 0);
                  particles[p].position[1] =
                    particles[p].f_position->value(particles[p].position, 1);
                  particles[p].velocity[0] =
                    particles[p].f_velocity->value(particles[p].position, 0);
                  particles[p].velocity[1] =
                    particles[p].f_velocity->value(particles[p].position, 1);
                  particles[p].omega[0] =
                    particles[p].f_omega->value(particles[p].position, 0);
                  particles[p].omega[1] =
                    particles[p].f_omega->value(particles[p].position, 1);
                  particles[p].omega[2] =
                    particles[p].f_omega->value(particles[p].position, 2);
                  particles[p].orientation[0] =
                    particles[p].f_orientation->value(particles[p].position, 0);
                  particles[p].orientation[1] =
                    particles[p].f_orientation->value(particles[p].position, 1);
                  particles[p].orientation[2] =
                    particles[p].f_orientation->value(particles[p].position, 2);
                  if (dim == 3)
                    {
                      particles[p].position[2] =
                        particles[p].f_position->value(particles[p].position,
                                                       2);
                      particles[p].velocity[2] =
                        particles[p].f_velocity->value(particles[p].position,
                                                       2);
                    }
                }
              else
                {
                  particles[p].position[0] += particles[p].velocity[0] * dt;
                  particles[p].position[1] += particles[p].velocity[1] * dt;
                  particles[p].orientation[0] = particles[p].omega[0] * dt;
                  particles[p].orientation[1] = particles[p].omega[1] * dt;
                  particles[p].orientation[2] = particles[p].omega[2] * dt;

                  if (dim == 3)
                    {
                      particles[p].position[2] += particles[p].velocity[2] * dt;
                    }
                }
              if (particles[p].position != particles[p].previous_positions[0] ||
                  particles[p].orientation !=
                    particles[p].previous_orientation[0])
                {
                  particles[p].clear_shape_cache();
                }
              particles[p].set_position(particles[p].position);
              particles[p].set_orientation(particles[p].orientation);
            }
        }


      if (this->simulation_parameters.non_linear_solver.verbosity !=
          Parameters::Verbosity::quiet)
        {
          this->pcout << "L_inf particle residual : " << particle_residual
                      << " particle id : " << worst_residual_particle_id
                      << std::endl;
          this->pcout << "L2 particle residual L2 "
                      << particles_residual_vect.l2_norm() << std::endl;
        }
    }
  else
    {
      // Direct application of the function for the velocity and position if the
      // particle dynamics is not integrated.
      for (unsigned int p = 0; p < particles.size(); ++p)
        {
          if (this->simulation_parameters.particlesParameters
                ->load_particles_from_file == false)
            {
              particles[p].f_position->set_time(time);
              particles[p].f_velocity->set_time(time);
              particles[p].f_omega->set_time(time);
              particles[p].f_orientation->set_time(time);

              particles[p].position[0] =
                particles[p].f_position->value(particles[p].position, 0);
              particles[p].position[1] =
                particles[p].f_position->value(particles[p].position, 1);
              particles[p].velocity[0] =
                particles[p].f_velocity->value(particles[p].position, 0);
              particles[p].velocity[1] =
                particles[p].f_velocity->value(particles[p].position, 1);
              particles[p].omega[0] =
                particles[p].f_omega->value(particles[p].position, 0);
              particles[p].omega[1] =
                particles[p].f_omega->value(particles[p].position, 1);
              particles[p].omega[2] =
                particles[p].f_omega->value(particles[p].position, 2);
              particles[p].orientation[0] =
                particles[p].f_orientation->value(particles[p].position, 0);
              particles[p].orientation[1] =
                particles[p].f_orientation->value(particles[p].position, 1);
              particles[p].orientation[2] =
                particles[p].f_orientation->value(particles[p].position, 2);
              if (dim == 3)
                {
                  particles[p].position[2] =
                    particles[p].f_position->value(particles[p].position, 2);
                  particles[p].velocity[2] =
                    particles[p].f_velocity->value(particles[p].position, 2);
                }
            }
          else
            {
              particles[p].position[0] += particles[p].velocity[0] * dt;
              particles[p].position[1] += particles[p].velocity[1] * dt;
              particles[p].orientation[0] = particles[p].omega[0] * dt;
              particles[p].orientation[1] = particles[p].omega[1] * dt;
              particles[p].orientation[2] = particles[p].omega[2] * dt;

              if (dim == 3)
                {
                  particles[p].position[2] += particles[p].velocity[2] * dt;
                }
            }
          if (particles[p].position != particles[p].previous_positions[0] ||
              particles[p].orientation != particles[p].previous_orientation[0])
            {
              particles[p].clear_shape_cache();
            }

          particles[p].set_position(particles[p].position);
          particles[p].set_orientation(particles[p].orientation);
        }
      particle_residual = 0;
    }
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::Visualization_IB::build_patches(
  std::vector<IBParticle<dim>> particles)
{
  properties_to_write = particles[0].get_properties_name();
  /**
   * A list of field names for all data components stored in patches.
   */
  vector_datasets.clear();
  dataset_names.clear();
  // Defining property field position
  int field_position = 0;
  // Iterating over properties
  for (auto properties_iterator = properties_to_write.begin();
       properties_iterator != properties_to_write.end();
       ++properties_iterator, ++field_position)
    {
      // Get the property field name
      const std::string field_name = properties_iterator->first;

      // Number of components of the corresponding property
      const unsigned components_number = properties_iterator->second;

      // Check to see if the property is a vector
      if (components_number == 3)
        {
          vector_datasets.push_back(std::make_tuple(
            field_position,
            field_position + components_number - 1,
            field_name,
            DataComponentInterpretation::component_is_part_of_vector));
        }
      dataset_names.push_back(field_name);
    }

  // Building the patch data
  patches.resize(particles.size());

  // Looping over particle to get the properties from the particle_handler
  for (unsigned int p = 0; p < particles.size(); ++p)
    {
      // Particle location
      patches[p].vertices[0] = particles[p].position;
      patches[p].patch_index = p;
      patches[p].data.reinit(particles[p].get_number_properties(), 1);

      // ID and other properties

      // Calculating force for visualization
      auto particle_properties = particles[p].get_properties();

      for (unsigned int property_index = 0;
           property_index < particles[p].get_number_properties();
           ++property_index)
        patches[p].data(property_index, 0) =
          particle_properties[property_index];
    }
}

template <int dim>
const std::vector<DataOutBase::Patch<0, dim>> &
GLSSharpNavierStokesSolver<dim>::Visualization_IB::get_patches() const
{
  return patches;
}

template <int dim>
std::vector<std::string>
GLSSharpNavierStokesSolver<dim>::Visualization_IB::get_dataset_names() const
{
  return dataset_names;
}

template <int dim>
std::vector<
  std::tuple<unsigned int,
             unsigned int,
             std::string,
             DataComponentInterpretation::DataComponentInterpretation>>
GLSSharpNavierStokesSolver<dim>::Visualization_IB::get_nonscalar_data_ranges()
  const
{
  return vector_datasets;
}

template <int dim>
GLSSharpNavierStokesSolver<dim>::Visualization_IB::~Visualization_IB()
{}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::finish_time_step_particles()
{
  // Store information about the particle used for the integration and print the
  // results if requested.

  const std::string folder =
    this->simulation_parameters.simulation_control.output_folder;
  const std::string particles_solution_name =
    this->simulation_parameters.particlesParameters->ib_particles_pvd_file;
  const unsigned int iter = this->simulation_control->get_step_number();
  const double       time = this->simulation_control->get_current_time();
  const unsigned int group_files =
    this->simulation_parameters.simulation_control.group_files;

  // If the processor id is id=0 we write the particles pvd.
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      Visualization_IB ib_particles_data;
      ib_particles_data.build_patches(particles);
      write_vtu_and_pvd<0, dim>(ib_particles_pvdhandler,
                                ib_particles_data,
                                folder,
                                particles_solution_name,
                                time,
                                iter,
                                group_files,
                                this->mpi_communicator);
    }
  else
    {
      // If the processor id is not id=0 we add an empty particle vector.
      Visualization_IB             ib_particles_data;
      std::vector<IBParticle<dim>> empty_particle_vector(0);
      ib_particles_data.build_patches(empty_particle_vector);
      write_vtu_and_pvd<0, dim>(ib_particles_pvdhandler,
                                ib_particles_data,
                                folder,
                                particles_solution_name,
                                time,
                                iter,
                                group_files,
                                this->mpi_communicator);
    }


  table_all_p.clear();
  for (unsigned int p = 0; p < particles.size(); ++p)
    {
      for (unsigned int i = particles[p].previous_velocity.size() - 1; i > 0;
           --i)
        {
          particles[p].previous_positions[i] =
            particles[p].previous_positions[i - 1];
          particles[p].previous_velocity[i] =
            particles[p].previous_velocity[i - 1];
          particles[p].previous_omega[i] = particles[p].previous_omega[i - 1];
        }

      particles[p].previous_positions[0] = particles[p].position;
      particles[p].previous_velocity[0]  = particles[p].velocity;
      particles[p].previous_omega[0]     = particles[p].omega;

      particles[p].previous_fluid_forces = particles[p].fluid_forces;
      particles[p].previous_fluid_viscous_forces =
        particles[p].fluid_viscous_forces;
      particles[p].previous_fluid_pressure_forces =
        particles[p].fluid_pressure_forces;
      particles[p].previous_fluid_torque = particles[p].fluid_torque;

      particles[p].velocity_iter        = particles[p].velocity;
      particles[p].impulsion_iter       = particles[p].impulsion;
      particles[p].omega_iter           = particles[p].omega;
      particles[p].omega_impulsion_iter = particles[p].omega_impulsion;
      particles[p].residual_velocity    = DBL_MAX;
      particles[p].residual_omega       = DBL_MAX;



      if (some_particles_are_coupled &&
          this->simulation_parameters.particlesParameters->print_dem)
        {
          this->pcout << "particle " << p << " position "
                      << particles[p].position << std::endl;
          if (dim == 2)
            {
              this->pcout << "particle " << p << " velocity "
                          << tensor_nd_to_2d(particles[p].velocity)
                          << std::endl;
            }
          else
            {
              this->pcout << "particle " << p << " velocity "
                          << particles[p].velocity << std::endl;
            }
        }
      table_p[p].add_value("particle_ID", p);
      table_all_p.add_value("particle_ID", p);
      if (this->simulation_parameters.simulation_control.method !=
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        table_p[p].add_value("time",
                             this->simulation_control->get_current_time());
      table_all_p.add_value("time",
                            this->simulation_control->get_current_time());
      if (dim == 3)
        {
          table_p[p].add_value("T_x", particles[p].fluid_torque[0]);
          table_p[p].set_precision(
            "T_x",
            this->simulation_parameters.simulation_control.log_precision);
          table_all_p.add_value("T_x", particles[p].fluid_torque[0]);
          table_all_p.set_precision(
            "T_x",
            this->simulation_parameters.simulation_control.log_precision);
          if (some_particles_are_coupled)
            {
              table_p[p].add_value("omega_x", particles[p].omega[0]);
              table_p[p].set_precision(
                "omega_x",
                this->simulation_parameters.simulation_control.log_precision);
              table_all_p.add_value("omega_x", particles[p].omega[0]);
              table_all_p.set_precision(
                "omega_x",
                this->simulation_parameters.simulation_control.log_precision);
            }

          table_p[p].add_value("T_y", particles[p].fluid_torque[1]);
          table_p[p].set_precision(
            "T_y",
            this->simulation_parameters.simulation_control.log_precision);
          table_all_p.add_value("T_y", particles[p].fluid_torque[1]);
          table_all_p.set_precision(
            "T_y",
            this->simulation_parameters.simulation_control.log_precision);
          if (some_particles_are_coupled)
            {
              table_p[p].add_value("omega_y", particles[p].omega[1]);
              table_p[p].set_precision(
                "omega_y",
                this->simulation_parameters.simulation_control.log_precision);
              table_all_p.add_value("omega_y", particles[p].omega[1]);
              table_all_p.set_precision(
                "omega_y",
                this->simulation_parameters.simulation_control.log_precision);
            }
        }

      table_p[p].add_value("T_z", particles[p].fluid_torque[2]);
      table_p[p].set_precision(
        "T_z", this->simulation_parameters.simulation_control.log_precision);
      table_all_p.add_value("T_z", particles[p].fluid_torque[2]);
      table_all_p.set_precision(
        "T_z", this->simulation_parameters.simulation_control.log_precision);
      if (some_particles_are_coupled)
        {
          table_p[p].add_value("omega_z", particles[p].omega[2]);
          table_p[p].set_precision(
            "omega_z",
            this->simulation_parameters.simulation_control.log_precision);
          table_all_p.add_value("omega_z", particles[p].omega[2]);
          table_all_p.set_precision(
            "omega_z",
            this->simulation_parameters.simulation_control.log_precision);
        }

      table_p[p].add_value("f_x", particles[p].fluid_forces[0]);
      table_all_p.add_value("f_x", particles[p].fluid_forces[0]);
      table_p[p].set_precision(
        "f_x", this->simulation_parameters.simulation_control.log_precision);
      table_all_p.set_precision(
        "f_x", this->simulation_parameters.simulation_control.log_precision);
      if (some_particles_are_coupled)
        {
          table_p[p].add_value("v_x", particles[p].velocity[0]);
          table_p[p].add_value("p_x", particles[p].position[0]);
          table_all_p.add_value("v_x", particles[p].velocity[0]);
          table_all_p.add_value("p_x", particles[p].position[0]);
          table_p[p].set_precision(
            "v_x",
            this->simulation_parameters.simulation_control.log_precision);

          table_p[p].set_precision(
            "p_x",
            this->simulation_parameters.simulation_control.log_precision);
          table_all_p.set_precision(
            "v_x",
            this->simulation_parameters.simulation_control.log_precision);
          table_all_p.set_precision(
            "p_x",
            this->simulation_parameters.simulation_control.log_precision);
        }

      table_p[p].add_value("f_y", particles[p].fluid_forces[1]);
      table_p[p].set_precision(
        "f_y", this->simulation_parameters.simulation_control.log_precision);
      table_all_p.add_value("f_y", particles[p].fluid_forces[1]);
      table_all_p.set_precision(
        "f_y", this->simulation_parameters.simulation_control.log_precision);
      if (some_particles_are_coupled)
        {
          table_p[p].add_value("v_y", particles[p].velocity[1]);
          table_p[p].add_value("p_y", particles[p].position[1]);
          table_all_p.add_value("v_y", particles[p].velocity[1]);
          table_all_p.add_value("p_y", particles[p].position[1]);
          table_p[p].set_precision(
            "v_y",
            this->simulation_parameters.simulation_control.log_precision);
          table_p[p].set_precision(
            "p_y",
            this->simulation_parameters.simulation_control.log_precision);
          table_all_p.set_precision(
            "v_y",
            this->simulation_parameters.simulation_control.log_precision);
          table_all_p.set_precision(
            "p_y",
            this->simulation_parameters.simulation_control.log_precision);
        }

      if (dim == 3)
        {
          table_p[p].add_value("f_z", particles[p].fluid_forces[2]);
          table_p[p].set_precision(
            "f_z",
            this->simulation_parameters.simulation_control.log_precision);
          table_all_p.add_value("f_z", particles[p].fluid_forces[2]);
          table_all_p.set_precision(
            "f_z",
            this->simulation_parameters.simulation_control.log_precision);
          if (some_particles_are_coupled)
            {
              table_p[p].add_value("v_z", particles[p].velocity[2]);
              table_p[p].add_value("p_z", particles[p].position[2]);
              table_all_p.add_value("v_z", particles[p].velocity[2]);
              table_all_p.add_value("p_z", particles[p].position[2]);
              table_p[p].set_precision(
                "v_z",
                this->simulation_parameters.simulation_control.log_precision);
              table_p[p].set_precision(
                "p_z",
                this->simulation_parameters.simulation_control.log_precision);
              table_all_p.set_precision(
                "v_z",
                this->simulation_parameters.simulation_control.log_precision);
              table_all_p.set_precision(
                "p_z",
                this->simulation_parameters.simulation_control.log_precision);
            }
        }

      table_p[p].add_value("f_xv", particles[p].fluid_viscous_forces[0]);
      table_all_p.add_value("f_xv", particles[p].fluid_viscous_forces[0]);
      table_p[p].add_value("f_yv", particles[p].fluid_viscous_forces[1]);
      table_all_p.add_value("f_yv", particles[p].fluid_viscous_forces[1]);

      table_p[p].set_precision(
        "f_xv", this->simulation_parameters.simulation_control.log_precision);
      table_p[p].set_precision(
        "f_yv", this->simulation_parameters.simulation_control.log_precision);
      table_all_p.set_precision(
        "f_xv", this->simulation_parameters.simulation_control.log_precision);
      table_all_p.set_precision(
        "f_yv", this->simulation_parameters.simulation_control.log_precision);
      if (dim == 3)
        {
          table_p[p].add_value("f_zv", particles[p].fluid_viscous_forces[2]);
          table_all_p.add_value("f_zv", particles[p].fluid_viscous_forces[2]);
          table_p[p].set_precision(
            "f_zv",
            this->simulation_parameters.simulation_control.log_precision);
          table_all_p.set_precision(
            "f_zv",
            this->simulation_parameters.simulation_control.log_precision);
        }

      table_p[p].add_value("f_xp", particles[p].fluid_pressure_forces[0]);
      table_all_p.add_value("f_xp", particles[p].fluid_pressure_forces[0]);
      table_p[p].add_value("f_yp", particles[p].fluid_pressure_forces[1]);
      table_all_p.add_value("f_yp", particles[p].fluid_pressure_forces[1]);
      table_p[p].set_precision(
        "f_xp", this->simulation_parameters.simulation_control.log_precision);
      table_p[p].set_precision(
        "f_yp", this->simulation_parameters.simulation_control.log_precision);
      table_all_p.set_precision(
        "f_xp", this->simulation_parameters.simulation_control.log_precision);
      table_all_p.set_precision(
        "f_yp", this->simulation_parameters.simulation_control.log_precision);
      if (dim == 3)
        {
          table_p[p].add_value("f_zp", particles[p].fluid_pressure_forces[2]);
          table_all_p.add_value("f_zp", particles[p].fluid_pressure_forces[2]);
          table_p[p].set_precision(
            "f_zp",
            this->simulation_parameters.simulation_control.log_precision);
          table_all_p.set_precision(
            "f_zp",
            this->simulation_parameters.simulation_control.log_precision);
        }
    }
  if (this->this_mpi_process == 0)
    {
      if (this->simulation_parameters.forces_parameters.verbosity ==
          Parameters::Verbosity::verbose)
        {
          std::cout << "+------------------------------------------+"
                    << std::endl;
          std::cout << "|  Force  summary particles  "
                    << "              |" << std::endl;
          std::cout << "+------------------------------------------+"
                    << std::endl;
          table_all_p.write_text(std::cout);
        }
    }
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::sharp_edge()
{
  // This function defines an Immersed Boundary based on the sharp edge method
  // on a solid of dim=2 or dim=3

  TimerOutput::Scope t(this->computing_timer, "assemble_sharp");
  using numbers::PI;
  Point<dim>                                                  center_immersed;
  Point<dim>                                                  pressure_bridge;
  std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors;
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
    active_neighbors_set;
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
                                   active_neighbors_2;
  const FEValuesExtractors::Scalar pressure(dim);

  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  // Define a map to all dofs and their support points
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(*this->mapping,
                                       this->dof_handler,
                                       support_points);

  // Initalize fe value objects in order to do calculation with it later
  QGauss<dim>        q_formula(this->number_quadrature_points);
  FEValues<dim>      fe_values(*this->fe,
                          q_formula,
                          update_quadrature_points | update_JxW_values);
  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;

  int    order = this->simulation_parameters.particlesParameters->order;
  double length_ratio =
    this->simulation_parameters.particlesParameters->length_ratio;



  IBStencil<dim>      stencil;
  std::vector<double> ib_coef = stencil.coefficients(order, length_ratio);

  unsigned int n_q_points = q_formula.size();

  // Define multiple local_dof_indices one for the cell iterator one for the
  // cell with the second point for the sharp edge stencil and one for
  // manipulation on the neighbour’s cell.

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_2(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_3(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_4(dofs_per_cell);
  std::set<unsigned int>               clear_line;

  // Define minimal cell length
  double dr = GridTools::minimal_cell_diameter(*this->triangulation) / sqrt(2);

  // Define cell iterator
  const auto &cell_iterator = this->dof_handler.active_cell_iterators();
  double      dt            = time_steps_vector[0];
  if (Parameters::SimulationControl::TimeSteppingMethod::steady ==
      this->simulation_parameters.simulation_control.method)
    dt = 1;

  if (this->simulation_parameters.particlesParameters
        ->assemble_navier_stokes_inside == true)
    {
      // impose pressure reference in each of the particle
      for (unsigned int p = 0; p < particles.size(); ++p)
        {
          Point<dim> pressure_reference_location =
            particles[p].pressure_location + particles[p].position;


          const auto &cell = LetheGridTools::find_cell_around_point_with_tree(
            this->dof_handler, pressure_reference_location);

          if (cell->is_locally_owned())
            {
              cell->get_dof_indices(local_dof_indices);
              double sum_line = 0;
              double volume   = 0;
              fe_values.reinit(cell);
              std::vector<int> set_pressure_cell;
              set_pressure_cell.resize(particles.size());

              // Define the order of magnitude for the stencil.
              for (unsigned int qf = 0; qf < n_q_points; ++qf)
                volume += fe_values.JxW(qf);

              sum_line = volume / dt;
              // Clear the line in the matrix
              unsigned int inside_index = local_dof_indices[dim];
              // Check on which DOF of the cell to impose the pressure. If the
              // dof is on a hanging node, it is already constrained and the
              // pressure cannot be imposed there. So we just go to the next
              // pressure DOF of the cell.

              for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
                {
                  const unsigned int component_i =
                    this->fe->system_to_component_index(i).first;
                  if (this->zero_constraints.is_constrained(
                        local_dof_indices[i]) == false &&
                      this->locally_owned_dofs.is_element(
                        local_dof_indices[i]) &&
                      component_i == dim)
                    {
                      inside_index = local_dof_indices[i];
                      break;
                    }
                }

              this->system_matrix.clear_row(inside_index);
              // this->system_matrix.clear_row(inside_index)
              // is not reliable on edge case

              // Set the new equation for the first pressure dofs of the
              // cell. this is the new reference pressure inside a
              // particle

              this->system_matrix.add(inside_index, inside_index, sum_line);
              this->system_rhs(inside_index) =
                0 - this->evaluation_point(inside_index) * sum_line;
            }
        }
    }

  ib_done.clear();
  // Loop on all the cell to define if the sharp edge cut them
  for (const auto &cell_cut : cell_iterator)
    {
      if (cell_cut->is_locally_owned() || cell_cut->is_ghost())
        {
          // Check if the cell_cut is cut or not by the IB and what the particle
          // the cut the cell_cut. If the particle is cut
          bool cell_is_cut;
          // The id of the particle that cut the cell_cut. Returns 0 if the
          // cell_cut is not cut.We also check the number of particles that cut
          // the cell_cut. If multiple particles cut the cell_cut, the dummy
          // dofs of pressure will be treated differently to avoid
          // self-reference.
          unsigned int ib_particle_id;
          unsigned int count_particles;
          std::tie(cell_is_cut, ib_particle_id, count_particles) =
            cut_cells_map[cell_cut];
          bool cell_is_overconstrained;
          std::tie(cell_is_overconstrained, std::ignore, std::ignore) =
            overconstrained_fluid_cell_map[cell_cut];
          if (cell_is_overconstrained)
            std::tie(cell_is_overconstrained, ib_particle_id, std::ignore) =
              overconstrained_fluid_cell_map[cell_cut];

          if (cell_is_cut || cell_is_overconstrained)
            {
              double sum_line = 0;
              fe_values.reinit(cell_cut);

              double volume = 0;
              // Define the order of magnitude for the stencil.
              for (unsigned int qf = 0; qf < n_q_points; ++qf)
                volume += fe_values.JxW(qf);

              sum_line = volume / dt;
              cell_cut->get_dof_indices(local_dof_indices);
              // If we are here, the cell_cut is cut by the IB.
              // Loops on the dof that represents the velocity  component
              // and pressure separately
              for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
                {
                  unsigned int global_index_overwrite = local_dof_indices[i];
                  if (ib_done[global_index_overwrite].first == false)
                    {
                      const unsigned int component_i =
                        this->fe->system_to_component_index(i).first;
                      bool dof_is_inside =
                        particles[ib_particle_id].get_levelset(
                          support_points[local_dof_indices[i]], cell_cut) <= 0;

                      // If multiple particles cut the cell_cut, we treat the
                      // dof of pressure as a dummy dof. We don't use them to
                      // set the boundary condition for the Poisson problem
                      // inside the particle.
                      bool use_ib_for_pressure =
                        (dof_is_inside) && (component_i == dim) &&
                        (this->simulation_parameters.particlesParameters
                           ->assemble_navier_stokes_inside == false) &&
                        count_particles < 2;

                      // Reapply the zero and nonzero constraints.
                      if ((this->zero_constraints.is_constrained(
                             local_dof_indices[i]) ||
                           this->nonzero_constraints.is_constrained(
                             local_dof_indices[i])) &&
                          this->locally_owned_dofs.is_element(
                            global_index_overwrite))
                        {
                          if (this->zero_constraints.is_constrained(
                                local_dof_indices[i]))
                            {
                              // Clear the line if there is something on it
                              this->system_matrix.clear_row(
                                global_index_overwrite);
                              // Get the constraint equations
                              auto local_entries =
                                *this->zero_constraints.get_constraint_entries(
                                  local_dof_indices[i]);
                              double interpolation = 0;
                              // Write the equation
                              for (unsigned int j = 0; j < local_entries.size();
                                   ++j)
                                {
                                  unsigned int col = local_entries[j].first;
                                  double entries   = local_entries[j].second;

                                  interpolation +=
                                    this->evaluation_point(col) * entries;

                                  try
                                    {
                                      this->system_matrix.add(
                                        local_dof_indices[i],
                                        col,
                                        entries * sum_line);
                                    }
                                  catch (...)
                                    {
                                      //  If we are here, an error happens
                                      //  when trying to fill the line in
                                      //  the matrix.
                                      // For example, this can occur if a
                                      // particle is close to a wall and we
                                      // are trying to impose the equation
                                      // on a DOF that as a boundary
                                      // condition applied to it. As such,
                                      // we discard these errors.
                                    }
                                }
                              this->system_matrix.add(local_dof_indices[i],
                                                      local_dof_indices[i],
                                                      sum_line);

                              // Write the RHS
                              this->system_rhs(local_dof_indices[i]) =
                                -this->evaluation_point(local_dof_indices[i]) *
                                  sum_line +
                                interpolation * sum_line +
                                this->zero_constraints.get_inhomogeneity(
                                  local_dof_indices[i]) *
                                  sum_line;
                            }
                          if (this->nonzero_constraints.is_constrained(
                                local_dof_indices[i]))
                            {
                              // Clear the line if there is something on it
                              this->system_matrix.clear_row(
                                global_index_overwrite);
                              // Get the constraint equations

                              auto local_entries_non_zero =
                                *this->nonzero_constraints
                                   .get_constraint_entries(
                                     local_dof_indices[i]);

                              double interpolation = 0;
                              // Write the equation
                              for (unsigned int j = 0;
                                   j < local_entries_non_zero.size();
                                   ++j)
                                {
                                  unsigned int col_non_zero =
                                    local_entries_non_zero[j].first;
                                  double entries_non_zero =
                                    local_entries_non_zero[j].second;

                                  interpolation +=
                                    this->evaluation_point(col_non_zero) *
                                    entries_non_zero;
                                  try
                                    {
                                      this->system_matrix.add(
                                        local_dof_indices[i],
                                        col_non_zero,
                                        entries_non_zero * sum_line);
                                    }
                                  catch (...)
                                    {
                                      //  If we are here, an error happens
                                      //  when trying to fill the line in
                                      //  the matrix.
                                      // For example, this can occur if a
                                      // particle is close to a wall and we
                                      // are trying to impose the equation
                                      // on a DOF that as a boundary
                                      // condition applied to it. As such,
                                      // we discard these errors.
                                    }
                                }
                              this->system_matrix.add(local_dof_indices[i],
                                                      local_dof_indices[i],
                                                      sum_line);

                              // Write the RHS
                              this->system_rhs(local_dof_indices[i]) =
                                -this->evaluation_point(local_dof_indices[i]) *
                                  sum_line +
                                interpolation * sum_line +
                                this->nonzero_constraints.get_inhomogeneity(
                                  local_dof_indices[i]) *
                                  sum_line;
                            }
                        }

                      // Check if the DOfs is owned and if it's not a hanging
                      // node.
                      if (((component_i < dim) || use_ib_for_pressure) &&
                          this->locally_owned_dofs.is_element(
                            global_index_overwrite) &&
                          (this->zero_constraints.is_constrained(
                             local_dof_indices[i]) ||
                           this->nonzero_constraints.is_constrained(
                             local_dof_indices[i])) == false)
                        {
                          // We are working on the DOFs of cells cut by a
                          // particle. We clear the current equation associated
                          // with the current DOF and replace it with the sharp
                          // IB constraints.

                          // Clear the current line of this dof
                          this->system_matrix.clear_row(global_index_overwrite);

                          // Define the points for the IB stencil based on the
                          // order, particle position, and DOF position. The
                          // definition of the output variable "point" changes
                          // depending on the order. In the case of stencil
                          // orders 1 to 4, the variable point returns the
                          // position of the DOF directly. In the case of higher
                          // order stencil (5 or more), it returns the position
                          // of the point that is on the IB. This is because
                          // stencil orders higher than four are not
                          // implemented. The function extrapolates the element
                          // at the particle's surface in these cases. To do so,
                          // we use the point at the surface of the particle.
                          // The variable "interpolation points" return the
                          // points used to define the cell_cut used for the
                          // stencil definition and the locations of the points
                          // used in the stencil calculation.

                          auto [point, interpolation_points] =
                            stencil.support_points_for_interpolation(
                              order,
                              length_ratio,
                              particles[ib_particle_id],
                              support_points[local_dof_indices[i]],
                              cell_cut);
                          // Find the cell_cut used for the stencil definition.
                          auto point_to_find_cell =
                            stencil.point_for_cell_detection(
                              particles[ib_particle_id],
                              support_points[local_dof_indices[i]],
                              cell_cut);
                          typename DoFHandler<dim>::active_cell_iterator
                               stencil_cell;
                          bool particle_close_to_wall = false;
                          (void)particle_close_to_wall;
                          try
                            {
                              stencil_cell = LetheGridTools::
                                find_cell_around_point_with_neighbors<dim>(
                                  this->dof_handler,
                                  vertices_to_cell,
                                  cell_cut,
                                  point_to_find_cell);
                            }
                          catch (...)
                            {
                              // If we are here, the DOF is on a boundary.
                              particle_close_to_wall = true;
                              stencil_cell           = cell_cut;
                              // If a boundary condition is already applied to
                              // this DOF we skip it otherwise we impose a value
                              // base on the velocity of the particle.
                              if (this->zero_constraints.is_constrained(
                                    global_index_overwrite) ||
                                  this->nonzero_constraints.is_constrained(
                                    global_index_overwrite))
                                {
                                  continue;
                                }
                            }

                          stencil_cell->get_dof_indices(local_dof_indices_2);
                          ib_done[global_index_overwrite] =
                            std::make_pair(true, stencil_cell);

                          bool skip_stencil = false;

                          // Check if the DOF intersect the IB
                          bool dof_on_ib = false;

                          // Check if this dof is a dummy dof or directly on IB
                          // and Check if the point used to define the cell_cut
                          // used for the definition of the stencil
                          // ("stencil_cell") is on a face between the cell_cut
                          // that is cut ("cell_cut") and the "stencil_cell".
                          bool point_in_cell = cell_cut->point_inside(
                            interpolation_points
                              [stencil.number_of_interpolation_support_points(
                                 order) -
                               1]);

                          bool         dof_is_dummy = false;
                          bool         cell2_is_cut;
                          unsigned int ib_particle_id_2;
                          std::tie(cell2_is_cut,
                                   ib_particle_id_2,
                                   std::ignore) = cut_cells_map[stencil_cell];
                          bool cell2_is_overconstrained;
                          std::tie(cell2_is_overconstrained,
                                   std::ignore,
                                   std::ignore) =
                            overconstrained_fluid_cell_map[stencil_cell];
                          if (cell2_is_cut || point_in_cell)
                            {
                              dof_is_dummy = true;
                            }

                          if (dof_is_dummy || use_ib_for_pressure ||
                              cell2_is_overconstrained ||
                              particle_close_to_wall)
                            {
                              // Give the DOF an approximated value. This help
                              // with pressure shock when the DOF passes from
                              // part of the boundary to the fluid.

                              this->system_matrix.add(global_index_overwrite,
                                                      global_index_overwrite,
                                                      sum_line);
                              skip_stencil = true;

                              // Tolerance to define an intersection of
                              // the DOF and IB
                              if (abs(particles[ib_particle_id].get_levelset(
                                    support_points[local_dof_indices[i]],
                                    cell_cut)) <= 1e-12 * dr)
                                {
                                  dof_on_ib = true;
                                }
                            }
                          // Define the variable used for the
                          // extrapolation of the actual solution at the
                          // boundaries in order to define the correction

                          // Define the unit cell_cut points for the points
                          // used in the stencil.
                          std::vector<Point<dim>>
                            unite_cell_interpolation_points(ib_coef.size());
                          unite_cell_interpolation_points[0] =
                            this->mapping->transform_real_to_unit_cell(
                              stencil_cell, point);
                          for (unsigned int j = 1; j < ib_coef.size(); ++j)
                            {
                              unite_cell_interpolation_points[j] =
                                this->mapping->transform_real_to_unit_cell(
                                  stencil_cell, interpolation_points[j - 1]);
                            }

                          std::vector<double> local_interp_sol(ib_coef.size());

                          // Define the new matrix entry for this dof
                          if (skip_stencil == false)
                            {
                              for (unsigned int j = 0;
                                   j < local_dof_indices_2.size();
                                   ++j)
                                {
                                  const unsigned int component_j =
                                    this->fe->system_to_component_index(j)
                                      .first;
                                  if (component_j == component_i)
                                    {
                                      //  Define the solution at each point used
                                      //  for the stencil and applied the
                                      //  stencil for the specific DOF. For
                                      //  stencils of order 4 or higher, the
                                      //  stencil is defined through direct
                                      //  extrapolation of the cell_cut. This
                                      //  can only be done when using a
                                      //  structured mesh as this required a
                                      //  mapping of a point outside of a
                                      //  cell_cut.

                                      // Define the local matrix entries of this
                                      // DOF based on its contribution of each
                                      // of the points used in the stencil
                                      // definition and the coefficient
                                      // associated with this point. This loop
                                      // defined the current solution at the
                                      // boundary using the same stencil. This
                                      // is needed to define the residual.
                                      double local_matrix_entry = 0;
                                      for (unsigned int k = 0;
                                           k < ib_coef.size();
                                           ++k)
                                        {
                                          local_matrix_entry +=
                                            this->fe->shape_value(
                                              j,
                                              unite_cell_interpolation_points
                                                [k]) *
                                            ib_coef[k];
                                          local_interp_sol[k] +=
                                            this->fe->shape_value(
                                              j,
                                              unite_cell_interpolation_points
                                                [k]) *
                                            this->evaluation_point(
                                              local_dof_indices_2[j]);
                                        }
                                      // update the matrix.
                                      try
                                        {
                                          this->system_matrix.add(
                                            global_index_overwrite,
                                            local_dof_indices_2[j],
                                            local_matrix_entry * sum_line);
                                        }
                                      catch (...)
                                        {
                                          //  If we are here, an error happens
                                          //  when trying to fill the line in
                                          //  the matrix.
                                          // For example, this can occur if a
                                          // particle is close to a wall and we
                                          // are trying to impose the equation
                                          // on a DOF that as a boundary
                                          // condition applied to it. As such,
                                          // we discard these errors.
                                        }
                                    }
                                }
                            }

                          // Define the RHS of the equation.

                          double rhs_add = 0;
                          // Different boundary conditions depending
                          // on the component index of the DOF and
                          // the dimension.
                          double v_ib = stencil.ib_velocity(
                            particles[ib_particle_id],
                            support_points[local_dof_indices[i]],
                            component_i,
                            cell_cut);

                          //  If the pressure is imposed trough IB inside the
                          //  particle we use an approximation of the pressure
                          //  outside of the IB in this cell_cut.
                          if (component_i == dim)
                            {
                              for (unsigned int k = 0;
                                   k < local_dof_indices.size();
                                   ++k)
                                {
                                  // Check if the dof is inside or outside of
                                  // the particle.
                                  bool dof_is_inside_p =
                                    particles[ib_particle_id].get_levelset(
                                      support_points[local_dof_indices[k]],
                                      cell_cut) <= 0;
                                  const unsigned int component_k =
                                    this->fe->system_to_component_index(k)
                                      .first;
                                  if (component_k == dim &&
                                      dof_is_inside_p == false)
                                    {
                                      v_ib = this->evaluation_point(
                                        local_dof_indices[k]);
                                      try
                                        {
                                          this->system_matrix.add(
                                            global_index_overwrite,
                                            local_dof_indices[k],
                                            -1 * sum_line);
                                        }
                                      catch (...)
                                        {
                                          //  If we are here, an error happens
                                          //  when trying to fill the line in
                                          //  the matrix.
                                          // For example, this can occur if a
                                          // particle is close to a wall and we
                                          // are trying to impose the equation
                                          // on a DOF that as a boundary
                                          // condition applied to it. As such,
                                          // we discard these errors.
                                        }
                                      break;
                                    }
                                }
                            }

                          for (unsigned int k = 0; k < ib_coef.size(); ++k)
                            {
                              rhs_add +=
                                -local_interp_sol[k] * ib_coef[k] * sum_line;
                            }
                          // The rhs contribution where the component would
                          // correspond to a pressure dof are divided by the
                          // pressure_scaling_factor. This is because when the
                          // pressure dofs are rescaled in the newton
                          // correction, there is no distinction between
                          // dummy and non-dummy pressure dofs. Since this
                          // case is only applicable to this solver, the
                          // choice was made to modify the rhs itself.
                          const double pressure_scaling_factor =
                            this->simulation_parameters.stabilization
                              .pressure_scaling_factor;
                          this->system_rhs(global_index_overwrite) =
                            component_i == dim ? (v_ib * sum_line + rhs_add) /
                                                   pressure_scaling_factor :
                                                 v_ib * sum_line + rhs_add;

                          if (dof_on_ib)
                            // Dof is on the immersed boundary
                            this->system_rhs(global_index_overwrite) =
                              component_i == dim ?
                                (v_ib * sum_line - this->evaluation_point(
                                                     global_index_overwrite) *
                                                     sum_line) /
                                  pressure_scaling_factor :
                                (v_ib * sum_line - this->evaluation_point(
                                                     global_index_overwrite) *
                                                     sum_line);

                          if (skip_stencil && dof_on_ib == false)
                            // Impose the value of the dummy dof. This help
                            // with pressure variation when the IB is
                            // moving.
                            this->system_rhs(global_index_overwrite) =
                              component_i == dim ?
                                (sum_line * v_ib - this->evaluation_point(
                                                     global_index_overwrite) *
                                                     sum_line) /
                                  pressure_scaling_factor :
                                sum_line * v_ib - this->evaluation_point(
                                                    global_index_overwrite) *
                                                    sum_line;
                        }

                      if (component_i == dim &&
                          this->locally_owned_dofs.is_element(
                            global_index_overwrite))
                        {
                          // Applied equation on dof that have no equation
                          // defined for them. those DOF become Dummy dof. This
                          // is useful for high order cells or when a dof is
                          // only element of cells that are cut.
                          unsigned int global_index_overwrite =
                            local_dof_indices[i];
                          bool dummy_dof = true;

                          // To check if the pressure dof is a dummy. first
                          // check if the matrix entry is close to 0.
                          if (abs(this->system_matrix.el(
                                global_index_overwrite,
                                global_index_overwrite)) <= 1e-16 * dr)
                            {
                              // If the matrix entry on the diagonal of this DOF
                              // is close to zero, check if all the cells close
                              // are cut. If it's the case, the DOF is a dummy
                              // DOF.
                              active_neighbors_set =
                                LetheGridTools::find_cells_around_cell<dim>(
                                  vertices_to_cell, cell_cut);
                              for (unsigned int m = 0;
                                   m < active_neighbors_set.size();
                                   m++)
                                {
                                  const auto &neighbor_cell =
                                    active_neighbors_set[m];
                                  neighbor_cell->get_dof_indices(
                                    local_dof_indices_3);
                                  for (unsigned int o = 0;
                                       o < local_dof_indices_3.size();
                                       ++o)
                                    {
                                      if (global_index_overwrite ==
                                          local_dof_indices_3[o])
                                        {
                                          // neighbor_cell contain the same dof
                                          // check if this cell_cut is cut if
                                          // it's not cut this dof must not
                                          // be overwritten
                                          bool cell_is_cut;
                                          std::tie(cell_is_cut,
                                                   std::ignore,
                                                   std::ignore) =
                                            cut_cells_map[neighbor_cell];
                                          bool cell_is_overconstrained;
                                          std::tie(cell_is_overconstrained,
                                                   std::ignore,
                                                   std::ignore) =
                                            overconstrained_fluid_cell_map
                                              [neighbor_cell];
                                          if (cell_is_cut == false &&
                                              cell_is_overconstrained == false)
                                            {
                                              dummy_dof = false;
                                              break;
                                            }
                                        }
                                    }
                                  if (dummy_dof == false)
                                    break;
                                }

                              if (dummy_dof)
                                {
                                  // The DOF is dummy
                                  this->system_matrix.add(
                                    global_index_overwrite,
                                    global_index_overwrite,
                                    sum_line);
                                  auto &system_rhs = this->system_rhs;
                                  system_rhs(global_index_overwrite) = 0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  this->system_rhs.compress(VectorOperation::insert);
  this->system_matrix.compress(VectorOperation::add);
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::setup_assemblers()
{
  this->assemblers.clear();
  assemblers_inside_ib.clear();

  if (this->check_existance_of_bc(
        BoundaryConditions::BoundaryType::function_weak))
    {
      this->assemblers.push_back(
        std::make_shared<WeakDirichletBoundaryCondition<dim>>(
          this->simulation_control,
          this->simulation_parameters.boundary_conditions));
    }
  if (this->check_existance_of_bc(BoundaryConditions::BoundaryType::pressure))
    {
      this->assemblers.push_back(
        std::make_shared<PressureBoundaryCondition<dim>>(
          this->simulation_control,
          this->simulation_parameters.boundary_conditions));
    }
  if (this->check_existance_of_bc(BoundaryConditions::BoundaryType::outlet))
    {
      this->assemblers.push_back(std::make_shared<OutletBoundaryCondition<dim>>(
        this->simulation_control,
        this->simulation_parameters.boundary_conditions));
    }
  if (this->check_existance_of_bc(
        BoundaryConditions::BoundaryType::partial_slip))
    {
      this->assemblers.push_back(
        std::make_shared<PartialSlipDirichletBoundaryCondition<dim>>(
          this->simulation_control,
          this->simulation_parameters.boundary_conditions));
    }
  if (this->simulation_parameters.multiphysics.VOF)
    {
      // Time-stepping schemes
      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesVOFAssemblerBDF<dim>>(
              this->simulation_control));
        }

      // Core assemblers
      if (this->simulation_parameters.physical_properties_manager
            .is_non_newtonian())
        {
          // Core assembler with Non newtonian viscosity
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesVOFAssemblerNonNewtonianCore<dim>>(
              this->simulation_control, this->simulation_parameters));
        }
      else
        {
          // Core assembler
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesVOFAssemblerCore<dim>>(
              this->simulation_control, this->simulation_parameters));
        }
    }
  else
    {
      // Time-stepping schemes
      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerBDF<dim>>(
              this->simulation_control));
        }
      else if (is_sdirk(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerSDIRK<dim>>(
              this->simulation_control));
        }

      // Velocity sources term
      if (this->simulation_parameters.velocity_sources.type ==
          Parameters::VelocitySource::VelocitySourceType::srf)
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerSRF<dim>>(
              this->simulation_parameters.velocity_sources));
        }

      // Core assemblers
      if (this->simulation_parameters.physical_properties_manager
            .is_non_newtonian())
        {
          // Core assembler with Non newtonian viscosity
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerNonNewtonianCore<dim>>(
              this->simulation_control));
        }
      else
        {
          // Core assembler
          this->assemblers.push_back(
            std::make_shared<PSPGSUPGNavierStokesAssemblerCore<dim>>(
              this->simulation_control));
        }
    }

  assemblers_inside_ib.push_back(
    std::make_shared<LaplaceAssembly<dim>>(this->simulation_control));
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim> &                        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &                copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();

  if (!cell->is_locally_owned())
    return;

  // We check if the cell is cut by the IB and if it is overconstrained. If the
  // cell is cut or overconstrained, we adjust the copy data so that the
  // assemblers don't execute work on this cell.
  bool cell_is_cut                                = false;
  std::tie(cell_is_cut, std::ignore, std::ignore) = cut_cells_map[cell];
  bool cell_is_overconstrained;
  std::tie(cell_is_overconstrained, std::ignore, std::ignore) =
    overconstrained_fluid_cell_map[cell];

  copy_data.cell_is_cut = cell_is_cut || cell_is_overconstrained;

  if (copy_data.cell_is_cut)
    return;
  scratch_data.reinit(
    cell,
    this->evaluation_point,
    this->previous_solutions,
    this->solution_stages,
    this->forcing_function,
    this->flow_control.get_beta(),
    this->simulation_parameters.stabilization.pressure_scaling_factor);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_vof);

      scratch_data.reinit_vof(
        phase_cell,
        *this->multiphysics->get_solution(PhysicsID::VOF),
        *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
        *this->multiphysics->get_previous_solutions(PhysicsID::VOF),
        std::vector<TrilinosWrappers::MPI::Vector>());
    }

  scratch_data.calculate_physical_properties();
  copy_data.reset();

  // check if we assemble the NS equation inside the particle or the Laplacian
  // of the variables
  bool cell_is_inside;
  std::tie(cell_is_inside, std::ignore) = cells_inside_map[cell];
  if (cell_is_inside && this->simulation_parameters.particlesParameters
                            ->assemble_navier_stokes_inside == false)
    {
      for (auto &assembler : this->assemblers_inside_ib)
        {
          assembler->assemble_matrix(scratch_data, copy_data);
        }
    }
  else
    {
      for (auto &assembler : this->assemblers)
        {
          assembler->assemble_matrix(scratch_data, copy_data);
        }
    }


  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::copy_local_matrix_to_global_matrix(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local || copy_data.cell_is_cut)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_matrix,
                                              copy_data.local_dof_indices,
                                              this->system_matrix);
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim> &                        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &                copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();

  if (!cell->is_locally_owned())
    return;

  // We check if the cell is cut by the IB and if it is overconstrained. If the
  // cell is cut or overconstrained, we adjust the copy data so that the
  // assemblers don't execute work on this cell.
  bool cell_is_cut                                = false;
  std::tie(cell_is_cut, std::ignore, std::ignore) = cut_cells_map[cell];

  bool cell_is_overconstrained;
  std::tie(cell_is_overconstrained, std::ignore, std::ignore) =
    overconstrained_fluid_cell_map[cell];

  copy_data.cell_is_cut = cell_is_cut || cell_is_overconstrained;

  if (copy_data.cell_is_cut)
    return;
  scratch_data.reinit(
    cell,
    this->evaluation_point,
    this->previous_solutions,
    this->solution_stages,
    this->forcing_function,
    this->flow_control.get_beta(),
    this->simulation_parameters.stabilization.pressure_scaling_factor);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_vof);

      scratch_data.reinit_vof(
        phase_cell,
        *this->multiphysics->get_solution(PhysicsID::VOF),
        *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
        *this->multiphysics->get_previous_solutions(PhysicsID::VOF),
        std::vector<TrilinosWrappers::MPI::Vector>());
    }

  scratch_data.calculate_physical_properties();
  copy_data.reset();

  // check if we assemble the NS equation inside the particle or the Laplacian
  // of the variables
  bool cell_is_inside;
  std::tie(cell_is_inside, std::ignore) = cells_inside_map[cell];
  if (cell_is_inside && this->simulation_parameters.particlesParameters
                            ->assemble_navier_stokes_inside == false)
    {
      for (auto &assembler : this->assemblers_inside_ib)
        {
          assembler->assemble_rhs(scratch_data, copy_data);
        }
    }
  else
    {
      for (auto &assembler : this->assemblers)
        {
          assembler->assemble_rhs(scratch_data, copy_data);
        }
    }

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::copy_local_rhs_to_global_rhs(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local || copy_data.cell_is_cut)
    return;


  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_rhs,
                                              copy_data.local_dof_indices,
                                              this->system_rhs);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::write_checkpoint()
{
  this->GLSNavierStokesSolver<dim>::write_checkpoint();


  // Write a table with all the relevant properties of the particle in a table.
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      std::string prefix =
        this->simulation_parameters.simulation_control.output_folder +
        this->simulation_parameters.restart_parameters.filename;

      TableHandler particles_information_table;
      std::string  filename =
        this->simulation_parameters.simulation_control.output_folder + prefix +
        ".ib_particles";
      std::ofstream output(filename.c_str());
      this->simulation_control->save(prefix);
      ib_particles_pvdhandler.save(prefix + ".ib_particles");

      this->pvdhandler.save(prefix);
      for (unsigned int i_particle = 0; i_particle < particles.size();
           ++i_particle)
        {
          for (unsigned int j = 0;
               j < particles[i_particle].previous_positions.size();
               ++j)
            {
              particles_information_table.add_value("ID", i_particle);
              particles_information_table.add_value(
                "p_x", particles[i_particle].previous_positions[j][0]);
              particles_information_table.set_precision("p_x", 16);
              particles_information_table.add_value(
                "p_y", particles[i_particle].previous_positions[j][1]);
              particles_information_table.set_precision("p_y", 16);
              if (dim == 3)
                {
                  particles_information_table.add_value(
                    "p_z", particles[i_particle].previous_positions[j][2]);
                  particles_information_table.set_precision("p_z", 16);
                }

              particles_information_table.add_value(
                "v_x", particles[i_particle].previous_velocity[j][0]);
              particles_information_table.set_precision("v_x", 16);
              particles_information_table.add_value(
                "v_y", particles[i_particle].previous_velocity[j][1]);
              particles_information_table.set_precision("v_y", 16);
              if (dim == 3)
                {
                  particles_information_table.add_value(
                    "v_z", particles[i_particle].previous_velocity[j][2]);
                  particles_information_table.set_precision("v_z", 16);
                }

              particles_information_table.add_value(
                "f_x", particles[i_particle].previous_fluid_forces[0]);
              particles_information_table.set_precision("f_x", 16);
              particles_information_table.add_value(
                "f_y", particles[i_particle].previous_fluid_forces[1]);
              particles_information_table.set_precision("f_y", 16);

              if (dim == 3)
                {
                  particles_information_table.add_value(
                    "f_z", particles[i_particle].previous_fluid_forces[2]);
                  particles_information_table.set_precision("f_z", 16);
                }

              particles_information_table.add_value(
                "f_xv", particles[i_particle].previous_fluid_viscous_forces[0]);
              particles_information_table.set_precision("f_xv", 16);
              particles_information_table.add_value(
                "f_yv", particles[i_particle].previous_fluid_viscous_forces[1]);
              particles_information_table.set_precision("f_yv", 16);

              if (dim == 3)
                {
                  particles_information_table.add_value(
                    "f_zv",
                    particles[i_particle].previous_fluid_viscous_forces[2]);
                  particles_information_table.set_precision("f_zv", 16);
                }

              particles_information_table.add_value(
                "f_xp",
                particles[i_particle].previous_fluid_pressure_forces[0]);
              particles_information_table.set_precision("f_xp", 16);
              particles_information_table.add_value(
                "f_yp",
                particles[i_particle].previous_fluid_pressure_forces[1]);
              particles_information_table.set_precision("f_yp", 16);

              if (dim == 3)
                {
                  particles_information_table.add_value(
                    "f_zp",
                    particles[i_particle].previous_fluid_pressure_forces[2]);
                  particles_information_table.set_precision("f_zp", 16);
                }

              if (dim == 3)
                {
                  particles_information_table.add_value(
                    "omega_x", particles[i_particle].previous_omega[j][0]);
                  particles_information_table.set_precision("omega_x", 16);
                  particles_information_table.add_value(
                    "omega_y", particles[i_particle].previous_omega[j][1]);
                  particles_information_table.set_precision("omega_y", 16);
                }
              particles_information_table.add_value(
                "omega_z", particles[i_particle].previous_omega[j][2]);
              particles_information_table.set_precision("omega_z", 16);
              if (dim == 3)
                {
                  particles_information_table.add_value(
                    "T_x", particles[i_particle].previous_fluid_torque[0]);
                  particles_information_table.set_precision("T_x", 16);
                  particles_information_table.add_value(
                    "T_y", particles[i_particle].previous_fluid_torque[1]);
                  particles_information_table.set_precision("T_y", 16);
                }
              particles_information_table.add_value(
                "T_z", particles[i_particle].previous_fluid_torque[2]);
              particles_information_table.set_precision("T_z", 16);
            }
        }
      // Write the table in the checkpoint file.
      particles_information_table.write_text(output);
    }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::read_checkpoint()
{
  this->GLSNavierStokesSolver<dim>::read_checkpoint();

  TimerOutput::Scope t(this->computing_timer,
                       "Reset Sharp-Edge particle information");

  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;

  std::string filename =
    this->simulation_parameters.simulation_control.output_folder + prefix +
    ".ib_particles";

  ib_particles_pvdhandler.read(prefix + ".ib_particles");
  // refill the table from checkpoint
  for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
    {
      std::string filename_table =
        this->simulation_parameters.simulation_control.output_folder +
        this->simulation_parameters.particlesParameters->ib_force_output_file +
        "." + Utilities::int_to_string(p_i, 2) + ".dat";
      fill_table_from_file(table_p[p_i], filename_table);
    }

  // Read the data of each particle and put the relevant information in a
  // vector.
  std::map<std::string, std::vector<double>> restart_data;
  fill_vectors_from_file(restart_data, filename);

  // Implement the data  in the particles.
  if (dim == 2)
    {
      unsigned int row = 0;
      for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
        {
          unsigned int j = 0;
          while (restart_data["ID"][row] == p_i and
                 row < restart_data["ID"].size())
            {
              if (j == 0)
                {
                  particles[p_i].position[0]     = restart_data["p_x"][row];
                  particles[p_i].position[1]     = restart_data["p_y"][row];
                  particles[p_i].velocity[0]     = restart_data["v_x"][row];
                  particles[p_i].velocity[1]     = restart_data["v_y"][row];
                  particles[p_i].fluid_forces[0] = restart_data["f_x"][row];
                  particles[p_i].fluid_forces[1] = restart_data["f_y"][row];
                  particles[p_i].fluid_viscous_forces[0] =
                    restart_data["f_xv"][row];
                  particles[p_i].fluid_viscous_forces[1] =
                    restart_data["f_yv"][row];
                  particles[p_i].fluid_pressure_forces[0] =
                    restart_data["f_xp"][row];
                  particles[p_i].fluid_pressure_forces[1] =
                    restart_data["f_yp"][row];
                  particles[p_i].omega[2]        = restart_data["omega_z"][row];
                  particles[p_i].fluid_torque[2] = restart_data["T_z"][row];

                  // fill previous time step
                  particles[p_i].previous_positions[j][0] =
                    restart_data["p_x"][row];
                  particles[p_i].previous_positions[j][1] =
                    restart_data["p_y"][row];
                  particles[p_i].previous_velocity[j][0] =
                    restart_data["v_x"][row];
                  particles[p_i].previous_velocity[j][1] =
                    restart_data["v_y"][row];
                  particles[p_i].previous_omega[j][2] =
                    restart_data["omega_z"][row];

                  particles[p_i].set_position(particles[p_i].position);
                  particles[p_i].set_orientation(particles[p_i].orientation);
                }
              else
                {
                  // fill previous time step
                  particles[p_i].previous_positions[j][0] =
                    restart_data["p_x"][row];
                  particles[p_i].previous_positions[j][1] =
                    restart_data["p_y"][row];
                  particles[p_i].previous_velocity[j][0] =
                    restart_data["v_x"][row];
                  particles[p_i].previous_velocity[j][1] =
                    restart_data["v_y"][row];
                  particles[p_i].previous_omega[j][2] =
                    restart_data["omega_z"][row];
                }
              row += 1;
              j += 1;
            }
        }
    }
  if (dim == 3)
    {
      unsigned int row = 0;
      for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
        {
          unsigned int j = 0;
          while (restart_data["ID"][row] == p_i and
                 row < restart_data["ID"].size())
            {
              if (j == 0)
                {
                  particles[p_i].position[0]     = restart_data["p_x"][row];
                  particles[p_i].position[1]     = restart_data["p_y"][row];
                  particles[p_i].position[2]     = restart_data["p_z"][row];
                  particles[p_i].velocity[0]     = restart_data["v_x"][row];
                  particles[p_i].velocity[1]     = restart_data["v_y"][row];
                  particles[p_i].velocity[2]     = restart_data["v_z"][row];
                  particles[p_i].fluid_forces[0] = restart_data["f_x"][row];
                  particles[p_i].fluid_forces[1] = restart_data["f_y"][row];
                  particles[p_i].fluid_forces[2] = restart_data["f_z"][row];
                  particles[p_i].fluid_viscous_forces[0] =
                    restart_data["f_xv"][row];
                  particles[p_i].fluid_viscous_forces[1] =
                    restart_data["f_yv"][row];
                  particles[p_i].fluid_viscous_forces[2] =
                    restart_data["f_zv"][row];
                  particles[p_i].fluid_pressure_forces[0] =
                    restart_data["f_xp"][row];
                  particles[p_i].fluid_pressure_forces[1] =
                    restart_data["f_yp"][row];
                  particles[p_i].fluid_pressure_forces[2] =
                    restart_data["f_zp"][row];
                  particles[p_i].omega[0]        = restart_data["omega_x"][row];
                  particles[p_i].omega[1]        = restart_data["omega_y"][row];
                  particles[p_i].omega[2]        = restart_data["omega_z"][row];
                  particles[p_i].fluid_torque[0] = restart_data["T_x"][row];
                  particles[p_i].fluid_torque[1] = restart_data["T_y"][row];
                  particles[p_i].fluid_torque[2] = restart_data["T_z"][row];

                  // fill previous time step
                  particles[p_i].previous_positions[j][0] =
                    restart_data["p_x"][row];
                  particles[p_i].previous_positions[j][1] =
                    restart_data["p_y"][row];
                  particles[p_i].previous_positions[j][2] =
                    restart_data["p_z"][row];
                  particles[p_i].previous_velocity[j][0] =
                    restart_data["v_x"][row];
                  particles[p_i].previous_velocity[j][1] =
                    restart_data["v_y"][row];
                  particles[p_i].previous_velocity[j][2] =
                    restart_data["v_z"][row];
                  particles[p_i].previous_omega[j][0] =
                    restart_data["omega_x"][row];
                  particles[p_i].previous_omega[j][1] =
                    restart_data["omega_y"][row];
                  particles[p_i].previous_omega[j][2] =
                    restart_data["omega_z"][row];

                  particles[p_i].set_position(particles[p_i].position);
                  particles[p_i].set_orientation(particles[p_i].orientation);
                }
              else
                {
                  // fill previous time step
                  particles[p_i].previous_positions[j][0] =
                    restart_data["p_x"][row];
                  particles[p_i].previous_positions[j][1] =
                    restart_data["p_y"][row];
                  particles[p_i].previous_positions[j][2] =
                    restart_data["p_z"][row];
                  particles[p_i].previous_velocity[j][0] =
                    restart_data["v_x"][row];
                  particles[p_i].previous_velocity[j][1] =
                    restart_data["v_y"][row];
                  particles[p_i].previous_velocity[j][2] =
                    restart_data["v_z"][row];
                  particles[p_i].previous_omega[j][0] =
                    restart_data["omega_x"][row];
                  particles[p_i].previous_omega[j][1] =
                    restart_data["omega_y"][row];
                  particles[p_i].previous_omega[j][2] =
                    restart_data["omega_z"][row];
                }
              row += 1;
              j += 1;
            }
        }
    }
  // Check if particles are spheres
  check_whether_all_particles_are_sphere();

  // Create the list of contact candidates
  ib_dem.update_contact_candidates();
  // Finish the time step of the particles.
  for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
    {
      particles[p_i].velocity_iter        = particles[p_i].velocity;
      particles[p_i].impulsion_iter       = particles[p_i].impulsion;
      particles[p_i].omega_iter           = particles[p_i].omega;
      particles[p_i].omega_impulsion_iter = particles[p_i].omega_impulsion;
      particles[p_i].residual_velocity    = DBL_MAX;
      particles[p_i].residual_omega       = DBL_MAX;
    }
  // Finish the time step of the particle.
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::load_particles_from_file()
{
  using numbers::PI;
  TimerOutput::Scope t(this->computing_timer,
                       "Reset Sharp-Edge particle information");
  this->pcout << "Loading particles from a file" << std::endl;
  std::string filename =
    this->simulation_parameters.particlesParameters->particles_file;

  // Read the data of each particle and put the relevant information in a
  // vector.
  std::map<std::string, std::vector<double>> restart_data;
  fill_vectors_from_file(restart_data, filename);
  particles.resize(restart_data["type"].size());

  this->pcout << "Particles found: " << particles.size() << std::endl;
  // Implement the data  in the particles.
  if (dim == 2)
    {
      // Loop over each line of the file and filling the particles properties.
      for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
        {
          particles[p_i].initialize_all();

          particles[p_i].particle_id = p_i;

          particles[p_i].position[0]    = restart_data["p_x"][p_i];
          particles[p_i].position[1]    = restart_data["p_y"][p_i];
          particles[p_i].velocity[0]    = restart_data["v_x"][p_i];
          particles[p_i].velocity[1]    = restart_data["v_y"][p_i];
          particles[p_i].orientation[2] = restart_data["orientation_z"][p_i];
          // Initialize shape of particles according to the type parameter in
          // the file
          if (restart_data["type"][p_i] == Shape<dim>::ShapeType::sphere)
            {
              std::vector<double> shape_argument(1);
              shape_argument[0] = restart_data["shape_argument_0"][p_i];
              particles[p_i].initialize_shape("sphere", shape_argument);
            }
          else if (restart_data["type"][p_i] ==
                   Shape<dim>::ShapeType::hyper_rectangle)
            {
              std::vector<double> shape_argument(2);
              shape_argument[0] = restart_data["shape_argument_0"][p_i];
              shape_argument[1] = restart_data["shape_argument_1"][p_i];
              particles[p_i].initialize_shape("hyper rectangle",
                                              shape_argument);
            }
          else if (restart_data["type"][p_i] ==
                   Shape<dim>::ShapeType::ellipsoid)
            {
              std::vector<double> shape_argument(2);
              shape_argument[0] = restart_data["shape_argument_0"][p_i];
              shape_argument[1] = restart_data["shape_argument_1"][p_i];
              particles[p_i].initialize_shape("ellipsoid", shape_argument);
            }


          particles[p_i].radius = particles[p_i].shape->effective_radius;
          particles[p_i].mass   = PI * particles[p_i].radius *
                                particles[p_i].radius *
                                restart_data["density"][p_i];

          particles[p_i].omega[2]      = restart_data["omega_z"][p_i];
          particles[p_i].inertia[0][0] = restart_data["inertia"][p_i];
          particles[p_i].inertia[1][1] = restart_data["inertia"][p_i];
          particles[p_i].inertia[2][2] = restart_data["inertia"][p_i];

          particles[p_i].pressure_location[0] = restart_data["pressure_x"][p_i];
          particles[p_i].pressure_location[1] = restart_data["pressure_y"][p_i];

          particles[p_i].youngs_modulus = restart_data["youngs_modulus"][p_i];
          particles[p_i].restitution_coefficient =
            restart_data["restitution_coefficient"][p_i];
          particles[p_i].friction_coefficient =
            restart_data["friction_coefficient"][p_i];
          particles[p_i].poisson_ratio = restart_data["poisson_ratio"][p_i];
          particles[p_i].rolling_friction_coefficient =
            restart_data["rolling_friction_coefficient"][p_i];
          particles[p_i].initialize_previous_solution();
          if (restart_data["integrate_motion"][p_i] == 0.0)
            {
              particles[p_i].integrate_motion = false;
            }
          else
            {
              particles[p_i].integrate_motion = true;
            }
        }
    }
  if (dim == 3)
    {
      // Loop over each line of the file and filling the particles properties.
      for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
        {
          particles[p_i].initialize_all();

          particles[p_i].particle_id    = p_i;
          particles[p_i].position[0]    = restart_data["p_x"][p_i];
          particles[p_i].position[1]    = restart_data["p_y"][p_i];
          particles[p_i].position[2]    = restart_data["p_z"][p_i];
          particles[p_i].velocity[0]    = restart_data["v_x"][p_i];
          particles[p_i].velocity[1]    = restart_data["v_y"][p_i];
          particles[p_i].velocity[2]    = restart_data["v_z"][p_i];
          particles[p_i].orientation[0] = restart_data["orientation_x"][p_i];
          particles[p_i].orientation[1] = restart_data["orientation_y"][p_i];
          particles[p_i].orientation[2] = restart_data["orientation_z"][p_i];

          if (restart_data["type"][p_i] == Shape<dim>::ShapeType::sphere)
            {
              std::vector<double> shape_argument(1);
              shape_argument[0] = restart_data["shape_argument_0"][p_i];
              particles[p_i].initialize_shape("sphere", shape_argument);
            }
          else if (restart_data["type"][p_i] ==
                   Shape<dim>::ShapeType::hyper_rectangle)
            {
              std::vector<double> shape_argument(3);
              shape_argument[0] = restart_data["shape_argument_0"][p_i];
              shape_argument[1] = restart_data["shape_argument_1"][p_i];
              shape_argument[2] = restart_data["shape_argument_2"][p_i];
              particles[p_i].initialize_shape("hyper rectangle",
                                              shape_argument);
            }
          else if (restart_data["type"][p_i] ==
                   Shape<dim>::ShapeType::ellipsoid)
            {
              std::vector<double> shape_argument(3);
              shape_argument[0] = restart_data["shape_argument_0"][p_i];
              shape_argument[1] = restart_data["shape_argument_1"][p_i];
              shape_argument[2] = restart_data["shape_argument_2"][p_i];
              particles[p_i].initialize_shape("ellipsoid", shape_argument);
            }
          else if (restart_data["type"][p_i] == Shape<dim>::ShapeType::torus)
            {
              std::vector<double> shape_argument(2);
              shape_argument[0] = restart_data["shape_argument_0"][p_i];
              shape_argument[1] = restart_data["shape_argument_1"][p_i];
              particles[p_i].initialize_shape("torus", shape_argument);
            }
          else if (restart_data["type"][p_i] == Shape<dim>::ShapeType::cone)
            {
              std::vector<double> shape_argument(2);
              shape_argument[0] = restart_data["shape_argument_0"][p_i];
              shape_argument[1] = restart_data["shape_argument_1"][p_i];
              particles[p_i].initialize_shape("cone", shape_argument);
            }
          else if (restart_data["type"][p_i] ==
                   Shape<dim>::ShapeType::cut_hollow_sphere)
            {
              std::vector<double> shape_argument(3);
              shape_argument[0] = restart_data["shape_argument_0"][p_i];
              shape_argument[1] = restart_data["shape_argument_1"][p_i];
              shape_argument[2] = restart_data["shape_argument_2"][p_i];
              particles[p_i].initialize_shape("cut hollow sphere",
                                              shape_argument);
            }
          else if (restart_data["type"][p_i] ==
                   Shape<dim>::ShapeType::death_star)
            {
              std::vector<double> shape_argument(3);
              shape_argument[0] = restart_data["shape_argument_0"][p_i];
              shape_argument[1] = restart_data["shape_argument_1"][p_i];
              shape_argument[2] = restart_data["shape_argument_2"][p_i];
              particles[p_i].initialize_shape("death star", shape_argument);
            }
          else if (restart_data["type"][p_i] ==
                   Shape<dim>::ShapeType::rbf_shape)
            {
              particles[p_i].shape =
                std::dynamic_pointer_cast<RBFShape<dim>>(
                  this->simulation_parameters.particlesParameters
                    ->particles[p_i]
                    .shape)
                  ->static_copy();
            }

          particles[p_i].radius = particles[p_i].shape->effective_radius;
          particles[p_i].mass   = 4.0 / 3.0 * PI * particles[p_i].radius *
                                particles[p_i].radius * particles[p_i].radius *
                                restart_data["density"][p_i];

          particles[p_i].omega[0] = restart_data["omega_x"][p_i];
          particles[p_i].omega[1] = restart_data["omega_y"][p_i];
          particles[p_i].omega[2] = restart_data["omega_z"][p_i];

          particles[p_i].inertia[0][0] = restart_data["inertia"][p_i];
          particles[p_i].inertia[1][1] = restart_data["inertia"][p_i];
          particles[p_i].inertia[2][2] = restart_data["inertia"][p_i];

          particles[p_i].pressure_location[0] = restart_data["pressure_x"][p_i];
          particles[p_i].pressure_location[1] = restart_data["pressure_y"][p_i];
          particles[p_i].pressure_location[2] = restart_data["pressure_z"][p_i];

          particles[p_i].youngs_modulus = restart_data["youngs_modulus"][p_i];
          particles[p_i].restitution_coefficient =
            restart_data["restitution_coefficient"][p_i];
          particles[p_i].friction_coefficient =
            restart_data["friction_coefficient"][p_i];
          particles[p_i].poisson_ratio = restart_data["poisson_ratio"][p_i];
          particles[p_i].rolling_friction_coefficient =
            restart_data["rolling_friction_coefficient"][p_i];
          particles[p_i].initialize_previous_solution();
          if (restart_data["integrate_motion"][p_i] == 0.0)
            {
              particles[p_i].integrate_motion = false;
            }
          else
            {
              particles[p_i].integrate_motion = true;
            }
        }
    }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::update_precalculations_for_ib()
{
  TimerOutput::Scope t(this->computing_timer,
                       "update_precalculations_shape_distance");
  for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
    {
      particles[p_i].update_precalculations(
        this->dof_handler,
        this->simulation_parameters.particlesParameters
          ->levels_not_precalculated);
    }
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::solve()
{
  MultithreadInfo::set_thread_limit(1);
  read_mesh_and_manifolds(
    *this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  define_particles();
  this->setup_dofs();

  if (this->simulation_parameters.restart_parameters.restart == false)
    {
      refinement_control(true);
      vertices_cell_mapping();

      if (all_spheres)
        optimized_generate_cut_cells_map();
      else
        generate_cut_cells_map();
    }

  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      this->forcing_function->set_time(
        this->simulation_control->get_current_time());

      if ((this->simulation_control->get_step_number() %
               this->simulation_parameters.mesh_adaptation.frequency !=
             0 ||
           this->simulation_parameters.mesh_adaptation.type ==
             Parameters::MeshAdaptation::Type::none ||
           this->simulation_control->is_at_start()) &&
          this->simulation_parameters.boundary_conditions.time_dependent)
        {
          this->update_boundary_conditions();
        }

      if (some_particles_are_coupled == false)
        integrate_particles();


      if (this->simulation_control->is_at_start())
        {
          vertices_cell_mapping();
          update_precalculations_for_ib();
          if (all_spheres)
            optimized_generate_cut_cells_map();
          else
            generate_cut_cells_map();

          ib_dem.update_particles_boundary_contact(this->particles,
                                                   this->dof_handler,
                                                   *this->face_quadrature,
                                                   *this->mapping);
          ib_dem.update_contact_candidates();
          this->iterate();
        }
      else
        {
          ib_done.clear();
          refinement_control(false);
          vertices_cell_mapping();

          if (all_spheres)
            optimized_generate_cut_cells_map();
          else
            generate_cut_cells_map();

          ib_dem.update_particles_boundary_contact(this->particles,
                                                   this->dof_handler,
                                                   *this->face_quadrature,
                                                   *this->mapping);
          if (this->simulation_control->get_step_number() == 0 ||
              this->simulation_control->get_step_number() %
                  this->simulation_parameters.particlesParameters
                    ->contact_search_frequency ==
                0)
            ib_dem.update_contact_candidates();

          // add initialization
          this->iterate();
        }

      this->postprocess_fd(false);

      if (this->simulation_parameters.particlesParameters->calculate_force_ib)
        force_on_ib();
      finish_time_step_particles();
      write_force_ib();
      this->finish_time_step();
    }

  if (this->simulation_parameters.particlesParameters->calculate_force_ib)
    this->finish_simulation();
}

// Pre-compile the 2D and 3D versopm solver to ensure that the library is
// valid before we actually compile the final solver
template class GLSSharpNavierStokesSolver<2>;
template class GLSSharpNavierStokesSolver<3>;
