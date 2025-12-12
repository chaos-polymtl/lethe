// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/boundary_conditions.h>
#include <core/grids.h>
#include <core/periodic_hills_grid.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <numbers>


template <int dim, int spacedim>
void
attach_grid_to_triangulation(Triangulation<dim, spacedim> &triangulation,
                             const Parameters::Mesh       &mesh_parameters)

{
  // GMSH input
  if (mesh_parameters.type == Parameters::Mesh::Type::gmsh)
    {
      if (mesh_parameters.simplex)
        {
          auto        comm      = triangulation.get_mpi_communicator();
          std::string file_name = mesh_parameters.file_name;

          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation_in_groups<dim, spacedim>(
              [file_name](dealii::Triangulation<dim, spacedim> &basetria) {
                GridIn<dim, spacedim> grid_in;
                grid_in.attach_triangulation(basetria);
                std::ifstream input_file(file_name);

                grid_in.read_msh(input_file);
              },
              [](dealii::Triangulation<dim, spacedim> &basetria,
                 const MPI_Comm                        comm,
                 const unsigned int /*group_size*/) {
#if defined(DEAL_II_WITH_METIS)
                GridTools::partition_triangulation(
                  Utilities::MPI::n_mpi_processes(comm),
                  basetria,
                  SparsityTools::Partitioner::metis);
#elif defined(DEAL_II_WITH_ZOLTAN)
                GridTools::partition_triangulation(
                  Utilities::MPI::n_mpi_processes(comm),
                  basetria,
                  SparsityTools::Partitioner::zoltan);
#else
                AssertThrow(
                  false,
                  ExcMessage(
                    "Parallel simulation with simplex meshes require that deal.II be compiled with either Metis or Zoltan"));
#endif
              },
              comm,
              Utilities::MPI::n_mpi_processes(comm) /* group size */,
              dealii::Triangulation<dim, spacedim>::none);

          triangulation.create_triangulation(construction_data);
          GridTools::scale(mesh_parameters.scale, triangulation);
        }
      else
        {
          GridIn<dim, spacedim> grid_in;
          grid_in.attach_triangulation(triangulation);
          std::ifstream input_file(mesh_parameters.file_name);
          grid_in.read_msh(input_file);
          GridTools::scale(mesh_parameters.scale, triangulation);
        }
    }
  // Dealii grids
  else if (mesh_parameters.type == Parameters::Mesh::Type::dealii)
    {
      if (mesh_parameters.simplex)
        {
          Triangulation<dim, spacedim> temporary_quad_triangulation;
          GridGenerator::generate_from_name_and_arguments(
            temporary_quad_triangulation,
            mesh_parameters.grid_type,
            mesh_parameters.grid_arguments);

          GridTools::scale(mesh_parameters.scale, temporary_quad_triangulation);

          // initial refinement
          const unsigned int initial_refinement =
            mesh_parameters.initial_refinement;
          temporary_quad_triangulation.refine_global(initial_refinement);

          // flatten the triangulation
          Triangulation<dim, spacedim> flat_temp_quad_triangulation;
          GridGenerator::flatten_triangulation(temporary_quad_triangulation,
                                               flat_temp_quad_triangulation);

          Triangulation<dim, spacedim> temporary_tri_triangulation(
            Triangulation<dim, spacedim>::limit_level_difference_at_vertices);
          GridGenerator::convert_hypercube_to_simplex_mesh(
            flat_temp_quad_triangulation, temporary_tri_triangulation);

          GridTools::partition_triangulation_zorder(
            Utilities::MPI::n_mpi_processes(
              triangulation.get_mpi_communicator()),
            temporary_tri_triangulation);
          GridTools::partition_multigrid_levels(temporary_tri_triangulation);

          // extract relevant information from distributed triangulation
          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation(
              temporary_tri_triangulation,
              triangulation.get_mpi_communicator(),
              TriangulationDescription::Settings::
                construct_multigrid_hierarchy);
          triangulation.create_triangulation(construction_data);
        }
      else
        {
          GridGenerator::generate_from_name_and_arguments(
            triangulation,
            mesh_parameters.grid_type,
            mesh_parameters.grid_arguments);

          GridTools::scale(mesh_parameters.scale, triangulation);

          if constexpr (dim == 3)
            {
              // Initial mesh translation
              GridTools::shift(mesh_parameters.translation, triangulation);

              // Initial mesh rotation
              GridTools::rotate(mesh_parameters.rotation_axis,
                                mesh_parameters.rotation_angle,
                                triangulation);
            }
        }
    }
  // Customizable cylinder mesh
  else if (mesh_parameters.type == Parameters::Mesh::Type::cylinder)
    {
      if (mesh_parameters.simplex)
        {
          throw std::runtime_error(
            "Unsupported mesh type - custom cylinder mesh with simplex is not supported. Use a dealii cylinder to use simplex mesh.");
        }
      else if (dim != 3)
        {
          throw std::runtime_error(
            "Unsupported mesh type - custom cylinder mesh is only supported in 3d.");
        }

      if constexpr (dim == 3)
        {
          // Separate arguments of the string
          std::vector<std::string> arguments;
          std::stringstream        s_stream(mesh_parameters.grid_arguments);
          while (s_stream.good())
            {
              std::string substr;
              getline(s_stream, substr, ':');
              arguments.push_back(substr);
            }

          // Arguments declaration
          unsigned int subdivisions;
          double       radius, half_height;
          if (arguments.size() != 3)
            {
              throw std::runtime_error(
                "Mandatory cylinder parameters are (x subdivisions: radius : half height)");
            }
          else
            {
              std::vector<double> arguments_double =
                dealii::Utilities::string_to_double(arguments);
              subdivisions = static_cast<int>(arguments_double[0]);
              radius       = arguments_double[1];
              half_height  = arguments_double[2];
            }

          if (mesh_parameters.grid_type == "classic")
            {
              // Create a subdivided cylinder from deal.ii
              GridGenerator::subdivided_cylinder(triangulation,
                                                 subdivisions,
                                                 radius,
                                                 half_height);

              GridTools::scale(mesh_parameters.scale, triangulation);
            }
          else
            {
              // Create a temporary 2d mesh
              Triangulation<2, spacedim - 1> temporary_triangulation;

              // Create a spherical manifold for 2d mesh
              Point<2>                                 center(0.0, 0.0);
              const SphericalManifold<2, spacedim - 1> m0(center);

              if (mesh_parameters.grid_type == "regularized" ||
                  mesh_parameters.grid_type == "squared")
                {
                  // Create a square mesh
                  double real_radius = radius * std::sin(M_PI_4);
                  GridGenerator::hyper_cube(temporary_triangulation,
                                            -real_radius,
                                            real_radius,
                                            true);

                  // Assign boundary 0 to perimeter as for cylinder
                  for (const auto &cell :
                       temporary_triangulation.active_cell_iterators())
                    {
                      if (cell->is_locally_owned())
                        {
                          // Looping through all the faces of the cell
                          for (const auto &face : cell->face_iterators())
                            {
                              // Check to see if the face is located at boundary
                              if (face->at_boundary())
                                {
                                  face->set_boundary_id(0);
                                }
                            }
                        }
                    }
                }
              else if (mesh_parameters.grid_type == "balanced")
                {
                  GridGenerator::hyper_ball_balanced(temporary_triangulation,
                                                     center,
                                                     radius);
                }
              else
                {
                  throw std::runtime_error(
                    "Unknown grid type. Choices are <classic|balanced|squared|regularized>.");
                }

              temporary_triangulation.reset_all_manifolds();
              temporary_triangulation.set_all_manifold_ids_on_boundary(0);
              temporary_triangulation.set_manifold(0, m0);

              if (mesh_parameters.grid_type == "regularized")
                {
                  // Pre-refinement to reduce mesh size at corners before
                  // regularization
                  temporary_triangulation.refine_global(2);
                  GridTools::regularize_corner_cells(temporary_triangulation);

                  // Flatten the triangulation
                  Triangulation<2, spacedim - 1> flat_temporary_triangulation;
                  flat_temporary_triangulation.copy_triangulation(
                    temporary_triangulation);
                  temporary_triangulation.clear();
                  GridGenerator::flatten_triangulation(
                    flat_temporary_triangulation, temporary_triangulation);
                }

              // Extrude the 2d temporary mesh to 3d cylinder
              GridGenerator::extrude_triangulation(temporary_triangulation,
                                                   subdivisions + 1,
                                                   2.0 * half_height,
                                                   triangulation,
                                                   true);

              // Rotate mesh in x-axis and set the (0,0,0) at the barycenter
              // to be comparable to dealii cylinder meshes
              Tensor<1, spacedim> axis_vector({0.0, 1.0, 0.0});
              GridTools::rotate(axis_vector, M_PI_2, triangulation);
              Tensor<1, spacedim> shift_vector({-half_height, 0.0, 0.0});
              GridTools::shift(shift_vector, triangulation);

              // Force the manifold id to be zero in the case of the balanced
              // cylinder
              if (mesh_parameters.grid_type == "balanced")
                triangulation.reset_manifold(1);

              // Add a cylindrical manifold on the final unrefined mesh
              const CylindricalManifold<3, spacedim> m1(0);
              triangulation.set_manifold(0, m1);

              GridTools::scale(mesh_parameters.scale, triangulation);
            }
        }
    }

  // Periodic Hills grid
  else if (mesh_parameters.type == Parameters::Mesh::Type::periodic_hills &&
           !mesh_parameters.simplex)
    {
      if (mesh_parameters.simplex)
        {
          throw std::runtime_error(
            "Unsupported mesh type - periodic hills mesh with simplex is not supported");
        }
      else
        {
          PeriodicHillsGrid<dim, spacedim> grid(mesh_parameters.grid_arguments);
          grid.make_grid(triangulation);

          GridTools::scale(mesh_parameters.scale, triangulation);
        }
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");
}


template <int dim, int spacedim>
void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<dim, spacedim> &triangulation,
  const BoundaryConditions::BoundaryConditions          &boundary_conditions)

{
  // Setup periodic boundary conditions
  for (auto const &[id, type] : boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          auto triangulation_ptr = &triangulation;

          std::vector<GridTools::PeriodicFacePair<
            typename Triangulation<dim, spacedim>::cell_iterator>>
            periodicity_vector;
          GridTools::collect_periodic_faces(
            *dynamic_cast<Triangulation<dim, spacedim> *>(triangulation_ptr),
            id,
            boundary_conditions.periodic_neighbor_id.at(id),
            boundary_conditions.periodic_direction.at(id),
            periodicity_vector);
          triangulation.add_periodicity(periodicity_vector);
        }
    }
}

template <int dim, int spacedim>
void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<dim, spacedim> &triangulation,
  const Parameters::Mesh                                &mesh_parameters,
  const Parameters::Manifolds                           &manifolds_parameters,
  const bool                                             restart,
  const BoundaryConditions::BoundaryConditions          &boundary_conditions)
{
  attach_grid_to_triangulation(triangulation, mesh_parameters);

  setup_periodic_boundary_conditions(triangulation, boundary_conditions);

  // If the mesh is a GMSH mesh, we need to manually set the manifold id
  // of the face to be that of the boundary
  // ID, otherwise, we will be unable to use external manifold.
  // This is done manually by looping through all faces and giving them a
  // manifold id if there is a manifold associated with this number.
  if (mesh_parameters.type == Parameters::Mesh::Type::gmsh)
    {
      // Gather all the manifold ids within a set
      std::set<int> manifold_ids;
      for (unsigned int i = 0; i < manifolds_parameters.size; ++i)
        manifold_ids.insert(manifolds_parameters.id[i]);

      // Reset all the manifolds manually and force them to zero
      triangulation.reset_all_manifolds();

      // If the parameter file forces the occurrence of manifold,
      // loop over the faces of the triangulation. If the face of the
      // triangulation has a boundary id which corresponds to a manifold id
      // identified within the parameter file, then fix the manifold id of this
      // face manually to be that of the boundary id. In the past, this was done
      // by default for every face, but since 2023-12 this throws (rightfully)
      // an error in deal.II
      if (manifolds_parameters.size > 0)
        {
          for (const auto &face : triangulation.active_face_iterators())
            {
              if (face->at_boundary() &&
                  manifold_ids.contains(face->boundary_id()))
                face->set_all_manifold_ids(face->boundary_id());
            }
        }
    }

  // Finally attach the manifolds to the triangulation
  // Right now this should only occur for GMSH mesh, but the function is generic
  // enough.
  attach_manifolds_to_triangulation(triangulation, manifolds_parameters);

  if (mesh_parameters.simplex)
    {
      // Refinement isn't possible yet
    }
  else
    {
      if (mesh_parameters.refine_until_target_size)
        {
          double minimal_cell_size =
            GridTools::minimal_cell_diameter(triangulation);
          double       target_size       = mesh_parameters.target_size;
          unsigned int number_refinement = static_cast<unsigned int>(floor(
            std::log(minimal_cell_size / target_size) / std::numbers::ln2));
          triangulation.refine_global(number_refinement);
        }
      else if (!restart)
        {
          const int initial_refinement = mesh_parameters.initial_refinement;
          triangulation.refine_global(initial_refinement);
          refine_triangulation_at_boundaries(
            mesh_parameters.boundaries_to_refine,
            mesh_parameters.initial_refinement_at_boundaries,
            triangulation);
        }
    }
}

template <int dim, int spacedim>
void
read_mesh_and_manifolds_for_stator_and_rotor(
  parallel::DistributedTriangulationBase<dim, spacedim> &triangulation,
  const Parameters::Mesh                                &mesh_parameters,
  const Parameters::Manifolds                           &manifolds_parameters,
  const bool                                             restart,
  const BoundaryConditions::BoundaryConditions          &boundary_conditions,
  const Parameters::Mortar<dim>                         &mortar_parameters)
{
  // First check if stator and rotor meshes are of the same type
  AssertThrow(
    mesh_parameters.type == mortar_parameters.rotor_mesh->type,
    ExcMessage(
      "The mesh types for the rotor and stator geometries must be the same."));

  // Since the rotor and stator meshes are read separately, a dummy
  // triangulation is created for each part of the domain and then merged

  // Stator triangulation
  Triangulation<dim> stator_temp_tria;
  attach_grid_to_triangulation(stator_temp_tria, mesh_parameters);

  // Rotor triangulation
  Triangulation<dim> rotor_temp_tria;
  attach_grid_to_triangulation(rotor_temp_tria, *mortar_parameters.rotor_mesh);

  if (mesh_parameters.type == Parameters::Mesh::Type::dealii)
    {
      // Get stator manifold ids without flat id
      unsigned int stator_ids_no_flat = 0;
      for (const auto &id : stator_temp_tria.get_manifold_ids())
        {
          if (id != numbers::flat_manifold_id)
            stator_ids_no_flat++;
        }

      // Get rotor manifold ids without flat id
      unsigned int rotor_ids_no_flat = 0;
      for (const auto &id : rotor_temp_tria.get_manifold_ids())
        {
          if (id != numbers::flat_manifold_id)
            rotor_ids_no_flat++;
        }

      // Shift rotor boundary IDs #
      unsigned int n_boundary_ids_stator =
        stator_temp_tria.get_boundary_ids().size();
      for (const auto &face : rotor_temp_tria.active_face_iterators())
        if (face->at_boundary())
          face->set_boundary_id(face->boundary_id() + n_boundary_ids_stator);

      // Keep track of modified IDs in faces and lines
      std::vector<bool> changed_manifold_ids_faces(
        rotor_temp_tria.n_active_faces(), false);
      std::vector<bool> changed_manifold_ids_lines(
        rotor_temp_tria.n_active_lines(), false);

      // Shift rotor manifold IDs #
      for (const auto &cell : rotor_temp_tria.active_cell_iterators())
        {
          if (cell->manifold_id() != numbers::flat_manifold_id)
            cell->set_manifold_id(cell->manifold_id() + stator_ids_no_flat);
          for (const auto &face : cell->face_iterators())
            {
              if (face->manifold_id() != numbers::flat_manifold_id &&
                  !changed_manifold_ids_faces[face->index()])
                {
                  face->set_manifold_id(face->manifold_id() +
                                        stator_ids_no_flat);
                  changed_manifold_ids_faces[face->index()] = true;
                }
              if constexpr (dim == 3)
                {
                  for (const auto line_index : cell->line_indices())
                    {
                      if (cell->line(line_index)->manifold_id() !=
                            numbers::flat_manifold_id &&
                          !changed_manifold_ids_lines[line_index])
                        cell->line(line_index)
                          ->set_manifold_id(
                            cell->line(line_index)->manifold_id() +
                            stator_ids_no_flat);
                      changed_manifold_ids_lines[line_index] = true;
                    }
                }
            }
        }

      // Store rotor manifolds in shifted IDs #
      for (unsigned int m = 0; m < rotor_temp_tria.get_manifold_ids().size();
           m++)
        {
          unsigned int temp = rotor_temp_tria.get_manifold_ids()[m];
          if (temp != numbers::flat_manifold_id)
            rotor_temp_tria.set_manifold(m + stator_ids_no_flat,
                                         rotor_temp_tria.get_manifold(m));
        }

      // Merge triangulations
      GridGenerator::merge_triangulations(
        stator_temp_tria, rotor_temp_tria, triangulation, 0.0, true, true);

      // Attach manifolds to merged triangulation
      unsigned int n = 0;
      for (unsigned int i = 0; i < stator_ids_no_flat; i++)
        {
          triangulation.set_manifold(n, stator_temp_tria.get_manifold(i));
          n++;
        }
      for (unsigned int j = stator_ids_no_flat;
           j < stator_ids_no_flat + rotor_ids_no_flat;
           j++)
        {
          triangulation.set_manifold(n, rotor_temp_tria.get_manifold(j));
          n++;
        }
    }
  else if (mesh_parameters.type == Parameters::Mesh::Type::gmsh)
    {
      // Merge triangulations
      GridGenerator::merge_triangulations(
        stator_temp_tria, rotor_temp_tria, triangulation, 0.0, true, true);

      // Gather all the manifold ids within a set
      std::set<int> manifold_ids;
      for (unsigned int i = 0; i < manifolds_parameters.size; ++i)
        manifold_ids.insert(manifolds_parameters.id[i]);

      // Reset all the manifolds manually and force them to zero
      triangulation.reset_all_manifolds();

      // Set manifolds
      if (manifolds_parameters.size > 0)
        {
          for (const auto &face : triangulation.active_face_iterators())
            {
              if (face->at_boundary() &&
                  manifold_ids.contains(face->boundary_id()))
                face->set_all_manifold_ids(face->boundary_id());
            }
        }

      // Attach manifolds to rotor-stator triangulation
      attach_manifolds_to_triangulation(triangulation, manifolds_parameters);
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");

  // Setup boundary conditions
  setup_periodic_boundary_conditions(triangulation, boundary_conditions);

  // Initial mesh refinement
  if (mesh_parameters.simplex)
    {
      // Refinement isn't possible yet
    }
  else
    {
      if (mesh_parameters.refine_until_target_size)
        {
          const double minimal_cell_size =
            GridTools::minimal_cell_diameter(triangulation);
          const double       target_size = mesh_parameters.target_size;
          const unsigned int number_refinement =
            static_cast<unsigned int>(floor(
              std::log(minimal_cell_size / target_size) / std::numbers::ln2));
          triangulation.refine_global(number_refinement);
        }
      else if (!restart)
        {
          const int initial_refinement = mesh_parameters.initial_refinement;
          triangulation.refine_global(initial_refinement);
          refine_triangulation_at_boundaries(
            mesh_parameters.boundaries_to_refine,
            mesh_parameters.initial_refinement_at_boundaries,
            triangulation);
        }
    }

  // Faces at rotor-stator interface
  unsigned int n_faces_rotor_interface  = 0;
  unsigned int n_faces_stator_interface = 0;

  // Check number of faces at the rotor-stator interface
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto face_no : cell->face_indices())
        if (cell->face(face_no)->at_boundary())
          {
            if (cell->face(face_no)->boundary_id() ==
                mortar_parameters.rotor_boundary_id)
              n_faces_rotor_interface++;
            if (cell->face(face_no)->boundary_id() ==
                mortar_parameters.stator_boundary_id)
              n_faces_stator_interface++;
          }

  // Total number of faces
  const unsigned int n_faces_rotor_interface_total =
    Utilities::MPI::sum(n_faces_rotor_interface,
                        triangulation.get_mpi_communicator());
  const unsigned int n_faces_stator_interface_total =
    Utilities::MPI::sum(n_faces_stator_interface,
                        triangulation.get_mpi_communicator());

  AssertThrow(
    n_faces_rotor_interface_total == n_faces_stator_interface_total,
    ExcMessage(
      "The number of faces at the rotor interface ID #" +
      std::to_string(mortar_parameters.rotor_boundary_id) + " (" +
      std::to_string(n_faces_rotor_interface_total) +
      ") is different from the number of faces at the stator interface ID #" +
      std::to_string(mortar_parameters.stator_boundary_id) + " (" +
      std::to_string(n_faces_rotor_interface_total) + ")."));
}

template void
attach_grid_to_triangulation(Triangulation<2>       &triangulation,
                             const Parameters::Mesh &mesh_parameters);
template void
attach_grid_to_triangulation(Triangulation<3>       &triangulation,
                             const Parameters::Mesh &mesh_parameters);
template void
attach_grid_to_triangulation(Triangulation<2, 3>    &triangulation,
                             const Parameters::Mesh &mesh_parameters);


template void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<2, 2> &triangulation,
  const BoundaryConditions::BoundaryConditions &boundary_conditions);
template void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<2, 3> &triangulation,
  const BoundaryConditions::BoundaryConditions &boundary_conditions);
template void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<3, 3> &triangulation,
  const BoundaryConditions::BoundaryConditions &boundary_conditions);

template void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<2>    &triangulation,
  const Parameters::Mesh                       &mesh_parameters,
  const Parameters::Manifolds                  &manifolds_parameters,
  const bool                                    restart,
  const BoundaryConditions::BoundaryConditions &boundary_conditions);
template void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<3>    &triangulation,
  const Parameters::Mesh                       &mesh_parameters,
  const Parameters::Manifolds                  &manifolds_parameters,
  const bool                                    restart,
  const BoundaryConditions::BoundaryConditions &boundary_conditions);
template void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<2, 3> &triangulation,
  const Parameters::Mesh                       &mesh_parameters,
  const Parameters::Manifolds                  &manifolds_parameters,
  const bool                                    restart,
  const BoundaryConditions::BoundaryConditions &boundary_conditions);

template void
read_mesh_and_manifolds_for_stator_and_rotor(
  parallel::DistributedTriangulationBase<2>    &triangulation,
  const Parameters::Mesh                       &mesh_parameters,
  const Parameters::Manifolds                  &manifolds_parameters,
  const bool                                    restart,
  const BoundaryConditions::BoundaryConditions &boundary_conditions,
  const Parameters::Mortar<2>                  &mortar_parameters);
template void
read_mesh_and_manifolds_for_stator_and_rotor(
  parallel::DistributedTriangulationBase<3>    &triangulation,
  const Parameters::Mesh                       &mesh_parameters,
  const Parameters::Manifolds                  &manifolds_parameters,
  const bool                                    restart,
  const BoundaryConditions::BoundaryConditions &boundary_conditions,
  const Parameters::Mortar<3>                  &mortar_parameters);
