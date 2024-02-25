
// Lethe includes
#include <core/boundary_conditions.h>
#include <core/grids.h>
#include <core/periodic_hills_grid.h>

// Deal.II includes
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

// Std
#include <fstream>
#include <iostream>


template <int dim, int spacedim>
void
attach_grid_to_triangulation(
  parallel::DistributedTriangulationBase<dim, spacedim> &triangulation,
  const Parameters::Mesh                                &mesh_parameters)

{
  // GMSH input
  if (mesh_parameters.type == Parameters::Mesh::Type::gmsh)
    {
      if (mesh_parameters.simplex)
        {
          auto        comm      = triangulation.get_communicator();
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
                GridTools::partition_triangulation(
                  Utilities::MPI::n_mpi_processes(comm), basetria);
              },
              comm,
              Utilities::MPI::n_mpi_processes(comm) /* group size */,
              dealii::Triangulation<dim, spacedim>::none);

          triangulation.create_triangulation(construction_data);
        }
      else
        {
          GridIn<dim, spacedim> grid_in;
          grid_in.attach_triangulation(triangulation);
          std::ifstream input_file(mesh_parameters.file_name);
          grid_in.read_msh(input_file);
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
            Utilities::MPI::n_mpi_processes(triangulation.get_communicator()),
            temporary_tri_triangulation);
          GridTools::partition_multigrid_levels(temporary_tri_triangulation);

          // extract relevant information from distributed triangulation
          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation(
              temporary_tri_triangulation,
              triangulation.get_communicator(),
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
            }
        }
    }

  // Colorized cylinder shell
  else if (mesh_parameters.type ==
           Parameters::Mesh::Type::colorized_cylinder_shell)
    {
      if (mesh_parameters.simplex)
        {
          throw std::runtime_error(
            "Unsupported mesh type - colorized cylinder shell with simplex is not supported.");
        }
      else if (dim != 3)
        {
          throw std::runtime_error(
            "Unsupported mesh type - colorized cylinder shell is only supported in 3d.");
        }

      // First generate the regular deal.II cylinder_shell using the
      // mesh_parameters arguments
      GridGenerator::generate_from_name_and_arguments(
        triangulation, "cylinder_shell", mesh_parameters.grid_arguments);

      // Now that we have our grid, colorize the boundary conditions.
      // Inner cylinder has boundary id 0
      // Outer cylinder has boundary id 1
      // Bottom (Z-) part of the cylinder has boundary id 2
      // Top (Z+) part of the cylinder has boundary id 3

      // Split the argument to extract the radius
      std::vector<std::string> split_arguments =
        Utilities::split_string_list(mesh_parameters.grid_arguments, ":");

      const double length = Utilities::string_to_double(split_arguments[0]);

      // Tolerance along z
      const double eps_z = 1e-6 * length;

      // Gather the inner radius from the faces instead of the argument, this is
      // more robust for some aspect ratios. First initialize the outer to 0 and
      // the inner to a large value
      double inner_radius = DBL_MAX;
      double outer_radius = 0.;

      // Loop over the cells once to acquire the min and max radius at the face
      // centers Otherwise, for some cell ratio, the center of the faces can be
      // at a radius which is significantly different from the one prescribed.
      for (const auto &cell : triangulation.active_cell_iterators())
        for (const unsigned int f : GeometryInfo<3>::face_indices())
          {
            if (!cell->face(f)->at_boundary())
              continue;

            const auto   face_center = cell->face(f)->center();
            const double z           = face_center[2];

            if ((std::fabs(z) > eps_z) &&
                (std::fabs(z - length) > eps_z)) // Not a zmin or zmax boundary
              {
                const double radius =
                  std::sqrt(face_center[0] * face_center[0] +
                            face_center[1] * face_center[1]);
                inner_radius = std::min(inner_radius, radius);
                outer_radius = std::max(outer_radius, radius);
              }
          }

      // Use the gathered radius to define the medium radial distance.
      double mid_radial_distance = 0.5 * (outer_radius - inner_radius);

      for (const auto &cell : triangulation.active_cell_iterators())
        for (const unsigned int f : GeometryInfo<3>::face_indices())
          {
            if (!cell->face(f)->at_boundary())
              continue;

            const auto face_center = cell->face(f)->center();

            const double radius = std::sqrt(face_center[0] * face_center[0] +
                                            face_center[1] * face_center[1]);

            const double z = face_center[2];

            if (std::fabs(z) < eps_z) // z = 0 set boundary 2
              {
                cell->face(f)->set_boundary_id(2);
              }
            else if (std::fabs(z - length) < eps_z) // z = length set boundary 3
              {
                cell->face(f)->set_boundary_id(3);
              }
            else if (std::fabs(radius - inner_radius) >
                     mid_radial_distance) // r =  outer_radius set boundary 1
              {
                cell->face(f)->set_boundary_id(1);
              }
            else if (std::fabs(radius - inner_radius) <
                     mid_radial_distance) // r =  inner_radius set boundary 0
              {
                cell->face(f)->set_boundary_id(0);
              }
            else
              throw std::runtime_error(
                "There was an error while setting up this mesh");
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
        }
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");
}


template <int dim, int spacedim>
void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<dim, spacedim>  &triangulation,
  const BoundaryConditions::BoundaryConditions<spacedim> &boundary_conditions)

{
  // Setup periodic boundary conditions
  for (unsigned int i_bc = 0; i_bc < boundary_conditions.size; ++i_bc)
    {
      if (boundary_conditions.type[i_bc] ==
          BoundaryConditions::BoundaryType::periodic)
        {
          auto triangulation_ptr = &triangulation;

          std::vector<GridTools::PeriodicFacePair<
            typename Triangulation<dim, spacedim>::cell_iterator>>
            periodicity_vector;
          GridTools::collect_periodic_faces(
            *dynamic_cast<Triangulation<dim, spacedim> *>(triangulation_ptr),
            boundary_conditions.id[i_bc],
            boundary_conditions.periodic_id[i_bc],
            boundary_conditions.periodic_direction[i_bc],
            periodicity_vector);
          triangulation.add_periodicity(periodicity_vector);
        }
    }
}

template <int dim, int spacedim>
void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<dim, spacedim>  &triangulation,
  const Parameters::Mesh                                 &mesh_parameters,
  const Parameters::Manifolds                            &manifolds_parameters,
  const bool                                             &restart,
  const BoundaryConditions::BoundaryConditions<spacedim> &boundary_conditions)
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
                  manifold_ids.find(face->boundary_id()) != manifold_ids.end())
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
          double       target_size = mesh_parameters.target_size;
          unsigned int number_refinement =
            floor(std::log(minimal_cell_size / target_size) / std::log(2));
          triangulation.refine_global(number_refinement);
        }
      else if (!restart)
        {
          const int initial_refinement = mesh_parameters.initial_refinement;
          triangulation.refine_global(initial_refinement);
        }
    }
}

template void
attach_grid_to_triangulation(
  parallel::DistributedTriangulationBase<2> &triangulation,
  const Parameters::Mesh                    &mesh_parameters);
template void
attach_grid_to_triangulation(
  parallel::DistributedTriangulationBase<3> &triangulation,
  const Parameters::Mesh                    &mesh_parameters);
template void
attach_grid_to_triangulation(
  parallel::DistributedTriangulationBase<2, 3> &triangulation,
  const Parameters::Mesh                       &mesh_parameters);


template void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<2, 2>    &triangulation,
  const BoundaryConditions::BoundaryConditions<2> &boundary_conditions);
template void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<2, 3>    &triangulation,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);
template void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<3, 3>    &triangulation,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);

template void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<2>       &triangulation,
  const Parameters::Mesh                          &mesh_parameters,
  const Parameters::Manifolds                     &manifolds_parameters,
  const bool                                      &restart,
  const BoundaryConditions::BoundaryConditions<2> &boundary_conditions);
template void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<3>       &triangulation,
  const Parameters::Mesh                          &mesh_parameters,
  const Parameters::Manifolds                     &manifolds_parameters,
  const bool                                      &restart,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);
template void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<2, 3>    &triangulation,
  const Parameters::Mesh                          &mesh_parameters,
  const Parameters::Manifolds                     &manifolds_parameters,
  const bool                                      &restart,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);
