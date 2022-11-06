#include <core/periodic_hills_grid.h>

#include <dem/read_mesh.h>


template <int dim, int spacedim>
void
read_mesh(const Parameters::Mesh &             mesh_parameters,
          const bool                           restart,
          const ConditionalOStream &           pcout,
          Triangulation<dim, spacedim> &       triangulation,
          double &                             triangulation_cell_diameter,
          const Parameters::Lagrangian::BCDEM &bc_params)
{
  pcout << "Reading triangulation " << std::endl;
  // GMSH input
  if (mesh_parameters.type == Parameters::Mesh::Type::gmsh)
    {
      GridIn<dim, spacedim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream input_file(mesh_parameters.file_name);
      grid_in.read_msh(input_file);
    }

  // Dealii grids
  else if (mesh_parameters.type == Parameters::Mesh::Type::dealii)
    {
      GridGenerator::generate_from_name_and_arguments(
        triangulation,
        mesh_parameters.grid_type,
        mesh_parameters.grid_arguments);
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

              // Add a cylindrical manifold on the final unrefined mesh
              const CylindricalManifold<3, spacedim> m1(0);
              triangulation.set_manifold(0, m1);
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
        }
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");

  triangulation_cell_diameter = 0.5 * GridTools::diameter(triangulation);

  if (bc_params.BC_type ==
      Parameters::Lagrangian::BCDEM::BoundaryType::periodic)
    {
      match_periodic_boundaries(triangulation, bc_params);
    }

  if (!restart)
    {
      if (mesh_parameters.refine_until_target_size)
        {
          double minimal_cell_size =
            GridTools::minimal_cell_diameter(triangulation);
          double       target_size = mesh_parameters.target_size;
          unsigned int number_refinement =
            floor(std::log(minimal_cell_size / target_size) / std::log(2));
          pcout << "Automatically refining grid until target size : "
                << target_size << std::endl;
          triangulation.refine_global(number_refinement);
          pcout << "Mesh was automatically refined : " << number_refinement
                << " times" << std::endl;
        }
      else
        {
          const int initial_refinement = mesh_parameters.initial_refinement;
          triangulation.refine_global(initial_refinement);
        }
    }

  pcout << std::endl << "Finished reading triangulation " << std::endl;
}

template <int dim, int spacedim>
void
match_periodic_boundaries(Triangulation<dim, spacedim> &       triangulation,
                          const Parameters::Lagrangian::BCDEM &bc_param)
{
  for (unsigned int i_bc = 0; i_bc < bc_param.DEM_BC_number; ++i_bc)
    {
      std::vector<GridTools::PeriodicFacePair<
        typename Triangulation<dim, spacedim>::cell_iterator>>
        periodicity_vector;

      GridTools::collect_periodic_faces(triangulation,
                                        bc_param.outlet_boundaries[i_bc],
                                        bc_param.periodic_boundaries[i_bc],
                                        bc_param.periodic_direction[i_bc],
                                        periodicity_vector);
      triangulation.add_periodicity(periodicity_vector);
    }
}


template void
read_mesh<1, 2>(const Parameters::Mesh &  mesh_parameters,
                const bool                restart,
                const ConditionalOStream &pcout,
                Triangulation<1, 2> &     triangulation,
                double &                  triangulation_cell_diameter,
                const Parameters::Lagrangian::BCDEM &bc_params);

template void
read_mesh<2, 2>(const Parameters::Mesh &  mesh_parameters,
                const bool                restart,
                const ConditionalOStream &pcout,
                Triangulation<2, 2> &     triangulation,
                double &                  triangulation_cell_diameter,
                const Parameters::Lagrangian::BCDEM &bc_params);

template void
read_mesh<2, 3>(const Parameters::Mesh &  mesh_parameters,
                const bool                restart,
                const ConditionalOStream &pcout,
                Triangulation<2, 3> &     triangulation,
                double &                  triangulation_cell_diameter,
                const Parameters::Lagrangian::BCDEM &bc_params);

template void
read_mesh<3, 3>(const Parameters::Mesh &  mesh_parameters,
                const bool                restart,
                const ConditionalOStream &pcout,
                Triangulation<3, 3> &     triangulation,
                double &                  triangulation_cell_diameter,
                const Parameters::Lagrangian::BCDEM &bc_params);
