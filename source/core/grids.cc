
// Lethe includes
#include <core/boundary_conditions.h>
#include <core/grids.h>
#include <core/periodic_hills_grid.h>

// Deal.II includes
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>


// Std
#include <fstream>
#include <iostream>


template <int dim, int spacedim>
void
attach_grid_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim, spacedim>>
                                                          triangulation,
  const Parameters::Mesh &                                mesh_parameters,
  const BoundaryConditions::BoundaryConditions<spacedim> &boundary_conditions)

{
  // GMSH input
  if (mesh_parameters.type == Parameters::Mesh::Type::gmsh)
    {
      if (mesh_parameters.simplex)
        {
          auto        comm      = triangulation->get_communicator();
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

          triangulation->create_triangulation(construction_data);
        }
      else
        {
          GridIn<dim, spacedim> grid_in;
          grid_in.attach_triangulation(*triangulation);
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
          const int initial_refinement = mesh_parameters.initial_refinement;
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
            Utilities::MPI::n_mpi_processes(triangulation->get_communicator()),
            temporary_tri_triangulation);
          GridTools::partition_multigrid_levels(temporary_tri_triangulation);

          // extract relevant information from distributed triangulation
          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation(
              temporary_tri_triangulation,
              triangulation->get_communicator(),
              TriangulationDescription::Settings::
                construct_multigrid_hierarchy);
          triangulation->create_triangulation(construction_data);
        }
      else
        {
          GridGenerator::generate_from_name_and_arguments(
            *triangulation,
            mesh_parameters.grid_type,
            mesh_parameters.grid_arguments);
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
          grid.make_grid(*triangulation);
        }
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");

  // Translate the mesh
  if (mesh_parameters.translate)
    throw std::runtime_error(
      "Main grid can't be moved: the solid mesh must be moved instead. Grid will not be moved.");
  // Rotate the mesh
  if (mesh_parameters.rotate)
    throw std::runtime_error(
      "Main grid can't be rotated: the solid mesh must be rotated instead. Grid will not be rotated.");


  // Setup periodic boundary conditions
  for (unsigned int i_bc = 0; i_bc < boundary_conditions.size; ++i_bc)
    {
      if (boundary_conditions.type[i_bc] ==
          BoundaryConditions::BoundaryType::periodic)
        {
          std::vector<GridTools::PeriodicFacePair<
            typename Triangulation<dim, spacedim>::cell_iterator>>
            periodicity_vector;
          GridTools::collect_periodic_faces(
            *dynamic_cast<Triangulation<dim, spacedim> *>(triangulation.get()),
            boundary_conditions.id[i_bc],
            boundary_conditions.periodic_id[i_bc],
            boundary_conditions.periodic_direction[i_bc],
            periodicity_vector);
          triangulation->add_periodicity(periodicity_vector);
        }
    }
}

template <int dim, int spacedim>
void
read_mesh_and_manifolds(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim, spacedim>>
                                                          triangulation,
  const Parameters::Mesh &                                mesh_parameters,
  const Parameters::Manifolds &                           manifolds_parameters,
  const bool &                                            restart,
  const BoundaryConditions::BoundaryConditions<spacedim> &boundary_conditions)
{
  attach_grid_to_triangulation(triangulation,
                               mesh_parameters,
                               boundary_conditions);
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
            GridTools::minimal_cell_diameter(*triangulation);
          double       target_size = mesh_parameters.target_size;
          unsigned int number_refinement =
            floor(std::log(minimal_cell_size / target_size) / std::log(2));
          triangulation->refine_global(number_refinement);
        }
      else if (!restart)
        {
          const int initial_refinement = mesh_parameters.initial_refinement;
          triangulation->refine_global(initial_refinement);
        }
    }
}

template void attach_grid_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<2>> triangulation,
  const Parameters::Mesh &                                   mesh_parameters,
  const BoundaryConditions::BoundaryConditions<2> &boundary_conditions);
template void attach_grid_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<3>> triangulation,
  const Parameters::Mesh &                                   mesh_parameters,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);
template void attach_grid_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<2, 3>> triangulation,
  const Parameters::Mesh &                                      mesh_parameters,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);

template void read_mesh_and_manifolds(
  std::shared_ptr<parallel::DistributedTriangulationBase<2>> triangulation,
  const Parameters::Mesh &                                   mesh_parameters,
  const Parameters::Manifolds &                    manifolds_parameters,
  const bool &                                     restart,
  const BoundaryConditions::BoundaryConditions<2> &boundary_conditions);
template void read_mesh_and_manifolds(
  std::shared_ptr<parallel::DistributedTriangulationBase<3>> triangulation,
  const Parameters::Mesh &                                   mesh_parameters,
  const Parameters::Manifolds &                    manifolds_parameters,
  const bool &                                     restart,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);
template void read_mesh_and_manifolds(
  std::shared_ptr<parallel::DistributedTriangulationBase<2, 3>> triangulation,
  const Parameters::Mesh &                                      mesh_parameters,
  const Parameters::Manifolds &                    manifolds_parameters,
  const bool &                                     restart,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);
