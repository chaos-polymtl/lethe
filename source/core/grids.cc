// Deal.II includes
#include <deal.II/grid/grid_tools.h>

// Lethe includes
#include "core/boundary_conditions.h"
#include "core/grids.h"

// Per_Hills needed includes
#include "core/per_hills_grid.h"

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
      GridIn<dim, spacedim> grid_in;
      grid_in.attach_triangulation(*triangulation);
      std::ifstream input_file(mesh_parameters.file_name);
      grid_in.read_msh(input_file);
    }

  // Dealii grids
  else if (mesh_parameters.type == Parameters::Mesh::Type::dealii)
    {
      GridGenerator::generate_from_name_and_arguments(
        *triangulation,
        mesh_parameters.grid_type,
        mesh_parameters.grid_arguments);
    }

  // Periodic Hills grid
  else if (mesh_parameters.type == Parameters::Mesh::Type::per_hills)
    {
      PeriodicHillsGrid<dim, spacedim> grid(mesh_parameters.grid_arguments);
      grid.make_grid(*triangulation);
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");


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
  const BoundaryConditions::BoundaryConditions<spacedim> &boundary_conditions)
{
  attach_grid_to_triangulation(triangulation,
                               mesh_parameters,
                               boundary_conditions);
  attach_manifolds_to_triangulation(triangulation, manifolds_parameters);

  // Refine the triangulation to its initial size
  const int initialSize = mesh_parameters.initial_refinement;
  triangulation->refine_global(initialSize);
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
  const BoundaryConditions::BoundaryConditions<2> &boundary_conditions);
template void read_mesh_and_manifolds(
  std::shared_ptr<parallel::DistributedTriangulationBase<3>> triangulation,
  const Parameters::Mesh &                                   mesh_parameters,
  const Parameters::Manifolds &                    manifolds_parameters,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);
template void read_mesh_and_manifolds(
  std::shared_ptr<parallel::DistributedTriangulationBase<2, 3>> triangulation,
  const Parameters::Mesh &                                      mesh_parameters,
  const Parameters::Manifolds &                    manifolds_parameters,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);
