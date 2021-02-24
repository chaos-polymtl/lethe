// Deal.II includes
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe includes
#include "core/boundary_conditions.h"
#include "core/grids.h"
#include "core/periodic_hills_grid.h"

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
  if (mesh_parameters.type == Parameters::Mesh::Type::gmsh &&
      !mesh_parameters.simplex)
    {
      GridIn<dim, spacedim> grid_in;
      grid_in.attach_triangulation(*triangulation);
      std::ifstream input_file(mesh_parameters.file_name);
      grid_in.read_msh(input_file);
    }
#ifdef DEAL_II_WITH_SIMPLEX_SUPPORT
  // Simplex && GMSH input
  else if (mesh_parameters.type == Parameters::Mesh::Type::gmsh &&
           mesh_parameters.simplex)
    {
      // create serial triangulation
      Triangulation<dim, spacedim> basetria(
        Triangulation<dim, spacedim>::limit_level_difference_at_vertices);
      GridIn<dim, spacedim> grid_in;
      grid_in.attach_triangulation(basetria);
      std::ifstream input_file(mesh_parameters.file_name);
      grid_in.read_msh(input_file);

      GridTools::partition_triangulation_zorder(
        Utilities::MPI::n_mpi_processes(triangulation->get_communicator()),
        basetria);
      GridTools::partition_multigrid_levels(basetria);

      // create parallel::fullydistributed::Triangulation from description
      auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(
          basetria,
          triangulation->get_communicator(),
          TriangulationDescription::Settings::construct_multigrid_hierarchy);
      triangulation->create_triangulation(construction_data);

      for (auto &cell : triangulation->active_cell_iterators())
        {
          CellId id        = cell->id();
          auto   cell_base = basetria.create_cell_iterator(id);
          // Assert(cell->center() == cell_base->center(),
          //       ExcMessage("Cells do not match"));
          for (unsigned int d = 0; d < dim; d++)
            Assert(std::abs(cell->center()[d] - cell_base->center()[d]) < 1e-9,
                   ExcMessage("Cells do not match"));
        }
    }
#endif
  // Dealii grids
  else if (mesh_parameters.type == Parameters::Mesh::Type::dealii)
    {
      GridGenerator::generate_from_name_and_arguments(
        *triangulation,
        mesh_parameters.grid_type,
        mesh_parameters.grid_arguments);
    }

  // Periodic Hills grid
  else if (mesh_parameters.type == Parameters::Mesh::Type::periodic_hills)
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
  const bool &                                            restart,
  const BoundaryConditions::BoundaryConditions<spacedim> &boundary_conditions)
{
  attach_grid_to_triangulation(triangulation,
                               mesh_parameters,
                               boundary_conditions);
  attach_manifolds_to_triangulation(triangulation, manifolds_parameters);

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
