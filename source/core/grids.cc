#include "core/grids.h"



template <int dim>
void
attach_grid_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  const Parameters::Mesh &                                     mesh_parameters)
{
  // GMSH input
  if (mesh_parameters.type == Parameters::Mesh::Type::gmsh)
    {
      GridIn<dim> grid_in;
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
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");

  // const int initialSize = mesh_parameters.initialRefinement;

  // Refine the triangulation to its initial size
  // triangulation->refine_global(initialSize);
}

// attach_grid_to_triangulation(
//  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
//  const Parameters::Mesh &mesh_parameters)
//{
//// Setup parallelism for periodic boundary conditions
// for (unsigned int i_bc = 0; i_bc < nsparam.boundaryConditions.size; ++i_bc)
//  {
//    if (nsparam.boundaryConditions.type[i_bc] ==
//        BoundaryConditions::BoundaryType::periodic)
//      {
//        std::vector<GridTools::PeriodicFacePair<
//          typename Triangulation<dim>::cell_iterator>>
//          periodicity_vector;
//        GridTools::collect_periodic_faces(
//          *dynamic_cast<Triangulation<dim> *>(this->triangulation.get()),
//          nsparam.boundaryConditions.id[i_bc],
//          nsparam.boundaryConditions.periodic_id[i_bc],
//          nsparam.boundaryConditions.periodic_direction[i_bc],
//          periodicity_vector);
//        this->triangulation->add_periodicity(periodicity_vector);
//      }
//  }
//}


template void attach_grid_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<2>> triangulation,
  const Parameters::Mesh &                                   mesh_parameters);

template void attach_grid_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<3>> triangulation,
  const Parameters::Mesh &                                   mesh_parameters);
