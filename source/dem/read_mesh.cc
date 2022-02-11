#include <dem/read_mesh.h>

template <int dim, int spacedim=dim>
void
read_mesh(const DEMSolverParameters<spacedim> &           parameters,
          const ConditionalOStream &                 pcout,
          Triangulation<dim, spacedim> &triangulation,
          double &triangulation_cell_diameter)
{
  pcout << "Reading triangulation " << std::endl;
  // GMSH input
  if (parameters.mesh.type == Parameters::Mesh::Type::gmsh)
    {
      GridIn<dim, spacedim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream input_file(parameters.mesh.file_name);
      grid_in.read_msh(input_file);
    }

  // Dealii grids
  else if (parameters.mesh.type == Parameters::Mesh::Type::dealii)
    {
      GridGenerator::generate_from_name_and_arguments(
        triangulation,
        parameters.mesh.grid_type,
        parameters.mesh.grid_arguments);
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");

  triangulation_cell_diameter = 0.5 * GridTools::diameter(triangulation);

  if (parameters.restart.restart == false)
    {
      if (parameters.mesh.refine_until_target_size)
        {
          double minimal_cell_size =
            GridTools::minimal_cell_diameter(triangulation);
          double       target_size = parameters.mesh.target_size;
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
          const int initial_refinement = parameters.mesh.initial_refinement;
          triangulation.refine_global(initial_refinement);
        }
    }

  pcout << std::endl << "Finished reading triangulation " << std::endl;
}

template void
read_mesh<1, 2>(const DEMSolverParameters<2> &           parameters,
                const ConditionalOStream &               pcout,
                Triangulation<1, 2> &triangulation,
                double &                                 triangulation_cell_diameter);

template void
read_mesh<2, 2>(const DEMSolverParameters<2> &           parameters,
                const ConditionalOStream &               pcout,
                Triangulation<2, 2> &triangulation,
                double &                                 triangulation_cell_diameter);

template void
read_mesh<2, 3>(const DEMSolverParameters<3> &           parameters,
                const ConditionalOStream &               pcout,
                Triangulation<2, 3> &triangulation,
                double &                                 triangulation_cell_diameter);

template void
read_mesh<3, 3>(const DEMSolverParameters<3> &           parameters,
                const ConditionalOStream &               pcout,
                Triangulation<3, 3> &triangulation,
                double &                                 triangulation_cell_diameter);
