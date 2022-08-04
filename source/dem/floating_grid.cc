#include "dem/floating_grid.h"

#include <core/solutions_output.h>

#include "dem/read_mesh.h"

#include <deal.II/grid/grid_out.h>

template <int dim, int spacedim>
FloatingGrid<dim, spacedim>::FloatingGrid(
  const Parameters::Lagrangian::FloatingGrid<spacedim>
    &                       floating_grid_parameters,
  const bool                restart,
  const ConditionalOStream &pcout,
  const double &            dem_time_step)
  : gridMotion(floating_grid_parameters.motion, dem_time_step)
{
  if (floating_grid_parameters.mesh.file_name != "none")
    {
      // Default boundary condition parameters
      Parameters::Lagrangian::BCDEM bc_parameters;
      bc_parameters.BC_type = Parameters::Lagrangian::BCDEM::BoundaryType::fixed_wall;

      read_mesh(floating_grid_parameters.mesh,
                restart,
                pcout,
                triangulation,
                triangulation_cell_diameter,
                bc_parameters);
    }
}

template <int dim, int spacedim>
void
FloatingGrid<dim, spacedim>::iterate()
{
  gridMotion.move_grid(triangulation);
}

template <int dim, int spacedim>
void
FloatingGrid<dim, spacedim>::write(const std::string  folder,
                                   const std::string  file_prefix,
                                   const double       time,
                                   const unsigned int iter,
                                   const unsigned int group_files,
                                   const MPI_Comm &   mpi_communicator,
                                   const unsigned int digits)
{
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
  DataOut<dim, DoFHandler<dim, spacedim>> dataOut;
#else
  DataOut<dim, spacedim> dataOut;
#endif
  dataOut.attach_triangulation(triangulation);
  dataOut.build_patches();

  write_vtu_and_pvd<dim, spacedim>(pvdHandler,
                                   dataOut,
                                   folder,
                                   file_prefix,
                                   time,
                                   iter,
                                   group_files,
                                   mpi_communicator,
                                   digits);
}

template class FloatingGrid<1, 2>;
template class FloatingGrid<2, 2>;
template class FloatingGrid<2, 3>;
template class FloatingGrid<3, 3>;