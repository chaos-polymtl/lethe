//
// Created by Emile Bergeron on 2022-01-30.
//
#include <core/solutions_output.h>
#include "dem/floating_grid.h"
#include "dem/read_mesh.h"

#include <deal.II/grid/grid_out.h>

template <int dim, int spacedim>
FloatingGrid<dim, spacedim>::FloatingGrid(const DEMSolverParameters<spacedim> &dem_parameters,
                                          const ConditionalOStream &     pcout,
                                          const double                   &dem_time_step)
  : gridMotion(dem_parameters, dem_time_step)
{
  read_mesh(dem_parameters,
            pcout,
            triangulation,
            triangulation_cell_diameter);
}

template <int dim, int spacedim>
void
FloatingGrid<dim, spacedim>::iterate()
{
  gridMotion.move_grid(triangulation);
}

template <int dim, int spacedim>
void
FloatingGrid<dim, spacedim>::write(const std::string                      folder,
                                   const std::string                      file_prefix,
                                   const double                           time,
                                   const unsigned int                     iter,
                                   const unsigned int                     group_files,
                                   const MPI_Comm &                       mpi_communicator,
                                   const unsigned int                     digits)
{
  DataOut<dim, spacedim> dataOut;
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