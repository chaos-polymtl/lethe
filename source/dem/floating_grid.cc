//
// Created by Emile Bergeron on 2022-01-30.
//
#include <core/solutions_output.h>

#include "dem/floating_grid.h"
#include "dem/read_mesh.h"

#include <deal.II/grid/grid_out.h>

template <int dim>
FloatingGrid<dim>::FloatingGrid(const DEMSolverParameters<dim> &dem_parameters,
                                const ConditionalOStream &     pcout,
                                const double                   &dem_time_step)
{
  read_mesh(dem_parameters,
            pcout,
            triangulation,
            triangulation_cell_diameter);
  gridMotion = GridMotion<dim>(dem_parameters, dem_time_step);
}

template <int dim>
void
FloatingGrid<dim>::iterate()
{
  gridMotion.move_grid(triangulation);
}

template <int dim>
void
FloatingGrid<dim>::write(const std::string                      folder,
                         const std::string                      file_prefix,
                         const double                           time,
                         const unsigned int                     iter,
                         const unsigned int                     group_files,
                         const MPI_Comm &                       mpi_communicator,
                         const unsigned int                     digits)
{
  DataOut<dim> dataOut;
  dataOut.attach_triangulation(triangulation);
  dataOut.build_patches();

  write_vtu_and_pvd<dim>(pvdHandler,
                         dataOut,
                         folder,
                         file_prefix,
                         time,
                         iter,
                         group_files,
                         mpi_communicator,
                         digits);
}