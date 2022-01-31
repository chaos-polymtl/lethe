//
// Created by Emile Bergeron on 2022-01-30.
//

#include <core/pvd_handler.h>

#include <dem/dem_solver_parameters.h>
#include <dem/grid_motion.h>


using namespace dealii;

#ifndef floating_grid_h
#define floating_grid_h

template <int dim>
class FloatingGrid
{
public:
  FloatingGrid(const DEMSolverParameters<dim> &dem_parameters,
               const ConditionalOStream &      pcout,
               const double &                  dem_time_step);

  void iterate();

  void write(const std::string                      folder,
             const std::string                      file_prefix,
             const double                           time,
             const unsigned int                     iter,
             const unsigned int                     group_files,
             const MPI_Comm &                       mpi_communicator,
             const unsigned int                     digits = 4);

private:
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  GridMotion<dim> gridMotion;
  double triangulation_cell_diameter;
  PVDHandler pvdHandler;
};

#endif // floating_grid_h
