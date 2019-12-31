

#include <deal.II/base/array_view.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <math.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "dem/contact_search.h"
#include "dem/dem_iterator.h"
#include "dem/integration.h"
#include "dem/parameters_dem.h"
#include "dem/particle_insertion.h"
#include "dem/particle_wall_contact_detection.h"



using namespace dealii;


template <int dim, int spacedim>
void
initilization()
{
  std::string filename;
   filename = "dem.prm";
  ParameterHandler   prm;
  ParametersDEM<dim> DEMparam;
  DEMparam.declare(prm);
  prm.parse_input(filename);
  DEMparam.parse(prm);

  std::vector<std::tuple<std::string, int>> properties(
    DEMparam.outputProperties.numProperties);
  properties[0]  = std::make_tuple("id", 1);
  properties[1]  = std::make_tuple("type", 1);
  properties[2]  = std::make_tuple("Diam", 1);
  properties[3]  = std::make_tuple("Dens", 1);
  properties[4]  = std::make_tuple("Pos", dim);
  properties[5]  = std::make_tuple("Pos", 1);
  properties[6]  = std::make_tuple("Pos", 1);
  properties[7]  = std::make_tuple("Vel", dim);
  properties[8]  = std::make_tuple("Vel", 1);
  properties[9]  = std::make_tuple("Vel", 1);
  properties[10] = std::make_tuple("a", dim);
  properties[11] = std::make_tuple("a", 1);
  properties[12] = std::make_tuple("a", 1);
  properties[13] = std::make_tuple("f", dim);
  properties[14] = std::make_tuple("f", 1);
  properties[15] = std::make_tuple("f", 1);
  properties[16] = std::make_tuple("w", dim);
  properties[17] = std::make_tuple("w", 1);
  properties[18] = std::make_tuple("w", 1);
  properties[19] = std::make_tuple("mass", 1);
  properties[20] = std::make_tuple("MOI", 1);
  properties[21] = std::make_tuple("M", dim);
  properties[22] = std::make_tuple("M", 1);
  properties[23] = std::make_tuple("M", 1);


  // total number of particles in the system
  int nPart = 0;
  // DEM clock
  int   DEM_step = 0;
  float DEM_time = 0;
  parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

  Particles::Particle<dim> particle;
  GridGenerator::hyper_cube(tr, -0.1, 0.1, true);
  int numRef = 3;
  tr.refine_global(numRef);
  auto &mappinggg = StaticMappingQ1<dim, spacedim>::mapping;

  Particles::ParticleHandler<dim, spacedim> particle_handler(
    tr, mappinggg, DEMparam.outputProperties.numProperties);
  Particles::PropertyPool propPool(DEMparam.outputProperties.numProperties);

  std::pair<
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>,
    std::vector<typename Triangulation<dim>::active_cell_iterator>>
                               cellNeighbor;
  ContactSearch<dim, spacedim> cs1;
  cellNeighbor = cs1.findCellNeighbors(tr);


  DEM_iterator<dim, spacedim> iter1;
  std::vector<std::tuple<std::pair<Particles::ParticleIterator<dim, spacedim>,
                                   Particles::ParticleIterator<dim, spacedim>>,
                         double,
                         Point<dim>,
                         double,
                         Point<dim>,
                         double,
                         double>>
    contactInfo;
  std::vector<
    std::tuple<std::pair<Particles::ParticleIterator<dim, spacedim>, int>,
               Point<dim>,
               Point<dim>,
               double,
               double,
               double,
               Point<dim>,
               double>>
    pwContactInfo;

  ContactForce<dim, spacedim> cf1;

  std::vector<std::tuple<int,
                         typename Triangulation<dim>::active_cell_iterator,
                         int,
                         Point<dim>,
                         Point<dim>>>
                                              boundaryCellInfo;
  ParticleWallContactDetection<dim, spacedim> pw1;
  pw1.boundaryCellsAndFaces(tr, boundaryCellInfo);

  ParticleWallContactForce<dim, spacedim> pwcf1;
  Integration<dim, spacedim>              integ1;

  // dem engine iterator:
  while (DEM_step < DEMparam.simulationControl.tFinal)
    {
      iter1.engine(nPart,
                   particle_handler,
                   tr,
                   DEM_step,
                   DEM_time,
                   DEMparam,
                   cellNeighbor,
                   contactInfo,
                   boundaryCellInfo,
                   pwContactInfo,
                   properties,
                   propPool,
                   cs1,
                   pw1,
                   cf1,
                   pwcf1,
                   integ1);
    }
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  initilization<3, 3>();

  return 0;
}
