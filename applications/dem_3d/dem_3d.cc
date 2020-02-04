

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

#include "dem/contact_info_struct.h"
#include "dem/contact_search.h"
#include "dem/dem.h"
#include "dem/dem_iterator.h"
#include "dem/dem_properties.h"
#include "dem/dem_solver_parameters.h"
#include "dem/insertion_info_struct.h"
#include "dem/particle_wall_contact_detection.h"
#include "dem/velocity_verlet_integrator.h"

using namespace dealii;

template <int dim, int spacedim>
void
initilization()
{
  std::string filename;
  filename = "dem.prm";
  ParameterHandler         prm;
  DEMSolverParameters<dim> DEMparam;
  DEMparam.declare(prm);
  prm.parse_input(filename);
  DEMparam.parse(prm);

  DEM::DEMProperties<3>                     properties_class;
  std::vector<std::tuple<std::string, int>> properties =
    properties_class.get_properties_name();

  // total number of active particles in the system
  int active_particle_number = 0;
  // DEM clock
  int                                                 DEM_step = 0;
  float                                               DEM_time = 0;
  parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr, -0.1, 0.1, true);
  int numRef = 3;
  tr.refine_global(numRef);
  auto &mappinggg = StaticMappingQ1<dim, spacedim>::mapping;

  Particles::ParticleHandler<dim, spacedim> particle_handler(
    tr, mappinggg, DEMparam.outputProperties.numProperties);
  Particles::PropertyPool property_pool(
    DEMparam.outputProperties.numProperties);

  std::pair<
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>,
    std::vector<typename Triangulation<dim>::active_cell_iterator>>
                               cellNeighbor;
  ContactSearch<dim, spacedim> cs1;
  cellNeighbor = cs1.findCellNeighbors(tr);

  DEM_iterator<dim, spacedim> iter1;
  std::vector<std::map<int, Particles::ParticleIterator<dim, spacedim>>>
                                                               inContactPairs(DEMparam.simulationControl.nTotal);
  std::vector<std::map<int, ContactInfoStruct<dim, spacedim>>> inContactInfo(
    DEMparam.simulationControl.nTotal);
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
  VelocityVerletIntegrator<dim, spacedim> integ1;

  // reading parameters from parameter handler file:
  const int    numberOfSteps = DEMparam.simulationControl.tFinal;
  const double dt            = DEMparam.simulationControl.dt;
  const int    nTotal        = DEMparam.simulationControl.nTotal;
  const int    writeFreq     = DEMparam.simulationControl.writeFrequency;
  Point<dim>   g             = {DEMparam.physicalProperties.gx,
                  DEMparam.physicalProperties.gy,
                  DEMparam.physicalProperties.gz};
  const double dp            = DEMparam.physicalProperties.diameter;
  const int    rhop          = DEMparam.physicalProperties.density;
  const int    Yp            = DEMparam.physicalProperties.Yp;
  const int    Yw            = DEMparam.physicalProperties.Yw;
  const float  vp            = DEMparam.physicalProperties.vp;
  const float  vw            = DEMparam.physicalProperties.vw;
  const float  ep            = DEMparam.physicalProperties.ep;
  const float  ew            = DEMparam.physicalProperties.ew;
  const float  mup           = DEMparam.physicalProperties.mup;
  const float  muw           = DEMparam.physicalProperties.muw;
  const float  murp          = DEMparam.physicalProperties.murp;
  const float  murw          = DEMparam.physicalProperties.murw;

  InsertionInfoStruct<dim, spacedim> insertion_info_struct;
  insertion_info_struct.insertion_steps_number =
    DEMparam.insertionInfo.tInsertion;
  insertion_info_struct.inserted_number_at_step =
    DEMparam.insertionInfo.nInsert;
  insertion_info_struct.insertion_frequency =
    DEMparam.insertionInfo.insertFrequency;
  insertion_info_struct.distance_threshold =
    DEMparam.insertionInfo.distance_threshold;
  insertion_info_struct.x_min = DEMparam.insertionInfo.x_min;
  insertion_info_struct.y_min = DEMparam.insertionInfo.y_min;
  insertion_info_struct.z_min = DEMparam.insertionInfo.z_min;
  insertion_info_struct.x_max = DEMparam.insertionInfo.x_max;
  insertion_info_struct.y_max = DEMparam.insertionInfo.y_max;
  insertion_info_struct.z_max = DEMparam.insertionInfo.z_max;

  const int numFields     = DEMparam.outputProperties.numFields;
  const int numProperties = DEMparam.outputProperties.numProperties;

  // dem engine iterator:
  while (DEM_step < numberOfSteps)
    {
      iter1.engine(active_particle_number,
                   particle_handler,
                   tr,
                   DEM_step,
                   DEM_time,
                   cellNeighbor,
                   inContactPairs,
                   inContactInfo,
                   boundaryCellInfo,
                   pwContactInfo,
                   properties,
                   property_pool,
                   cs1,
                   pw1,
                   cf1,
                   pwcf1,
                   &integ1,
                   numberOfSteps,
                   dt,
                   nTotal,
                   writeFreq,
                   g,
                   dp,
                   rhop,
                   Yp,
                   Yw,
                   vp,
                   vw,
                   ep,
                   ew,
                   mup,
                   muw,
                   murp,
                   murw,
                   insertion_info_struct,
                   numFields,
                   numProperties);
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
