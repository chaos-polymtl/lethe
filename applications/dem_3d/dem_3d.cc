

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
#include "dem/dem.h"
#include "dem/dem_iterator.h"
#include "dem/dem_properties.h"
#include "dem/dem_solver_parameters.h"
#include "dem/find_boundary_cells_information.h"
#include "dem/find_cell_neighbors.h"
#include "dem/insertion_info_struct.h"
#include "dem/particle_wall_contact_detection.h"
#include "dem/physical_info_struct.h"
#include "dem/pp_broad_search.h"
#include "dem/pp_fine_search.h"
#include "dem/pp_linear_force.h"
#include "dem/pp_nonlinear_force.h"
#include "dem/pw_broad_search.h"
#include "dem/velocity_verlet_integrator.h"

#include "variant"
using namespace dealii;

template <int dim, int spacedim> void initilization() {
  std::string filename;
  filename = "dem.prm";
  ParameterHandler prm;
  DEMSolverParameters<dim> DEMparam;
  DEMparam.declare(prm);
  prm.parse_input(filename);
  DEMparam.parse(prm);

  DEM::DEMProperties<3> properties_class;
  std::vector<std::tuple<std::string, int>> properties =
      properties_class.get_properties_name();

  // DEM clock
  int DEM_step = 0;
  float DEM_time = 0;
  parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr, -0.1, 0.1, true);
  int numRef = 3;
  tr.refine_global(numRef);
  auto &mappinggg = StaticMappingQ1<dim, spacedim>::mapping;

  Particles::ParticleHandler<dim, spacedim> particle_handler(
      tr, mappinggg, DEMparam.outputProperties.numProperties);
  Particles::PropertyPool property_pool(
      DEMparam.outputProperties.numProperties);

  std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
      cellNeighbor;

  FindCellNeighbors<dim, spacedim> cn1;
  cellNeighbor = cn1.find_cell_neighbors(tr);

  DEM_iterator<dim, spacedim> iter1;
  std::vector<std::map<int, Particles::ParticleIterator<dim, spacedim>>>
      inContactPairs(DEMparam.simulationControl.nTotal);
  std::vector<std::map<int, contact_info_struct<dim, spacedim>>> inContactInfo(
      DEMparam.simulationControl.nTotal);
  std::vector<
      std::tuple<Particles::ParticleIterator<dim, spacedim>, Tensor<1, dim>,
                 Point<dim>, double, double, double, Point<dim>, double>>
      pwContactInfo;

  ParticleWallContactDetection<dim, spacedim> pw1;

  std::vector<boundary_cells_info_struct<dim>> boundary_cells_information;
  FindBoundaryCellsInformation<dim, dim> boundary_cell_object;
  boundary_cells_information =
      boundary_cell_object.find_boundary_cells_information(tr);

  ParticleWallContactForce<dim, spacedim> pwcf1;
  VelocityVerletIntegrator<dim, spacedim> integ1;
  PPBroadSearch<dim, spacedim> ppbs;
  PPFineSearch<dim, spacedim> ppfs;
  // PPLinearForce<dim, spacedim> pplf;
  PPNonLinearForce<dim, spacedim> ppnlf;
  PWBroadSearch<dim, spacedim> pwbs;

  // reading parameters from parameter handler file:
  const int numberOfSteps = DEMparam.simulationControl.tFinal;
  const double dt = DEMparam.simulationControl.dt;
  const int nTotal = DEMparam.simulationControl.nTotal;
  const int writeFreq = DEMparam.simulationControl.writeFrequency;
  Tensor<1, dim> g{{DEMparam.physicalProperties.gx,
                    DEMparam.physicalProperties.gy,
                    DEMparam.physicalProperties.gz}};

  physical_info_struct<dim> physical_info_struct;
  physical_info_struct.particle_diameter = DEMparam.physicalProperties.diameter;
  physical_info_struct.particle_density = DEMparam.physicalProperties.density;
  physical_info_struct.Young_modulus_particle = DEMparam.physicalProperties.Yp;
  physical_info_struct.Young_modulus_wall = DEMparam.physicalProperties.Yw;
  physical_info_struct.restitution_coefficient_particle =
      DEMparam.physicalProperties.ep;
  physical_info_struct.restitution_coefficient_wall =
      DEMparam.physicalProperties.ew;
  physical_info_struct.Poisson_ratio_particle = DEMparam.physicalProperties.vp;
  physical_info_struct.Poisson_ratio_wall = DEMparam.physicalProperties.vw;
  physical_info_struct.friction_coefficient_particle =
      DEMparam.physicalProperties.mup;
  physical_info_struct.friction_coefficient_wall =
      DEMparam.physicalProperties.muw;
  physical_info_struct.rolling_friction_coefficient_particle =
      DEMparam.physicalProperties.murp;
  physical_info_struct.rolling_friction_coefficient_wall =
      DEMparam.physicalProperties.murw;

  insertion_info_struct<dim, spacedim> insertion_info_struct;
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

  const int numFields = DEMparam.outputProperties.numFields;
  const int numProperties = DEMparam.outputProperties.numProperties;

  // dem engine iterator:
  while (DEM_step < numberOfSteps) {
    iter1.engine(particle_handler, tr, DEM_step, DEM_time, cellNeighbor,
                 inContactPairs, inContactInfo, boundary_cells_information,
                 pwContactInfo, properties, property_pool, pw1, &ppnlf, pwcf1,
                 &integ1, dt, nTotal, writeFreq, physical_info_struct,
                 insertion_info_struct, g, numFields, numProperties, ppbs, ppfs,
                 pwbs);
  }
}

int main(int argc, char *argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
      argc, argv, numbers::invalid_unsigned_int);

  initilization<3, 3>();

  return 0;
}
