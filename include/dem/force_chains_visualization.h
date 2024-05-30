#include <core/auxiliary_math_functions.h>
#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_particle_contact_force.h>
#include <dem/particle_particle_contact_info.h>
#include <dem/rolling_resistance_torque_models.h>

#include <deal.II/particles/particle_handler.h>

#include <boost/range/adaptor/map.hpp>

#include <vector>

using namespace dealii;
using namespace DEM;

#ifndef particles_force_chains_v2_h
#  define particles_force_chains__v2_h



template <
  int                                                       dim,
  Parameters::Lagrangian::ParticleParticleContactForceModel force_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
class ParticlesForceChains
  : public ParticleParticleContactForce<dim,
                                        force_model,
                                        rolling_friction_model>
{
public:
  ParticlesForceChains(const DEMSolverParameters<dim> &dem_parameters);

  virtual ~ParticlesForceChains()
  {}

  void
  multi_general_cell(Triangulation<1, 3>         &tria,
                     const std::vector<Point<3>> &vertices);

  void
  calculate_force_chains(
    DEMContactManager<dim>    &container_manager,
    const double               dt);


  void
  write_force_chains(MPI_Comm           mpi_communicator,
                     const std::string  folder,
                     const unsigned int iter);

private:
  std::vector<double>             force_normal;
  std::vector<Point<3>>           vertices;
  const DEMSolverParameters<dim> &dem_parameters;
};


#endif
