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

#ifndef particles_force_chains_h
#  define particles_force_chains_h

/**
 * Base class for the particles force chains contact force models.
 * This class does not implement any of the models, but ensures that
 * an interface without template specialization is available. All of the
 * actual implementation of the models are carried out in the
 * ParticlesForceChains class which is templated by the contact model
 * type.
 */
template <int dim>
class ParticlesForceChainsBase
{
public:
  /**
   * @brief Use a ParticleParticleContactForce object to calculate normal forces between
   * all touching particles. Store normal forces and particles position in
   * vectors.
   *
   * @param container_manager The container manager object that contains
   * containers to modify of contact pair periodic candidates with other
   * containers with periodic neighbors lists
   */
  virtual void
  calculate_force_chains(DEMContactManager<dim> &container_manager) = 0;
  /**
   * @brief Output the force chains in a single vtu file for each iteration.
   *
   * @param mpi_communicator The mpi communicator
   * @param folder a string that contains the path where the results are to be saved
   * @param iter the iteration number associated with the file
   */
  virtual void
  write_force_chains(const MPI_Comm     mpi_communicator,
                     const std::string  folder,
                     const unsigned int iter) = 0;
};

/**
 * @brief Class that carries out the calculation of
 * particle-particle contact force and the visualization of force chains
 * by writing vtu files. Instead of using a inheritance hiearchy to distinguish
 * between the contact model, the class is templated with the type of force
 * model and rolling friction model. Consequently, the code for each combination
 * of force model is generated at compile time.
 *
 * @tparam dim The dimension of the problem
 * @tparam force_model The particle-particle contact force model
 * @tparam rolling_friction_model The rolling resistance model used
 */
template <
  int                                                       dim,
  Parameters::Lagrangian::ParticleParticleContactForceModel force_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
class ParticlesForceChains
  : public ParticlesForceChainsBase<dim>,
    public ParticleParticleContactForce<dim,
                                        force_model,
                                        rolling_friction_model>
{
public:
  ParticlesForceChains(const DEMSolverParameters<dim> &dem_parameters);

  virtual ~ParticlesForceChains()
  {}
  /**
   * @brief Create a triangulation with a unique cell (represented by a line)
   * for each normal force between two particles.
   *
   * @param tria Empty triangulation used to create a triangulation with all
   * the vertices needed.
   * @param vertices Vector of points used to create a triangulation
   */
  void
  multi_general_cell(Triangulation<1, 3>         &tria,
                     const std::vector<Point<3>> &vertices);
  /**
   * @brief Calculate normal forces between all touching particles with
   * ParticleParticleContactForce class' methods. Stock normal forces and
   * particles position in vectors.
   *
   * @param container_manager The container manager object that contains
   * containers to modify of contact pair periodic candidates with other
   * containers with periodic neighbors lists
   * @param dt DEM time step
   */
  void
  calculate_force_chains(DEMContactManager<dim> &container_manager) override;
  /**
   * @brief Output the force chains in a single vtu file for each iteration.
   *
   * @param mpi_communicator The mpi communicator
   * @param folder a string that contains the path where the results are to be saved
   * @param iter the iteration number associated with the file
   */
  void
  write_force_chains(const MPI_Comm     mpi_communicator,
                     const std::string  folder,
                     const unsigned int iter) override;

private:
  // vector of normal forces between each touching particles.
  std::vector<double> force_normal;
  // vector of positions of touching particles.
  std::vector<Point<3>> vertices;
};


#endif
