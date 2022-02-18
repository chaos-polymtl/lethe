//
// Created by lucka on 2021-10-11.
//

#ifndef LETHE_IB_PARTICLES_DEM_H
#define LETHE_IB_PARTICLES_DEM_H

#include <core/ib_particle.h>
#include <core/ib_stencil.h>
#include <core/lethegridtools.h>

#include <solvers/navier_stokes_base.h>

#include <dem/particle_particle_contact_force.h>
#include <dem/particle_particle_contact_info_struct.h>
#include <dem/particle_wall_contact_info_struct.h>

#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

using namespace dealii;


/**
 * A solver class for the DEM used in conjunction with IB particles and
 * gls_sharp_navier_stokes. This class defines and uses some functions  of the
 * DEM class that has been modified and simplified to use with IB_particles.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 */

template <int dim>
class IBParticlesDEM
{
public:
  /**
   * @brief
   * Initialize the IBParticlesDEM object with the parameters, the mpi
   * communicator and the particles.
   *
   * @param p_nsparam The parameters' object of the simulation.
   * @param dem_parameters
   * @param mpi_communicator_input The mpi communicator of the simulation.
   *
   * @param particles The particles vector containing all the IB particles.
   */
  void
  initialize(std::shared_ptr<Parameters::IBParticles<dim>> &p_nsparam,
             MPI_Comm &                   mpi_communicator_input,
             std::vector<IBParticle<dim>> particles);


  /**
   * @brief
   * Updates the boundary cells that are contact candidates for each of the
   * particles.
   *
   * @param particles The particles vector containing all the IB particles.
   *
   * @param time The current CFD time.
   */
  void
  update_particles(std::vector<IBParticle<dim>> particles, double &time);


  /**
   * @brief
   * Integrates the dynamics of the IB_particle taking into account the contact
   * between particles and between particles and walls.
   * @param dt The CFD time step.
   *
   */
  void
  particles_dem(double &dt);

  /**
   * @brief Calculates non-linear (Hertzian) particle-particle contact force
   *
   * @param dt_dem The sub time stepping time step.
   *
   * @param contact_force a vector containing the contact force between particles
   *
   * @param contact_force a vector containing the contact torques between particles
   */
  void
  calculate_pp_contact_force(const double &             dt_dem,
                             std::vector<Tensor<1, 3>> &contact_force,
                             std::vector<Tensor<1, 3>> &contact_torque);


  /**
   * @brief Calculates non-linear (Hertzian) particle-wall contact force
   *
   * @param dt_dem The sub time stepping time step.
   *
   * @param contact_force a vector containing the contact force between particles
   *
   * @param contact_force a vector containing the contact torques between particles
   */
  void
  calculate_pw_contact_force(const double &             dt_dem,
                             std::vector<Tensor<1, 3>> &contact_force,
                             std::vector<Tensor<1, 3>> &contact_torque);

  /**
   * @brief  Updates the boundary cells that are contact candidates for each of the particles.
   *
   * @param particles The particles vector containing all the IB particles.
   *
   * @param dof_handler The dof handler of the mesh used for the fluid simulation.
   *
   * @param face_quadrature_formula The face quadrature formula used in the elements.
   *
   * @param mapping The FEM mapping of the face element.
   */

  void
  update_particles_boundary_contact(
    std::vector<IBParticle<dim>> &particles,
    DoFHandler<dim> &             dof_handler,
    const Quadrature<dim - 1> &   face_quadrature_formula,
    const Mapping<dim> &          mapping);


  std::vector<IBParticle<dim>> dem_particles;

private:
  // A struct to store boundary cells' information
  struct BoundaryCellsInfo
  {
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
      for (unsigned int i = 0; i < dim; ++i)
        {
          ar &normal_vector;
          ar &point_on_boundary;
        }
    }

    Tensor<1, dim> normal_vector;
    Point<dim>     point_on_boundary;
  };

  std::shared_ptr<Parameters::IBParticles<dim>> parameters;
  DEMSolverParameters<dim>                      dem_parameters{};
  MPI_Comm                                      mpi_communicator;

  std::shared_ptr<ParticleParticleContactForce<dim>>
    particle_particle_contact_force_object;

  std::shared_ptr<ParticleWallContactForce<dim>>
    particle_wall_contact_force_object;

  // Empty parameters to initilize particle_wall_contact_force_object
  double triangulation_cell_diameter = 0.0;

  // Particles contact history
  std::map<unsigned int,
           std::map<unsigned int, particle_particle_contact_info_struct<dim>>>
    pp_contact_map;
  std::map<unsigned int,
           std::map<unsigned int, particle_wall_contact_info_struct<dim>>>
    pw_contact_map;

  // A vector of vectors of candidate cells for each of the particle.
  std::vector<std::vector<BoundaryCellsInfo>> boundary_cells;

private:
  double cfd_time;
};


#endif // LETHE_IB_PARTICLES_DEM_H
