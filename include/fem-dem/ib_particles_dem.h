// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_ib_particles_dem_h
#define lethe_ib_particles_dem_h

#include <core/ib_particle.h>

#include <dem/particle_particle_contact_force.h>
#include <dem/particle_wall_contact_force.h>

#include <deal.II/base/tensor.h>

using namespace dealii;

namespace Parameters
{
  template <int dim>
  class IBParticles;
  namespace Lagrangian
  {
    template <int dim>
    class FloatingWalls;
    struct BCDEM;
  } // namespace Lagrangian
} // namespace Parameters

/**
 * @brief A solver class for the DEM used in conjunction with IB particles and
 * fluid_dynamics_sharp. This class defines and uses some functions of the
 * DEM class that have been modified and simplified to be compatible with
 * IB_particles.
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
  /** @brief Struct that contains history of the contact between two objects. This is the commonly used
   * constructor since it houses all the information required to perform the
   * contact calculation.
   * */
  struct ContactInfo
  {
    Tensor<1, 3> normal_vector;
    Point<3>     contact_point;
    double       normal_overlap;
    double       normal_relative_velocity;
    Tensor<1, 3> tangential_displacement;
    Tensor<1, 3> rolling_resistance_spring_torque;
    std::vector<Tensor<1, 3>>
      tangential_relative_velocity; // keep each step of RK4 in memory
  };

  /** @brief Struct is used to reduce the number of parameters passed in the function calculate_force_model.
   * */
  struct ObjectProperties
  {
    double object_mass;
    double object_radius;
    double object_youngs_modulus;
    double object_poisson_ratio;
    double object_restitution_coefficient;
    double object_friction_coefficient;
    double object_rolling_friction_coefficient;
  };

  /**
   * @brief
   * Initialize the IBParticlesDEM object with the parameters, the mpi
   * communicator and the particles.
   *
   * @param p_nsparam The parameters for the immersed boundary particles
   * @param floating_walls_parameters The parameters for the floating walls.
   * @param boundary_conditions_parameters The parameters for the DEM boundary conditions.
   * @param mpi_communicator_input The mpi communicator of the simulation.
   * @param particles The particles vector containing all the IB particles.
   */
  void
  initialize(const std::shared_ptr<Parameters::IBParticles<dim>> &p_nsparam,
             const std::shared_ptr<Parameters::Lagrangian::FloatingWalls<dim>>
                                                 floating_walls_parameters,
             const std::shared_ptr<Parameters::Lagrangian::BCDEM>
               &boundary_conditions_parameters,
             const MPI_Comm                     &mpi_communicator_input,
             const std::vector<IBParticle<dim>> &particles);


  /**
   * @brief
   * Updates the boundary cells that are contact candidates for each of the
   * particles.
   *
   * @param particles The particles vector containing all the IB particles.
   * @param time The current CFD time.
   */
  void
  update_particles(const std::vector<IBParticle<dim>> &particles, double time);


  /**
   * @brief
   * Integrates the dynamics of the IB_particle taking into account the contact
   * between particles and between particles and walls.
   *
   * @param dt The CFD time step.
   * @param h_max The gap between the two particle.
   * @param h_min The minimal gap distance considered in the force calculation to avoid the singularity of the model.
   * @param rho The fluid density.
   * @param mu The fluid viscosity.
   */
  void
  integrate_particles_motion(const double dt,
                             const double h_max,
                             const double h_min,
                             const double rho,
                             const double mu);

  /**
   * @brief Calculates non-linear (Hertzian) particle-particle contact force
   *
   * @param dt_dem The sub time stepping time step.
   * @param contact_force a vector containing the contact force between particles
   * @param contact_torque a vector containing the contact torques between particles
   */
  void
  calculate_pp_contact_force(const double               dt_dem,
                             std::vector<Tensor<1, 3>> &contact_force,
                             std::vector<Tensor<1, 3>> &contact_torque);

  /**
   * @brief Updates the contact candidates of all the particles
   */
  void
  update_contact_candidates();


  /**
   * @brief Calculates non-linear (Hertzian) particle-wall contact force
   *
   * @param dt_dem The sub time stepping time step.
   * @param contact_force a vector containing the contact force between particles
   * @param contact_torque a vector containing the contact torques between particles
   */
  void
  calculate_pw_contact_force(const double               dt_dem,
                             std::vector<Tensor<1, 3>> &contact_force,
                             std::vector<Tensor<1, 3>> &contact_torque);


  /**
   * @brief Calculates non-linear (Hertzian) force between two objects. Take into input the physical properties and contact state of the two objects.
   *
   * @param normal_overlap Contact normal overlap
   * @param contact_info A container that contains the required information for
   * calculation of the contact force for a particle pair in contact
   * @param contact_point The contact point between the two objects
   * @param contact_normal Contact normal unit vector
   * @param normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param rolling_resistance_torque Contact rolling resistance torque
   * @param particle_one_position the position of particle one
   * @param particle_one_velocity the velocity of particle one
   * @param particle_one_omega the angular velocity of particle one
   * @param particle_one_properties the physical properties of particle one
   * @param particle_two_position the position of particle two
   * @param particle_two_velocity the velocity of particle two
   * @param particle_two_omega the angular velocity of particle two
   * @param particle_two_properties the physical properties of particle two
   * @param dt the dem time step.
   */
  void
  calculate_force_model(const double           normal_overlap,
                        ContactInfo           &contact_info,
                        Point<3>              &contact_point,
                        Tensor<1, 3>          &contact_normal,
                        Tensor<1, 3>          &normal_force,
                        Tensor<1, 3>          &tangential_force,
                        Tensor<1, 3>          &rolling_resistance_torque,
                        Point<dim>            &particle_one_position,
                        Tensor<1, 3>          &particle_one_velocity,
                        Tensor<1, 3>          &particle_one_omega,
                        const ObjectProperties particle_one_properties,
                        Point<dim>            &particle_two_position,
                        Tensor<1, 3>          &particle_two_velocity,
                        Tensor<1, 3>          &particle_two_omega,
                        const ObjectProperties particle_two_properties,
                        const double           dt);


  /**
   *  @brief Calculates particle-particle lubrication force. The force is based on the formula from
   *  Microhydrodynamics: Principles and Selected Applications by Kim, Sangtae;
   * Karrila, Seppo J. ISBN 13: 9780750691734
   *
   * @param dt_dem The sub time stepping time step.
   * @param h_max The gap between the two particle.
   * @param h_min The minimal gap distance considered in the force calculation to avoid the singularity of the model.
   * @param mu The fluid viscosity.
   * @param lubrication_force a vector containing the lubrication force on the particles.
   * @param lubrication_torque a vector containing the lubrication torques on the particles.
   */
  void
  calculate_pp_lubrication_force(const double               dt_dem,
                                 const double               h_max,
                                 const double               h_min,
                                 const double               mu,
                                 std::vector<Tensor<1, 3>> &lubrication_force,
                                 std::vector<Tensor<1, 3>> &lubrication_torque);


  /**
   * @brief Calculates particle-wall lubrication force
   *
   * @param dt_dem The sub time stepping time step.
   * @param h_max The gap between particle and the wall below which we evaluate the force.
   * @param h_min The minimal gap distance considered in the force calculation to avoid the singularity of the model.
   * @param mu The fluid viscosity.
   * @param lubrication_force a vector containing the lubrication force on the particles.
   * @param lubrication_torque a vector containing the lubrication torques on the particles.
   */
  void
  calculate_pw_lubrication_force(const double               dt_dem,
                                 const double               h_max,
                                 const double               h_min,
                                 const double               mu,
                                 std::vector<Tensor<1, 3>> &lubrication_force,
                                 std::vector<Tensor<1, 3>> &lubrication_torque);

  /**
   * @brief  Updates the boundary cells that are contact candidates for each of the particles.The force is based on the formula from
   *  Microhydrodynamics: Principles and Selected Applications by Kim, Sangtae;
   * Karrila, Seppo J. ISBN 13: 9780750691734
   *
   * @param particles The particles vector containing all the IB particles.
   * @param dof_handler The dof handler of the mesh used for the fluid simulation.
   * @param face_quadrature_formula The face quadrature formula used in the elements.
   * @param mapping The FEM mapping of the face element.
   */

  void
  update_particles_boundary_contact(
    std::vector<IBParticle<dim>> &particles,
    const DoFHandler<dim>        &dof_handler,
    const Quadrature<dim - 1>    &face_quadrature_formula,
    const Mapping<dim>           &mapping);


  std::vector<IBParticle<dim>> dem_particles;

private:
  bool
  is_boundary_excluded(const unsigned int boundary_id) const;

  void
  get_wall_motion(const unsigned int boundary_id,
                  const Point<dim>  &point_on_boundary,
                  Tensor<1, 3>      &wall_velocity,
                  Tensor<1, 3>      &wall_angular_velocity,
                  Point<dim>        &wall_center_of_rotation) const;

  // A struct to store boundary cells' information
  struct BoundaryCellsInfo
  {
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int /*version*/)
    {
      for (unsigned int i = 0; i < dim; ++i)
        {
          ar &normal_vector;
          ar &point_on_boundary;
          ar &boundary_index;
        }
    }

    Tensor<1, dim> normal_vector;
    Point<dim>     point_on_boundary;
    unsigned int   boundary_index;
  };

  // These structs are used to specify the default value of a variable in a map
  // or set. This is used to find the best particle wall contact candidate.
  struct DefaultDBL_MAX
  {
    double value = DBL_MAX;
  };
  struct DefaultUINT_MAX
  {
    int value = UINT_MAX;
  };

  // This enum defines the lowest index of a floating wall in the particle wall
  // contact. This prevents a wall floating wall from shearing the same index as
  // a standard boundary.
  enum lowest_floating_wall_indices
  {
    lowest_floating_wall_indices = 1000000
  };

  std::shared_ptr<Parameters::IBParticles<dim>> parameters;
  std::shared_ptr<Parameters::Lagrangian::FloatingWalls<dim>>
                           floating_walls_parameters;
  std::shared_ptr<Parameters::Lagrangian::BCDEM>
    boundary_conditions_parameters;
  DEMSolverParameters<dim> dem_parameters{};
  MPI_Comm                 mpi_communicator;

  std::vector<std::set<unsigned int>> particles_contact_candidates;

  // Store the previous contact point as initial guess for the next search
  std::map<unsigned int, std::map<unsigned int, Point<dim>>>
    previous_wall_contact_point;
  std::map<unsigned int, std::map<unsigned int, Point<dim>>>
    previous_particle_particle_contact_point;

  // Particles contact history
  std::map<unsigned int, std::map<unsigned int, ContactInfo>> pp_contact_map;
  std::map<unsigned int, std::map<unsigned int, ContactInfo>> pw_contact_map;

  // A vector of vectors of candidate cells for each of the particle.
  std::vector<std::map<unsigned int, BoundaryCellsInfo>> boundary_cells;

  double cfd_time;
};


#endif // LETHE_IB_PARTICLES_DEM_H
