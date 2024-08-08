#include <core/parameters_lagrangian.h>

#include <dem/particle_particle_contact_force.h>

using namespace DEM;
using namespace Parameters::Lagrangian;

template <int                               dim,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
ParticleParticleContactForce<dim, contact_model, rolling_friction_model>::
  ParticleParticleContactForce(const DEMSolverParameters<dim> &dem_parameters)
  : dmt_cut_off_threshold(dem_parameters.model_parameters.dmt_cut_off_threshold)
{
  set_effective_properties(dem_parameters);
}

template <int                               dim,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
void
ParticleParticleContactForce<dim, contact_model, rolling_friction_model>::
  calculate_particle_particle_contact_force(
    DEMContactManager<dim>    &contact_manager,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force)
{
  auto &local_adjacent_particles = contact_manager.local_adjacent_particles;
  auto &ghost_adjacent_particles = contact_manager.ghost_adjacent_particles;
  auto &local_periodic_adjacent_particles =
    contact_manager.local_periodic_adjacent_particles;
  auto &ghost_periodic_adjacent_particles =
    contact_manager.ghost_periodic_adjacent_particles;
  auto &ghost_local_periodic_adjacent_particles =
    contact_manager.ghost_local_periodic_adjacent_particles;

  // Contact forces calculations of local-local and local-ghost particle
  // pairs are performed in separate loops

  // Set the force_calculation_threshold_distance. This is useful for non-
  // contact cohesive force models such as the DMT force model. This is used in
  // every loop.
  const double force_calculation_threshold_distance =
    set_force_calculation_threshold_distance();

  // Looping over local_adjacent_particles values with iterator
  // adjacent_particles_list
  for (auto &&adjacent_particles_list :
       local_adjacent_particles | boost::adaptors::map_values)
    {
      if (!adjacent_particles_list.empty())
        {
          execute_contact_calculation<ContactType::local_particle_particle>(
              adjacent_particles_list, torque, force, dt);
        }
    }

  // Doing the same calculations for local-ghost particle pairs

  // Looping over ghost_adjacent_particles with iterator
  // adjacent_particles_iterator
  for (auto &&adjacent_particles_list :
       ghost_adjacent_particles | boost::adaptors::map_values)
    {
      if (!adjacent_particles_list.empty())
        {
          execute_contact_calculation<ContactType::ghost_particle_particle>(
              adjacent_particles_list, torque, force, dt);
        }
    }

  for (auto &&periodic_adjacent_particles_list :
       local_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      if (!periodic_adjacent_particles_list.empty())
        {
          execute_contact_calculation<
            ContactType::local_periodic_particle_particle>(
            periodic_adjacent_particles_list, torque, force, dt);
        }
    }

  for (auto &&periodic_adjacent_particles_list :
       ghost_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      if (!periodic_adjacent_particles_list.empty())
        {
          execute_contact_calculation<
            ContactType::ghost_periodic_particle_particle>(
            periodic_adjacent_particles_list, torque, force, dt);
        }
    }

  for (auto &&periodic_adjacent_particles_list :
       ghost_local_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      if (!periodic_adjacent_particles_list.empty())
        {
          execute_contact_calculation<
            ContactType::ghost_local_periodic_particle_particle>(
            periodic_adjacent_particles_list, torque, force, dt);
        }
    }
}

// No resistance
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
