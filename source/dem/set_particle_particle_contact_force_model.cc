#include <dem/set_particle_particle_contact_force_model.h>

using namespace dealii;

template <int dim>
std::shared_ptr<ParticleParticleContactForce<dim>>
set_particle_particle_contact_force_model(
  const DEMSolverParameters<dim> &dem_parameters)
{
  std::shared_ptr<ParticleParticleContactForce<dim>>
    particle_particle_contact_force_object;

  if (dem_parameters.model_parameters.particle_particle_contact_force_method ==
      Parameters::Lagrangian::ModelParameters::
        ParticleParticleContactForceModel::linear)
    {
      particle_particle_contact_force_object =
        std::make_shared<ParticleParticleLinearForce<dim>>(dem_parameters);
    }
  else if (dem_parameters.model_parameters
             .particle_particle_contact_force_method ==
           Parameters::Lagrangian::ModelParameters::
             ParticleParticleContactForceModel::hertz_mindlin_limit_overlap)
    {
      particle_particle_contact_force_object =
        std::make_shared<ParticleParticleHertzMindlinLimitOverlap<dim>>(
          dem_parameters);
    }
  else if (dem_parameters.model_parameters
             .particle_particle_contact_force_method ==
           Parameters::Lagrangian::ModelParameters::
             ParticleParticleContactForceModel::hertz_mindlin_limit_force)
    {
      particle_particle_contact_force_object =
        std::make_shared<ParticleParticleHertzMindlinLimitForce<dim>>(
          dem_parameters);
    }
  else if (dem_parameters.model_parameters
             .particle_particle_contact_force_method ==
           Parameters::Lagrangian::ModelParameters::
             ParticleParticleContactForceModel::hertz)
    particle_particle_contact_force_object =
      std::make_shared<ParticleParticleHertz<dim>>(dem_parameters);
  else
    {
      throw "The chosen particle-particle contact force model is invalid";
    }
  return particle_particle_contact_force_object;
}

template std::shared_ptr<ParticleParticleContactForce<2>>
set_particle_particle_contact_force_model(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<ParticleParticleContactForce<3>>
set_particle_particle_contact_force_model(
  const DEMSolverParameters<3> &dem_parameters);
