#include <dem/set_particle_particle_contact_force_model.h>

using namespace dealii;

template <int dim>
std::shared_ptr<ParticleParticleContactForceBase<dim>>
set_particle_particle_contact_force_model(
  const DEMSolverParameters<dim> &dem_parameters)
{
  std::shared_ptr<ParticleParticleContactForceBase<dim>>
    particle_particle_contact_force_object;

  if (dem_parameters.model_parameters.particle_particle_contact_force_model ==
      Parameters::Lagrangian::ParticleParticleContactForceModel::linear)
    {
      switch (dem_parameters.model_parameters.rolling_resistance_method)
        {
          case Parameters::Lagrangian::RollingResistanceMethod::no_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    linear,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    no_resistance>>(dem_parameters);
              break;
            }
          case Parameters::Lagrangian::RollingResistanceMethod::
            constant_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    linear,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    constant_resistance>>(dem_parameters);
              break;
            }
          case Parameters::Lagrangian::RollingResistanceMethod::
            viscous_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    linear,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    viscous_resistance>>(dem_parameters);
              break;
            }
          default:
            throw std::runtime_error(
              "Invalid contact force model and rolling resistance method combination");
        }
    }
  else if (dem_parameters.model_parameters
             .particle_particle_contact_force_model ==
           Parameters::Lagrangian::ParticleParticleContactForceModel::
             hertz_mindlin_limit_overlap)
    {
      switch (dem_parameters.model_parameters.rolling_resistance_method)
        {
          case Parameters::Lagrangian::RollingResistanceMethod::no_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz_mindlin_limit_overlap,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    no_resistance>>(dem_parameters);
              break;
            }
          case Parameters::Lagrangian::RollingResistanceMethod::
            constant_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz_mindlin_limit_overlap,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    constant_resistance>>(dem_parameters);
              break;
            }
          case Parameters::Lagrangian::RollingResistanceMethod::
            viscous_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz_mindlin_limit_overlap,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    viscous_resistance>>(dem_parameters);
              break;
            }
          default:
            throw std::runtime_error(
              "Invalid contact force model and rolling resistance method combination");
        }
    }
  else if (dem_parameters.model_parameters
             .particle_particle_contact_force_model ==
           Parameters::Lagrangian::ParticleParticleContactForceModel::
             hertz_mindlin_limit_force)
    {
      switch (dem_parameters.model_parameters.rolling_resistance_method)
        {
          case Parameters::Lagrangian::RollingResistanceMethod::no_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz_mindlin_limit_force,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    no_resistance>>(dem_parameters);
              break;
            }
          case Parameters::Lagrangian::RollingResistanceMethod::
            constant_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz_mindlin_limit_force,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    constant_resistance>>(dem_parameters);
              break;
            }
          case Parameters::Lagrangian::RollingResistanceMethod::
            viscous_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz_mindlin_limit_force,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    viscous_resistance>>(dem_parameters);
              break;
            }
          default:
            throw std::runtime_error(
              "Invalid contact force model and rolling resistance method combination");
        }
    }
  else if (dem_parameters.model_parameters
             .particle_particle_contact_force_model ==
           Parameters::Lagrangian::ParticleParticleContactForceModel::hertz)
    {
      switch (dem_parameters.model_parameters.rolling_resistance_method)
        {
          case Parameters::Lagrangian::RollingResistanceMethod::no_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    no_resistance>>(dem_parameters);
              break;
            }
          case Parameters::Lagrangian::RollingResistanceMethod::
            constant_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    constant_resistance>>(dem_parameters);
              break;
            }
          case Parameters::Lagrangian::RollingResistanceMethod::
            viscous_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    viscous_resistance>>(dem_parameters);
              break;
            }
          default:
            throw std::runtime_error(
              "Invalid contact force model and rolling resistance method combination");
        }
    }
  else if (dem_parameters.model_parameters
             .particle_particle_contact_force_model ==
           Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR)
    {
      switch (dem_parameters.model_parameters.rolling_resistance_method)
        {
          case Parameters::Lagrangian::RollingResistanceMethod::no_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz_JKR,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    no_resistance>>(dem_parameters);
              break;
            }
          case Parameters::Lagrangian::RollingResistanceMethod::
            constant_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz_JKR,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    constant_resistance>>(dem_parameters);
              break;
            }
          case Parameters::Lagrangian::RollingResistanceMethod::
            viscous_resistance:
            {
              particle_particle_contact_force_object =
                std::make_shared<ParticleParticleContactForce<
                  dim,
                  Parameters::Lagrangian::ParticleParticleContactForceModel::
                    hertz_JKR,
                  Parameters::Lagrangian::RollingResistanceMethod::
                    viscous_resistance>>(dem_parameters);
              break;
            }
          default:
            throw std::runtime_error(
              "Invalid contact force model and rolling resistance method combination");
        }
    }
  else
    {
      throw std::runtime_error(
        "The chosen particle-particle contact force model is invalid");
    }
  return particle_particle_contact_force_object;
}

template std::shared_ptr<ParticleParticleContactForceBase<2>>
set_particle_particle_contact_force_model(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<ParticleParticleContactForceBase<3>>
set_particle_particle_contact_force_model(
  const DEMSolverParameters<3> &dem_parameters);
