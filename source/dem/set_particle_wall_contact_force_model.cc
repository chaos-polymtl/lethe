#include <dem/set_particle_wall_contact_force_model.h>

using namespace dealii;

template <int dim>
std::shared_ptr<ParticleWallContactForce<dim>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<dim>                  &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation)
{
  std::shared_ptr<ParticleWallContactForce<dim>>
    particle_wall_contact_force_object;

  std::vector<types::boundary_id> boundary_index =
    triangulation.get_boundary_ids();
  if (dem_parameters.model_parameters.particle_wall_contact_force_method ==
      Parameters::Lagrangian::ModelParameters::ParticleWallContactForceModel::
        linear)
    {
      particle_wall_contact_force_object =
        std::make_shared<ParticleWallLinearForce<dim>>(dem_parameters,
                                                       boundary_index);
    }
  else if (dem_parameters.model_parameters.particle_wall_contact_force_method ==
           Parameters::Lagrangian::ModelParameters::
             ParticleWallContactForceModel::nonlinear)
    {
      particle_wall_contact_force_object =
        std::make_shared<ParticleWallNonLinearForce<dim>>(dem_parameters,
                                                          boundary_index);
    }
  else
    {
      throw "The chosen particle-wall contact force model is invalid";
    }
  return particle_wall_contact_force_object;
}

template std::shared_ptr<ParticleWallContactForce<2>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<2>                  &dem_parameters,
  const parallel::distributed::Triangulation<2> &triangulation);

template std::shared_ptr<ParticleWallContactForce<3>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<3>                  &dem_parameters,
  const parallel::distributed::Triangulation<3> &triangulation);
