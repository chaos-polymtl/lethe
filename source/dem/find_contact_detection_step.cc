#include <dem/find_contact_detection_step.h>

using namespace dealii;

template <int dim>
bool
find_particle_contact_detection_step(
  Particles::ParticleHandler<dim> &particle_handler,
  const double                     dt,
  const double                     smallest_contact_search_criterion,
  MPI_Comm &                       mpi_communicator,
  bool                             sorting_in_subdomains_step,
  std::vector<double> &            displacement)
{
  if (sorting_in_subdomains_step)
    for (auto &d : displacement)
      d = 0.;

  double max_displacement       = 0.;
  bool   contact_detection_step = false;

  // Updating displacement
  for (auto &particle : particle_handler)
    {
      auto &particle_properties   = particle.get_properties();
      auto &particle_displacement = displacement[particle.get_local_index()];

      // Finding displacement of each particle during last step
      particle_displacement +=
        dt * sqrt(particle_properties[DEM::PropertiesIndex::v_x] *
                    particle_properties[DEM::PropertiesIndex::v_x] +
                  particle_properties[DEM::PropertiesIndex::v_y] *
                    particle_properties[DEM::PropertiesIndex::v_y] +
                  particle_properties[DEM::PropertiesIndex::v_z] *
                    particle_properties[DEM::PropertiesIndex::v_z]);

      // Updating maximum displacement of particles
      max_displacement = std::max(max_displacement, particle_displacement);
    }

  // If the maximum displacement of particles exceeds criterion, this step
  // is a contact detection step
  if (max_displacement > smallest_contact_search_criterion)
    {
      contact_detection_step = true;
    }

  // Broadcasting contact detection step value to other processors
  contact_detection_step =
    Utilities::MPI::logical_or(contact_detection_step, mpi_communicator);

  return contact_detection_step;
}

template bool find_particle_contact_detection_step(
  Particles::ParticleHandler<2> &particle_handler,
  const double                   dt,
  const double                   smallest_contact_search_criterion,
  MPI_Comm &                     mpi_communicator,
  bool                           sorting_in_subdomains_step,
  std::vector<double> &          displacement);

template bool find_particle_contact_detection_step(
  Particles::ParticleHandler<3> &particle_handler,
  const double                   dt,
  const double                   smallest_contact_search_criterion,
  MPI_Comm &                     mpi_communicator,
  bool                           sorting_in_subdomains_step,
  std::vector<double> &          displacement);


template <int dim>
bool
find_floating_mesh_mapping_step(
  const double smallest_contact_search_criterion,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> solids)
{
  bool floating_mesh_requires_map = false;

  for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
    {
      double displacement =
        solids[i_solid]->get_max_displacement_since_mapped();
      floating_mesh_requires_map =
        floating_mesh_requires_map ||
        displacement > smallest_contact_search_criterion;
    }

  return floating_mesh_requires_map;
}

template bool
find_floating_mesh_mapping_step(
  const double smallest_contact_search_criterion,
  std::vector<std::shared_ptr<SerialSolid<2, 3>>> solids);

template bool
find_floating_mesh_mapping_step(
  const double smallest_contact_search_criterion,
  std::vector<std::shared_ptr<SerialSolid<1, 2>>> solids);
