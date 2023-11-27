#include <dem/find_maximum_particle_size.h>

double
find_maximum_particle_size(
  const Parameters::Lagrangian::LagrangianPhysicalProperties
    &physical_properties,
  const std::unordered_map<unsigned int, std::shared_ptr<Distribution>>
    &distribution_object_container)
{
  double maximum_particle_size = 0;

  for (unsigned int i = 0; i < physical_properties.particle_type_number; ++i)
    {
      maximum_particle_size =
        std::max(maximum_particle_size,
                 distribution_object_container.at(i)->find_max_diameter());
    }
  return maximum_particle_size;
}
