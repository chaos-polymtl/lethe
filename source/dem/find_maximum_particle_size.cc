#include <dem/find_maximum_particle_size.h>

using namespace dealii;

template <int dim>
double
find_maximum_particle_size(
  const Parameters::Lagrangian::PhysicalProperties<dim> &physical_properties)
{
  double maximum_particle_size = 0;

  for (unsigned int i = 0; i < physical_properties.particle_type_number; ++i)
    {
      maximum_particle_size =
        std::max(maximum_particle_size,
                 physical_properties.particle_average_diameter.at(i) +
                   standard_deviation_multiplier *
                     physical_properties.particle_size_std.at(i));
    }
  return maximum_particle_size;
}

template double
find_maximum_particle_size(
  const Parameters::Lagrangian::PhysicalProperties<2> &physical_properties);

template double
find_maximum_particle_size(
  const Parameters::Lagrangian::PhysicalProperties<3> &physical_properties);
