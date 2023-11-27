#include <dem/distributions.h>

NormalDistribution::NormalDistribution(const double d_averages,
                                       const double d_standard_deviations)
{
  diameter_averages   = d_averages;
  standard_deviations = d_standard_deviations;
}

void
NormalDistribution::particle_size_sampling(const unsigned int particle_number)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  std::random_device         rd{};
  std::mt19937               gen{rd()};
  std::normal_distribution<> distribution{diameter_averages,
                                          standard_deviations};

  for (unsigned int n = 0; n < particle_number; ++n)
    this->particle_sizes.push_back(distribution(gen));
}

double
NormalDistribution::find_min_diameter()
{
  double min_particle_size =
    diameter_averages -
    2.5 * standard_deviations; // 2.5 -> approx 99% of all diameters are smaller

  return min_particle_size;
}

double
NormalDistribution::find_max_diameter()
{
  double max_particle_size =
    diameter_averages +
    2.5 * standard_deviations; // 2.5 -> approx 99% of all diameters are bigger

  return max_particle_size;
}

UniformDistribution::UniformDistribution(const double d_values)
{
  diameter_values = d_values;
}

void
UniformDistribution::particle_size_sampling(const unsigned int particle_number)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  for (unsigned int n = 0; n < particle_number; ++n)
    this->particle_sizes.push_back(this->diameter_values);
}

double
UniformDistribution::find_min_diameter()
{
  return this->diameter_values;
}

double
UniformDistribution::find_max_diameter()
{
  return this->diameter_values;
}
