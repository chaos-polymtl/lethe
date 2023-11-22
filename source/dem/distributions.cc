#include <dem/distributions.h>

NormalDistribution::NormalDistribution(
  const std::unordered_map<unsigned int, double> d_averages,
  const std::unordered_map<unsigned int, double> d_standard_deviations)
{
  diameter_averages   = d_averages;
  standard_deviations = d_standard_deviations;
}

void
NormalDistribution::particle_size_sampling(const unsigned int particle_number,
                                           const unsigned int particle_type)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  std::random_device         rd{};
  std::mt19937               gen{rd()};
  std::normal_distribution<> distribution{
    diameter_averages.at(particle_type), standard_deviations.at(particle_type)};

  for (unsigned int n = 0; n < particle_number; ++n)
    this->particle_sizes.push_back(distribution(gen));
}


UniformDistribution::UniformDistribution(
  const std::unordered_map<unsigned int, double> d_values)
{
  diameter_values = d_values;
}

void
UniformDistribution::particle_size_sampling(const unsigned int particle_number,
                                            const unsigned int particle_type)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  for (unsigned int n = 0; n < particle_number; ++n)
    this->particle_sizes.push_back(this->diameter_values.at(particle_type));
}
