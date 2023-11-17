#include <dem/normal_distribution.h>

NormalDistribution::NormalDistribution(
  const std::unordered_map<unsigned int, double> normal_distribution_average,
  const std::unordered_map<unsigned int, double>
    normal_distribution_standard_deviation)
{
  diameter_averages   = normal_distribution_average;
  standard_deviations = normal_distribution_standard_deviation;
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
