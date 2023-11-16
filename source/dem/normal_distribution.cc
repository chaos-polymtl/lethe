#include <dem/normal_distribution.h>

NormalDistribution::NormalDistribution(
  const double normal_distribution_average,
  const double normal_distribution_standard_deviation)
{
  average            = normal_distribution_average;
  standard_deviation = normal_distribution_standard_deviation;
}

void
NormalDistribution::particle_size_sampling(const unsigned int particle_number)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  std::random_device         rd{};
  std::mt19937               gen{rd()};
  std::normal_distribution<> distribution{average, standard_deviation};

  for (unsigned int n = 0; n < particle_number; ++n)
    this->particle_sizes.push_back(distribution(gen));
}
