#include <core/utilities.h>

#include <dem/distributions.h>

#include <algorithm>

NormalDistribution::NormalDistribution(const double &d_average,
                                       const double &d_standard_deviation)
  : diameter_average(d_average)
  , standard_deviation(d_standard_deviation)
{}

void
NormalDistribution::particle_size_sampling(const unsigned int &particle_number)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  std::random_device         rd{};
  std::mt19937               gen{rd()};
  std::normal_distribution<> distribution{diameter_average, standard_deviation};

  for (unsigned int n = 0; n < particle_number; ++n)
    this->particle_sizes.push_back(distribution(gen));
}

double
NormalDistribution::find_min_diameter()
{
  double min_particle_size =
    diameter_average -
    2.5 * standard_deviation; // 2.5 -> approx 99% of all diameters are smaller

  return min_particle_size;
}

double
NormalDistribution::find_max_diameter()
{
  double max_particle_size =
    diameter_average +
    2.5 * standard_deviation; // 2.5 -> approx 99% of all diameters are bigger

  return max_particle_size;
}

UniformDistribution::UniformDistribution(const double &d_values)
  : diameter_value(d_values)
{}

void
UniformDistribution::particle_size_sampling(const unsigned int &particle_number)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  for (unsigned int n = 0; n < particle_number; ++n)
    this->particle_sizes.push_back(this->diameter_value);
}

double
UniformDistribution::find_min_diameter()
{
  return this->diameter_value;
}

double
UniformDistribution::find_max_diameter()
{
  return this->diameter_value;
}

CustomDistribution::CustomDistribution(
  const std::vector<double> &d_list,
  const std::vector<double> &d_probabilities)
  : diameter_custom_values(d_list)
{
  std::vector<double> cumulative_probability_vector, n_i_vector;
  double              n_tot = 0.;

  cumulative_probability_vector.reserve(d_probabilities.size());
  n_i_vector.reserve(d_probabilities.size());

  for (unsigned int i = 0; i < diameter_custom_values.size(); ++i)
    {
      n_i_vector.push_back(d_probabilities[i] / Utilities::fixed_power<3>(
                                                  diameter_custom_values[i]));
      n_tot += n_i_vector[i];
    }

  double i_probability;
  double cumulative_value = 0.;

  for (unsigned int i = 0; i < diameter_custom_values.size(); ++i)
    {
      i_probability = n_i_vector[i] / n_tot;
      cumulative_value += i_probability;
      cumulative_probability_vector.push_back(cumulative_value);
    }

  diameter_custom_cumm_prob = cumulative_probability_vector;
}

void
CustomDistribution::particle_size_sampling(const unsigned int &particle_number)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  std::random_device               rd{};
  std::mt19937                     gen{rd()};
  std::uniform_real_distribution<> dis(0.0, 1.0);

  for (unsigned int i = 0; i < particle_number; ++i)
    {
      // Search to find the appropriate diameter index
      auto it = std::upper_bound(diameter_custom_cumm_prob.begin(),
                                 diameter_custom_cumm_prob.end(),
                                 dis(gen));

      unsigned int index = std::distance(diameter_custom_cumm_prob.begin(), it);

      this->particle_sizes.push_back(diameter_custom_values[index]);
    }
}

double
CustomDistribution::find_min_diameter()
{
  return *std::min_element(diameter_custom_values.begin(),
                           diameter_custom_values.end());
}

double
CustomDistribution::find_max_diameter()
{
  return *std::max_element(diameter_custom_values.begin(),
                           diameter_custom_values.end());
}
