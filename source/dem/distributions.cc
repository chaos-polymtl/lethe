#include <dem/distributions.h>

NormalDistribution::NormalDistribution(const double &d_averages,
                                       const double &d_standard_deviations)
{
  diameter_averages   = d_averages;
  standard_deviations = d_standard_deviations;
}

void
NormalDistribution::particle_size_sampling(const unsigned int &particle_number)
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

UniformDistribution::UniformDistribution(const double &d_values)
{
  diameter_values = d_values;
}

void
UniformDistribution::particle_size_sampling(const unsigned int &particle_number)
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

HistogramDistribution::HistogramDistribution(
  const std::vector<double> &d_list,
  const std::vector<double> &d_probabilities)
{
  diameter_hist_values = d_list;
  std::vector<double> cumulative_probability_vector;
  cumulative_probability_vector.reserve(d_probabilities.size());

  double cumulative_value;
  for (double i_probability : d_probabilities)
    {
      cumulative_value += i_probability;
      cumulative_probability_vector.push_back(cumulative_value);
    }

  diameter_hist_cumm_prob = cumulative_probability_vector;
}

void
HistogramDistribution::particle_size_sampling(
  const unsigned int &particle_number)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  std::random_device               rd{};
  std::mt19937                     gen{rd()};
  std::uniform_real_distribution<> dis(0.0, 1.0);

  for (unsigned int i = 0; i < particle_number; ++i)
    {
      // Search to find the appropriate diameter index
      auto it = std::upper_bound(diameter_hist_cumm_prob.begin(),
                                 diameter_hist_cumm_prob.end(),
                                 dis(gen));

      unsigned int index = std::distance(diameter_hist_cumm_prob.begin(), it);

      this->particle_sizes.push_back(diameter_hist_values[index]);
    }
}

double
HistogramDistribution::find_min_diameter()
{
  return *std::min_element(diameter_hist_values.begin(),
                           diameter_hist_values.end());
}

double
HistogramDistribution::find_max_diameter()
{
  return *std::max_element(diameter_hist_values.begin(),
                           diameter_hist_values.end());
}
