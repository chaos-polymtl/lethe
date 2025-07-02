#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/logstream.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <iomanip>

#include <core/sdirk_stage_data.h>

using namespace dealii;

SDIRKStageData
sdirk_stage_data(const FullMatrix<double> &butcher_table,
                     const std::vector<double> &c,
                     const std::vector<double> &b,
                     const unsigned int         stage_i)
{
  const unsigned int n_stages = b.size();

  if (stage_i >= n_stages)
    throw std::invalid_argument("Invalid stage index.");

  SDIRKStageData data;
  data.a_ij.resize(stage_i + 1);

  for (unsigned int j = 0; j <= stage_i; ++j)
    data.a_ij[j] = butcher_table(stage_i, j);

  data.c_i = c[stage_i];
  data.b_i = b[stage_i];

  return data;
}

