
#ifndef LETHE_POSTPROCESSORS_H
#define LETHE_POSTPROCESSORS_H

// DEALII INCLUDES
// Base
#include <deal.II/base/tensor.h>

// Lac
#include <deal.II/lac/vector.h>

// Numerics
#include <deal.II/numerics/data_postprocessor.h>

// STD includes
#include <vector>

using namespace dealii;

// Calculates the vorticity within each element using the velocity vector
template <int dim>
class vorticity_postprocessor : public DataPostprocessorVector<dim> {
public:
  vorticity_postprocessor()
      : DataPostprocessorVector<dim>("vorticity", update_gradients) {}
  virtual void evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &input_data,
      std::vector<Vector<double>> &computed_quantities) const {
    AssertDimension(input_data.solution_gradients.size(),
                    computed_quantities.size());
    for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p) {
      AssertDimension(computed_quantities[p].size(), dim);
      if (dim == 3) {
        computed_quantities[p](0) = (input_data.solution_gradients[p][2][1] -
                                     input_data.solution_gradients[p][1][2]);
        computed_quantities[p](1) = (input_data.solution_gradients[p][0][2] -
                                     input_data.solution_gradients[p][2][0]);
        computed_quantities[p](2) = (input_data.solution_gradients[p][1][0] -
                                     input_data.solution_gradients[p][0][1]);
      } else {
        computed_quantities[p][0] = (input_data.solution_gradients[p][1][0] -
                                     input_data.solution_gradients[p][0][1]);
      }
    }
  }
};

// Calculate the Q criterion within each element
template <int dim>
class qcriterion_postprocessor : public DataPostprocessorScalar<dim> {
public:
  qcriterion_postprocessor()
      : DataPostprocessorScalar<dim>("q_criterion", update_gradients) {}
  virtual void evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &input_data,
      std::vector<Vector<double>> &computed_quantities) const {
    for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p) {
      AssertDimension(computed_quantities[p].size(), 1);
      double p1 = 0.0, r1 = 0.0;
      std::vector<Tensor<2, dim>> vorticity_vector(
          input_data.solution_gradients.size());
      std::vector<Tensor<2, dim>> strain_rate_tensor(
          input_data.solution_gradients.size());
      for (unsigned int j = 0; j < dim; j++) {
        for (unsigned int k = 0; k < dim; k++) {
          vorticity_vector[p][j][k] =
              0.5 * ((input_data.solution_gradients[p][j][k]) -
                     input_data.solution_gradients[p][k][j]);
          strain_rate_tensor[p][j][k] =
              0.5 * ((input_data.solution_gradients[p][j][k]) +
                     input_data.solution_gradients[p][k][j]);
        }
      }
      for (unsigned int m = 0; m < dim; m++) {
        for (unsigned int n = 0; n < dim; n++) {
          p1 += (vorticity_vector[p][m][n]) * (vorticity_vector[p][m][n]);
          r1 += (strain_rate_tensor[p][m][n]) * (strain_rate_tensor[p][m][n]);
        }
      }
      computed_quantities[p] = 0.5 * (p1 - r1);
    }
  }
};

#endif
