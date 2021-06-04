
#ifndef lethe_post_processors_h
#define lethe_post_processors_h

// DEALII INCLUDES
// Base
#include <deal.II/base/tensor.h>

// Numerics
#include <deal.II/numerics/data_postprocessor.h>

// standard library includes includes
#include <vector>

using namespace dealii;

/**
 * @brief VorticityPostprocessor Post-processor class used to calculate
 * the vorticity field within the domain. The vorticity is defined
 * as $\mathbf{\omega} = \nabla \times \mathbf{u}$ or, in other notation,
 * $\mathbf{\omega} = rot(\mathbf{u})$.
 */
template <int dim>
class VorticityPostprocessor : public DataPostprocessorVector<dim>
{
public:
  VorticityPostprocessor()
    : DataPostprocessorVector<dim>("vorticity", update_gradients)
  {}
  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                        std::vector<Vector<double>> &computed_quantities) const
  {
    AssertDimension(input_data.solution_gradients.size(),
                    computed_quantities.size());
    for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {
        AssertDimension(computed_quantities[p].size(), dim);
        if (dim == 3)
          {
            computed_quantities[p](0) =
              (input_data.solution_gradients[p][2][1] -
               input_data.solution_gradients[p][1][2]);
            computed_quantities[p](1) =
              (input_data.solution_gradients[p][0][2] -
               input_data.solution_gradients[p][2][0]);
            computed_quantities[p](2) =
              (input_data.solution_gradients[p][1][0] -
               input_data.solution_gradients[p][0][1]);
          }
        else
          {
            computed_quantities[p][0] =
              (input_data.solution_gradients[p][1][0] -
               input_data.solution_gradients[p][0][1]);
          }
      }
  }
};

/**
 * @brief QCriterionPostprocessor Post-processor class used to calculate
 * the QCriterion field within the domain. The Q Criterion is defined as
 * as $Q = ||\Omega||^2 -||S||^2$ where $S$ is the symetric deformation tensor
 */
template <int dim>
class QCriterionPostprocessor : public DataPostprocessorScalar<dim>
{
public:
  QCriterionPostprocessor()
    : DataPostprocessorScalar<dim>("q_criterion", update_gradients)
  {}
  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                        std::vector<Vector<double>> &computed_quantities) const
  {
    for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {
        AssertDimension(computed_quantities[p].size(), 1);
        double                      p1 = 0.0, r1 = 0.0;
        std::vector<Tensor<2, dim>> vorticity_vector(
          input_data.solution_gradients.size());
        std::vector<Tensor<2, dim>> strain_rate_tensor(
          input_data.solution_gradients.size());
        for (unsigned int j = 0; j < dim; j++)
          {
            for (unsigned int k = 0; k < dim; k++)
              {
                vorticity_vector[p][j][k] =
                  0.5 * ((input_data.solution_gradients[p][j][k]) -
                         input_data.solution_gradients[p][k][j]);
                strain_rate_tensor[p][j][k] =
                  0.5 * ((input_data.solution_gradients[p][j][k]) +
                         input_data.solution_gradients[p][k][j]);
              }
          }
        for (unsigned int m = 0; m < dim; m++)
          {
            for (unsigned int n = 0; n < dim; n++)
              {
                p1 += (vorticity_vector[p][m][n]) * (vorticity_vector[p][m][n]);
                r1 +=
                  (strain_rate_tensor[p][m][n]) * (strain_rate_tensor[p][m][n]);
              }
          }
        computed_quantities[p] = 0.5 * (p1 - r1);
      }
  }
};


/**
 * @brief Calculates the velocity in the Eulerian frame of reference
 * for simulations which are carried out in a rotating frame.
 * The velocity in the Eulerian frame is $\mathbf{u}_E=\mathbf{u}_L+
 * \mathbf{\Omega} \times \mathbf{R}$ where $\Omega$ is the velocity of the
 * frame of reference and $\mathbf{R}$ the position vector
 */
template <int dim>
class SRFPostprocessor : public DataPostprocessorVector<dim>
{
public:
  SRFPostprocessor(double p_omega_x, double p_omega_y, double p_omega_z)
    : DataPostprocessorVector<dim>(std::string("velocity_eulerian"),
                                   update_values | update_quadrature_points)
    , omega_x(p_omega_x)
    , omega_y(p_omega_y)
    , omega_z(p_omega_z)
  {}
  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &inputs,
                        std::vector<Vector<double>> &computed_quantities) const
  {
    Tensor<1, dim> omega;
    omega[0] = omega_x;
    omega[1] = omega_y;
    if (dim == 3)
      omega[2] = omega_z;

    const unsigned int n_quadrature_points = inputs.solution_values.size();
    for (unsigned int q = 0; q < n_quadrature_points; ++q)
      {
        Tensor<1, dim> velocity;
        velocity[0] = inputs.solution_values[q][0];
        velocity[1] = inputs.solution_values[q][1];
        if (dim == 2)
          velocity = velocity + omega_z * (-1.) *
                                  cross_product_2d(inputs.evaluation_points[q]);

        if (dim == 3)
          {
            velocity[2] = inputs.solution_values[q][2];
            velocity =
              velocity + cross_product_3d(omega, inputs.evaluation_points[q]);
          }
        for (int d = 0; d < dim; ++d)
          {
            computed_quantities[q][d] = velocity[d];
          }
      }
  }

private:
  double omega_x;
  double omega_y;
  double omega_z;
};

#endif
