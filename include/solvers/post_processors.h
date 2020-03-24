
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
class vorticity_postprocessor : public DataPostprocessorVector<dim>
{
public:
  vorticity_postprocessor()
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

// Calculate the Q criterion within each element
template <int dim>
class qcriterion_postprocessor : public DataPostprocessorScalar<dim>
{
public:
  qcriterion_postprocessor()
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


// Calculate the velocity field in the eulerian frame of reference
template <int dim>
class SRF_postprocessor : public DataPostprocessorVector<dim>
{
public:
  SRF_postprocessor(double p_omega_x, double p_omega_y, double p_omega_z)
    : DataPostprocessorVector<dim>(std::string("velocity_eulerian"), update_values | update_quadrature_points),
      omega_x(p_omega_x),
      omega_y(p_omega_y),
      omega_z(p_omega_z)
  {}
  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &inputs,
                        std::vector<Vector<double>> &computed_quantities) const
  {
    AssertDimension(input_data.solution_gradients.size(),
                    computed_quantities.size());

    Tensor<1,dim> omega;
    omega[0]=omega_x;
    omega[1]=omega_y;
    if (dim==3) omega[2]=omega_z;

    const unsigned int n_quadrature_points = inputs.solution_values.size();
    Assert(computed_quantities.size() == n_quadrature_points,
           ExcInternalError());
    Assert(inputs.solution_values[0].size() == dim, ExcInternalError());
    for (unsigned int q = 0; q < n_quadrature_points; ++q)
      {
        Tensor<1,dim> velocity;
        velocity[0] = inputs.solution_values[q][0];
        velocity[1] = inputs.solution_values[q][1];
        if (dim==2)
            velocity = velocity + omega_z * (-1.) * cross_product_2d(inputs.evaluation_points[q]);

        if (dim==3)
          {
            velocity[2] = inputs.solution_values[q][2];
            velocity = velocity + cross_product_3d(omega,inputs.evaluation_points[q]);
          }
        for (int d = 0 ; d < dim ; ++d)
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
