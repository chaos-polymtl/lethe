// SPDX-FileCopyrightText: Copyright (c) 2019-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_postprocessors_h
#define lethe_postprocessors_h

#include <core/rheological_model.h>

#include <deal.II/base/tensor.h>

#include <deal.II/numerics/data_postprocessor.h>

#include <vector>

using namespace dealii;

/**
 * @brief Compute the vorticity vector field within a domain.
 *
 * The vorticity is defined as \f$\mathbf{\omega} = \nabla \times \mathbf{u}
 * = \mathrm{rot}(\mathbf{u})\f$ where \f$\mathbf{u}\f$ is the velocity vector.
 */
template <int dim>
class VorticityPostprocessor : public DataPostprocessorVector<dim>
{
public:
  /**
   * @brief Constructor of the vorticity postprocessor.
   */
  VorticityPostprocessor()
    : DataPostprocessorVector<dim>("vorticity", update_gradients)
  {}

  /**
   * @brief Compute the vorticity vector from the velocity vector at given
   * evaluation points.
   *
   * @param[in] input_data Fluid dynamics velocity and pressure information.
   *
   * @param[out] computed_quantities Computed vorticity vector at evaluation
   * points.
   *
   * @note @p input_data contains both velocity and pressure data since fluid
   * dynamics in Lethe is solved in a monolithic way, but only the velocity
   * gradient is used to compute the vorticity.
   */
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &input_data,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension(input_data.solution_gradients.size(),
                    computed_quantities.size());
    for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {
        AssertDimension(computed_quantities[p].size(), dim);
        if constexpr (dim == 3)
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
        if constexpr (dim == 2)
          {
            computed_quantities[p][0] =
              (input_data.solution_gradients[p][1][0] -
               input_data.solution_gradients[p][0][1]);
          }
      }
  }
};

/**
 * @brief Compute the Q-criterion scalar field within a domain.
 *
 * The Q-criterion is defined as \f$Q = 0.5 (||\mathbf{\omega}||^2 -
 * ||\mathbf{S}||^2)\f$ where \f$\mathbf{\omega}\f$ is the vorticity vector and
 * \f$\mathbf{S}\f$ is the symmetric deformation tensor.
 */
template <int dim>
class QCriterionPostprocessor : public DataPostprocessorScalar<dim>
{
public:
  /**
   * @brief Constructor of the Q-criterion postprocessor.
   */
  QCriterionPostprocessor()
    : DataPostprocessorScalar<dim>("q_criterion", update_gradients)
  {}

  /**
   * @brief Computes the Q-criterion from the velocity vector at given
   * evaluation points.
   *
   * @param[in] input_data Fluid dynamics velocity and pressure information.
   *
   * @param[out] computed_quantities Computed Q-criterion value at evaluation
   * points.
   *
   * @note @p input_data contains both velocity and pressure data since fluid
   * dynamics in Lethe is solved in a monolithic way, but only the velocity
   * gradient is used to compute the Q-criterion.
   */
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &input_data,
    std::vector<Vector<double>> &computed_quantities) const override
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
 * @brief Compute the velocity divergence scalar field within a domain.
 *
 * The divergence is defined as \f$\nabla \cdot \mathbf{u}\f$ where
 * \f$\mathbf{u}\f$ is the velocity vector.
 */
template <int dim>
class DivergencePostprocessor : public DataPostprocessorScalar<dim>
{
public:
  /**
   * @brief Constructor of the velocity divergence postprocessor.
   */
  DivergencePostprocessor()
    : DataPostprocessorScalar<dim>("velocity_divergence", update_gradients)
  {}

  /**
   * @brief Compute the divergence of the velocity at given evaluation points.
   *
   * @param[in] input_data Fluid dynamics velocity and pressure information.
   *
   * @param[out] computed_quantities Computed velocity divergence value at
   * evaluation points.
   *
   * @note @p input_data contains both velocity and pressure data since fluid
   * dynamics in Lethe is solved in a monolithic way, but only the velocity
   * gradient is used to compute the velocity divergence.
   */
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &input_data,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension(input_data.solution_gradients.size(),
                    computed_quantities.size());
    for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {
        double velocity_divergence = 0;
        for (int d = 0; d < dim; ++d)
          {
            velocity_divergence += input_data.solution_gradients[p][d][d];
          }
        computed_quantities[p] = velocity_divergence;
      }
  }
};


/**
 * @brief Compute the velocity field in the Eulerian frame of reference for
 * simulations carried out in a single rotating frame (SRF).
 *
 * The velocity in the Eulerian frame is \f$\mathbf{u}_\mathrm{E}
 * =\mathbf{u}_\mathrm{L}+ \mathbf{\Omega} \times \mathbf{R}\f$ where
 * \f$\mathbf{\Omega}\f$ is the velocity of the frame of reference and
 * \f$\mathbf{R}\f$ the position vector.
 */
template <int dim>
class SRFPostprocessor : public DataPostprocessorVector<dim>
{
public:
  /**
   * @brief Constructor of the single rotating frame SRF postprocessor.
   *
   * @param[in] p_omega_x Value of the \f$x\f$ component of the SRF's velocity.
   *
   * @param[in] p_omega_y Value of the \f$y\f$ component of the SRF's velocity.
   *
   * @param[in] p_omega_z Value of the \f$z\f$ component of the SRF's velocity.
   */
  SRFPostprocessor(double p_omega_x, double p_omega_y, double p_omega_z)
    : DataPostprocessorVector<dim>(std::string("velocity_eulerian"),
                                   update_values | update_quadrature_points)
    , omega_x(p_omega_x)
    , omega_y(p_omega_y)
    , omega_z(p_omega_z)
  {}

  /**
   * @brief Compute the Eulerian velocity in a SRF at given evaluation points.
   *
   * @param[in] inputs Fluid dynamics velocity and pressure information.
   *
   * @param[out] computed_quantities Computed Eulerian velocity vector at
   * evaluation points.
   *
   * @note @p input_data contains both velocity and pressure data since fluid
   * dynamics in Lethe is solved in a monolithic way, but only the velocity
   * vector is used to compute the Eulerian velocity.
   */
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &computed_quantities) const override
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

        if constexpr (dim == 3)
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
  /// Value of the \f$x\f$ component of the SRF's velocity.
  double omega_x;
  /// Value of the \f$y\f$ component of the SRF's velocity.
  double omega_y;
  /// Value of the \f$z\f$ component of the SRF's velocity.
  double omega_z;
};


/**
 * @brief Compute the kinematic viscosity scalar field of a fluid within a
 * domain.
 *
 * The expression for the kinematic viscosity depends on the rheological
 * model of the fluid. See <a href="https://chaos-polymtl.github.io/lethe/
 * documentation/parameters/cfd/physical_properties.html#rheological-models"
 * target="_blank">documentation on rheological models</a>.
 */
template <int dim>
class KinematicViscosityPostprocessor : public DataPostprocessorScalar<dim>
{
public:
  /**
   * @brief Constructor of the kinematic viscosity postprocessor.
   *
   * @param[in] rheological_model Rheological model of the fluid.
   *
   * @param[in] fluid_id Identifier corresponding to the fluid.
   */
  KinematicViscosityPostprocessor(
    std::shared_ptr<RheologicalModel> rheological_model,
    const unsigned int                fluid_id = 0)
    : DataPostprocessorScalar<dim>("kinematic_viscosity" +
                                     Utilities::to_string(fluid_id, 2),
                                   update_gradients)
    , rheological_model(rheological_model)
  {}

  /**
   * @brief Compute the kinematic viscosity of a fluid at given evaluation
   * points.
   *
   * @param[in] inputs Fluid dynamics velocity and pressure information.
   *
   * @param[out] computed_quantities Computed kinematic viscosity value at
   * evaluation points.
   *
   * @note @p input_data contains both velocity and pressure data since fluid
   * dynamics in Lethe is solved in a monolithic way, but only the velocity
   * gradient is used to compute the kinematic viscosity when the rheological
   * model requires it.
   */
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    const unsigned int n_quadrature_points = inputs.solution_gradients.size();

    for (unsigned int q = 0; q < n_quadrature_points; ++q)
      {
        const auto     gradient = inputs.solution_gradients[q];
        Tensor<2, dim> shear_rate;
        for (int i = 0; i < dim; ++i)
          {
            for (int j = 0; j < dim; ++j)
              {
                shear_rate[i][j] += gradient[i][j] + gradient[j][i];
              }
          }
        const double shear_rate_magnitude =
          calculate_shear_rate_magnitude(shear_rate);

        std::map<field, double> field_values;
        field_values[field::shear_rate] = shear_rate_magnitude;

        computed_quantities[q] = rheological_model->value(field_values);
      }
  }

private:
  /// Model describing the fluid's rheology.
  std::shared_ptr<RheologicalModel> rheological_model;
};


/**
 * @brief Compute the dynamic viscosity scalar field of a fluid within a
 * domain.
 *
 * The dynamic viscosity is defined as \f$\mu = \frac{\nu}{\rho_\mathrm{ref}}\f$
 * where \f$\nu\f$ is the kinematic viscosity and \f$\rho_\mathrm{ref}\f$ is the
 * density of the fluid at reference state.
 *
 * The expression for the kinematic viscosity depends on the rheological
 * model of the fluid. See <a href="https://chaos-polymtl.github.io/lethe/
 * documentation/parameters/cfd/physical_properties.html#rheological-models"
 * target="_blank">documentation on rheological models</a>.
 */
template <int dim>
class DynamicViscosityPostprocessor : public DataPostprocessorScalar<dim>
{
public:
  /**
   * @brief Constructor of the dynamic viscosity postprocessor.
   *
   * @param[in] rheological_model Rheological model of the fluid.
   *
   * @param[in] density_ref Density of the fluid at reference state.
   *
   * @param[in] fluid_id Identifier corresponding to the fluid.
   */
  DynamicViscosityPostprocessor(
    std::shared_ptr<RheologicalModel> rheological_model,
    const double                      density_ref,
    const unsigned int                fluid_id = 0)
    : DataPostprocessorScalar<dim>("dynamic_viscosity_" +
                                     Utilities::to_string(fluid_id, 2),
                                   update_gradients)
    , rheological_model(rheological_model)
    , density_ref(density_ref)
  {}

  /**
   * @brief Compute the dynamic viscosity of a fluid at given evaluation
   * points.
   *
   * @param[in] inputs Fluid dynamics velocity and pressure information.
   *
   * @param[out] computed_quantities Computed dynamic viscosity value at
   * evaluation points.
   *
   * @note @p input_data contains both velocity and pressure data since fluid
   * dynamics in Lethe is solved in a monolithic way, but only the velocity
   * gradient is used to compute the dynamic viscosity when the rheological
   * model requires it.
   */
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    const unsigned int n_quadrature_points = inputs.solution_gradients.size();

    for (unsigned int q = 0; q < n_quadrature_points; ++q)
      {
        const auto     gradient = inputs.solution_gradients[q];
        Tensor<2, dim> shear_rate;
        for (int i = 0; i < dim; ++i)
          {
            for (int j = 0; j < dim; ++j)
              {
                shear_rate[i][j] += gradient[i][j] + gradient[j][i];
              }
          }
        const double shear_rate_magnitude =
          calculate_shear_rate_magnitude(shear_rate);

        std::map<field, double> field_values;
        field_values[field::shear_rate] = shear_rate_magnitude;

        computed_quantities[q] =
          rheological_model->get_dynamic_viscosity(density_ref, field_values);
      }
  }

private:
  /// Model describing the fluid's rheology.
  std::shared_ptr<RheologicalModel> rheological_model;
  /// Density of the fluid at reference state.
  double density_ref;
};


/**
 * @brief Compute the shear rate scalar field within a domain.
 *
 * The shear rate is defined as \f$\dot \gamma  = \frac{\partial u_i}
 * {\partial x_j} + \frac{\partial u_j}{\partial x_i}\f$ where \f$i\f$ and
 * \f$j\f$ are dimension component indexes.
 */
template <int dim>
class ShearRatePostprocessor : public DataPostprocessorScalar<dim>
{
public:
  /**
   * @brief Constructor of the shear rate postprocessor.
   */
  ShearRatePostprocessor()
    : DataPostprocessorScalar<dim>("shear_rate", update_gradients)
  {}

  /**
   * @brief Compute the shear rate at given evaluation points.
   *
   * @param[in] inputs Fluid dynamics velocity and pressure information.
   *
   * @param[out] computed_quantities Computed shear rate value at evaluation
   * points.
   *
   * @note @p input_data contains both velocity and pressure data since fluid
   * dynamics in Lethe is solved in a monolithic way, but only the velocity
   * gradient is used to compute the shear rate.
   */
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    const unsigned int n_quadrature_points = inputs.solution_gradients.size();

    for (unsigned int q = 0; q < n_quadrature_points; ++q)
      {
        const auto     gradient = inputs.solution_gradients[q];
        Tensor<2, dim> shear_rate;
        for (int i = 0; i < dim; ++i)
          {
            for (int j = 0; j < dim; ++j)
              {
                shear_rate[i][j] += gradient[i][j] + gradient[j][i];
              }
          }
        computed_quantities[q] = calculate_shear_rate_magnitude(shear_rate);
      }
  }
};


/**
 * @brief Compute a scalar field for a given function within a domain.
 */
template <int dim>
class ScalarFunctionPostprocessor : public DataPostprocessorScalar<dim>
{
public:
  /**
   * @brief Constructor of the scalar function postprocessor.
   *
   * @param[in] name Name of the outputted field.
   *
   * @param[in] scalar_function Postprocessing scalar function.
   *
   */
  ScalarFunctionPostprocessor(
    std::string                          name,
    const std::shared_ptr<Function<dim>> scalar_function)
    : DataPostprocessorScalar<dim>(name, update_quadrature_points)
    , scalar_function(scalar_function)
  {}

  /**
   * @brief Compute the scalar given by the postprocessing function at given
   * evaluation points.
   *
   * @param[in] input_data Independent variable of the postprocessing function.
   *
   * @param[out] computed_quantities Computed scalar value at evaluation
   * points.
   */
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &input_data,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    const unsigned int n_quadrature_points =
      input_data.evaluation_points.size();

    for (unsigned int q = 0; q < n_quadrature_points; ++q)
      {
        AssertDimension(computed_quantities[q].size(), 1);

        const Point<dim> evaluation_point = input_data.evaluation_points[q];
        computed_quantities[q] = scalar_function->value(evaluation_point);
      }
  }

private:
  /// Scalar postprocessing function
  std::shared_ptr<Function<dim>> scalar_function;
};

/**
 * @brief Compute level set field for a given shape within a domain.
 */
template <int dim>
class LevelsetPostprocessor : public DataPostprocessorScalar<dim>
{
public:
  /**
   * @brief Constuctor of the level set postprocessor.
   *
   * @param[in] shape Geometrical entity.
   */
  LevelsetPostprocessor(const std::shared_ptr<Shape<dim>> shape)
    : DataPostprocessorScalar<dim>("levelset", update_quadrature_points)
    , shape(shape)
  {}

  /**
   * @brief Compute signed distance value of evaluation points with respect to
   * the specified shape.
   *
   * @param[in] input_data Fluid dynamics velocity and pressure information.
   *
   * @param[out] computed_quantities Computed signed distance value at
   * evaluation points.
   *
   * @note @p input_data contains both velocity and pressure data since fluid
   * dynamics in Lethe is solved in a monolithic way, but here only the
   * positions of evaluation points are used.
   */
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &input_data,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    const unsigned int n_quadrature_points =
      input_data.evaluation_points.size();

    for (unsigned int q = 0; q < n_quadrature_points; ++q)
      {
        AssertDimension(computed_quantities[q].size(), 1);

        const Point<dim> evaluation_point = input_data.evaluation_points[q];
        const auto       cell             = input_data.template get_cell<dim>();
        computed_quantities[q] =
          shape->value_with_cell_guess(evaluation_point, cell);
      }
  }

private:
  /// Geometrical entity
  std::shared_ptr<Shape<dim>> shape;
};

/**
 * @brief Compute a level set gradient vector field for a given shape within a
 * domain.
 */
template <int dim>
class LevelsetGradientPostprocessor : public DataPostprocessorVector<dim>
{
public:
  /**
   * @brief Constructor for the level set gradient postprocessor
   *
   * @param[in] shape Geometrical entity.
   */
  LevelsetGradientPostprocessor(const std::shared_ptr<Shape<dim>> shape)
    : DataPostprocessorVector<dim>(std::string("levelset_gradient"),
                                   update_values | update_quadrature_points)
    , shape(shape)
  {}

  /**
   * @brief Compute signed distance gradient of evaluation points with respect
   * to the specified shape.
   *
   * @param[in] input_data Fluid dynamics velocity and pressure information.
   *
   * @param[out] computed_quantities Computed signed distance value at
   * evaluation points.
   *
   * @note @p input_data contains both velocity and pressure data since fluid
   * dynamics in Lethe is solved in a monolithic way, but here only the
   * positions of evaluation points are used.
   */
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    const unsigned int n_quadrature_points = inputs.evaluation_points.size();

    Tensor<1, dim> levelset_gradient{};
    for (unsigned int q = 0; q < n_quadrature_points; ++q)
      {
        AssertDimension(computed_quantities[q].size(), dim);

        const Point<dim> evaluation_point = inputs.evaluation_points[q];
        const auto       cell             = inputs.template get_cell<dim>();
        levelset_gradient =
          shape->gradient_with_cell_guess(evaluation_point, cell);
        for (int d = 0; d < dim; d++)
          computed_quantities[q][d] = levelset_gradient[d];
      }
  }

private:
  /// Geometrical entity
  std::shared_ptr<Shape<dim>> shape;
};


/**
 * @brief Compute the density scalar field of a material within a domain.
 *
 * The expression for the density depends on the density model of the material.
 * See <a
 * href="https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/
 * physical_properties.html#density-models" target="_blank">documentation on
 * density models</a>.
 */
template <int dim>
class DensityPostprocessor : public DataPostprocessorScalar<dim>
{
public:
  /**
   * @brief Constructor for the density postprocessor.
   *
   * @param[in] p_density_model Density model of the material.
   *
   * @param[in] material_id Identifier corresponding to the material (fluid or
   * solid).
   */
  DensityPostprocessor(std::shared_ptr<DensityModel> p_density_model,
                       const unsigned int            material_id = 0)
    : DataPostprocessorScalar<dim>("density_" +
                                     Utilities::to_string(material_id, 2),
                                   update_values)
    , density_model(p_density_model)
  {}

  /**
   * @brief Compute the density of a material at given evaluation points.
   *
   * @param[in] inputs Fluid dynamics velocity and pressure information.
   *
   * @param[out] computed_quantities Computed density value at evaluation
   * points.
   */
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    const unsigned int      n_quadrature_points = inputs.solution_values.size();
    std::map<field, double> field_values;

    for (unsigned int q = 0; q < n_quadrature_points; ++q)
      {
        AssertDimension(computed_quantities[q].size(), 1);

        field_values[field::pressure] = inputs.solution_values[q][dim];
        computed_quantities[q]        = density_model->value(field_values);
      }
  }

private:
  /// Model describing the material's density.
  std::shared_ptr<DensityModel> density_model;
};

/**
 * @brief Compute the local heat flux vector field within a domain.
 *
 * The local heat flux is defined as \f$\mathbf{q} = -k \nabla T\f$ where
 * \f$k\f$ is the thermal conductivity of the material and \f$\nabla T\f$ the
 * temperature gradient.
 */
template <int dim>
class HeatFluxPostprocessor : public DataPostprocessorVector<dim>
{
public:
  /**
   * @brief Constructor of the heat flux postprocessor.
   *
   * @param[in] p_thermal_conductivity_model Thermal conductivity model of the
   * material
   *
   * @param[in] p_material_string Type of material: "f" (fluid) or "s" (solid)
   *
   * @param[in] p_id Either fluid ID or solid ID (as initialized in the physical
   * properties). It is used to name the postprocessed variable vector.
   *
   * @param[in] p_mesh_material_id Identifier corresponding to the material as
   * identified in the mesh. By convention, fluids are given the material_id 0
   * and solids follow with ids 1, 2, etc. This is used to identify which cell
   * lies in which material.
   *
   * @param[in] p_solid_phase_present Boolean indicating if the simulation has
   * solids. If there is no solid, the @p material_id of cells is not checked.
   */
  HeatFluxPostprocessor(const std::shared_ptr<ThermalConductivityModel>
                                           &p_thermal_conductivity_model,
                        const std::string  &p_material_string     = "f",
                        const unsigned int &p_id                  = 0,
                        const unsigned int &p_mesh_material_id    = 0,
                        const bool         &p_solid_phase_present = false)
    : DataPostprocessorVector<dim>("heat_flux_" + p_material_string +
                                     Utilities::to_string(p_id, 2),
                                   update_values | update_gradients)
    , thermal_conductivity_model(p_thermal_conductivity_model)
    , mesh_material_id(p_mesh_material_id)
    , solid_phase_present(p_solid_phase_present)
  {}

  /**
   * @brief Compute local heat fluxes from temperature gradient vector field and
   * thermal conductivity at given evaluation points.
   *
   * @param[in] inputs Temperature information.
   *
   * @param[out] computed_quantities Computed local heat flux vectors at
   * evaluation points.
   */
  virtual void
  evaluate_scalar_field(
    const DataPostprocessorInputs::Scalar<dim> &inputs,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    const unsigned int n_evaluation_points = computed_quantities.size();
    AssertDimension(inputs.solution_gradients.size(), n_evaluation_points);

    //  If solid materials are present, identify each material's subdomain.
    if (solid_phase_present)
      {
        // Get Cell
        const auto cell = inputs.template get_cell<dim>();

        // Check if the cell lies in the material of interest
        if (cell->material_id() == mesh_material_id)
          {
            compute_quantities(inputs, computed_quantities);
          }
        // Apply null values to irrelevant subdomain
        else
          {
            Vector<double>              null_vector(dim);
            std::vector<Vector<double>> null_computed_quantities(
              n_evaluation_points, null_vector);
            computed_quantities = null_computed_quantities;
          }
      }
    // Only fluids are present
    else
      {
        compute_quantities(inputs, computed_quantities);
      }
  }

private:
  /**
   * @brief Computes local heat flux vectors.
   *
   * @note The function is used to improve readability of the
   * HeatFluxPostprocessor::evaluate_scalar_field function and avoid repeated
   * lines.
   *
   * @param[in] inputs Temperature information.
   *
   * @param[out] computed_quantities Computed local heat flux vectors at
   * evaluation points.
   */
  void
  compute_quantities(const DataPostprocessorInputs::Scalar<dim> &inputs,
                     std::vector<Vector<double>> &computed_quantities) const
  {
    std::map<field, double> field_values;
    for (unsigned int p = 0; p < computed_quantities.size(); ++p)
      {
        AssertDimension(computed_quantities[p].size(), dim);
        field_values[field::temperature] = inputs.solution_values[p];

        for (unsigned int i = 0; i < dim; ++i)
          {
            computed_quantities[p][i] =
              -thermal_conductivity_model->value(field_values) *
              inputs.solution_gradients[p][i];
          }
      }
  }

  /// Model describing of thermal conductivity of the material.
  std::shared_ptr<ThermalConductivityModel> thermal_conductivity_model;
  /// Identifier corresponding to the material as identified in the mesh.
  const unsigned int mesh_material_id;
  /// Boolean indicating if the simulation has solids.
  const bool solid_phase_present;
};
#endif