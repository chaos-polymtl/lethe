// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_specific_heat_model_h
#define lethe_specific_heat_model_h

#include <core/parameters.h>
#include <core/phase_change.h>
#include <core/physical_property_model.h>

using namespace dealii;

/**
 * @brief Abstract class for calculating the specific heat of materials on
 * quadrature points using the temperature of the fluid.
 */
class SpecificHeatModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a SpecificHeatModel object
   * by casting it to the proper child class.
   *
   * @param[in] material_properties Property parameters of a material (fluid or
   * solid).
   *
   * @return Casted SpecificHeatModel object.
   */
  static std::shared_ptr<SpecificHeatModel>
  model_cast(const Parameters::Material &material_properties);
};


/**
 * @brief Constant specific heat model.
 */
class ConstantSpecificHeat : public SpecificHeatModel
{
public:
  /**
   * @brief Constructor of the constant specific heat model.
   *
   * @param[in] p_specific_heat Constant specific heat value.
   */
  ConstantSpecificHeat(const double p_specific_heat)
    : specific_heat(p_specific_heat)
  {}

  /**
   * @brief Compute specific heat value.
   *
   * @param[in] field_values Value of the various field on which the specific
   * heat may depend. The map stores a single value per field.
   *
   * @note Here, the @p field_values parameter is ignored since the specific
   * heat remains constant.
   *
   * @return Specific heat value.
   */
  double
  value(const std::map<field, double> &field_values) override
  {
    (void)field_values;
    return specific_heat;
  };

  /**
   * @brief vector_value Compute specific heat values.
   *
   * @param[in] field_vectors Vectors of the fields on which the specific heat
   * may depend. The map stores a vector of values per field.
   *
   * @param[out] property_vector Vector of computed specific heat values.
   *
   * @note Here, the @p field_vectors parameter is ignored since the specific
   * heat remains constant.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    (void)field_vectors;
    std::fill(property_vector.begin(), property_vector.end(), specific_heat);
  }

  /**
   * @brief Compute the jacobian (the partial derivative) of the specific heat
   * with respect to a specified field.
   *
   * @param[in] field_values Values of the various fields on which the specific
   * heat may depend. The map stores a single value per field.
   *
   * @param[in] id Indicator of the field with respect to which the jacobian
   * should be computed.
   *
   * @note Here, the @p field_values and @p id parameters are ignored since the
   * specific heat remains constant.
   *
   * @return Value of the derivative of the specific heat with respect to the
   * specified field. Since the specific heat remains constant, this function
   * returns zero.
   */
  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    (void)field_values;
    (void)id;
    return 0;
  };

  /**
   * @brief Computes the derivative of the specific heat with respect to
   * specified fields.
   *
   * @param[in] field_vectors Vector of values of the fields used to evaluate
   * the specific heat. The map stores a vector of values per field.
   *
   * @param[in] id Identifier of the field with respect to which a derivative
   * should be computed.
   *
   * @param[out] jacobian Vector of computed derivative values of the specific
   * heat
   * with respect to the field of the specified @p id. In this case, it returns
   * a vector of zeros since the specific heat remains constant.
   *
   * @note Here, the @p field_vectors and @p id parameters are ignored since the
   * specific heat remains constant.
   *
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    (void)field_vectors;
    (void)id;
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  };

private:
  const double specific_heat;
};


/**
 * @brief Phase change specific heat model.
 *
 * The phase change specific heat \f$(c^{*}_\mathrm{p})\f$ model takes into
 * account the phase change of a material by considering the latent heat into
 * the specific heat over a phase change temperature interval determined by
 * \f$\left[T_\mathrm{s},T_\mathrm{l}\right]\f$ where \f$T_\mathrm{s}\f$ is the
 * solidus temperature and \f$T_\mathrm{l}\f$ is the liquidus temperature.
 *
 * The specific heat is evaluated as:
 * \f$ c^{*}_\mathrm{p}(T)  =
 *    \begin{cases}
 *      C_\mathrm{p,s} & \text{if} \; T<T_\mathrm{s}\\
 *      \frac{C_\mathrm{p,s}+C_\mathrm{p,l}}{2}+\frac{h_\mathrm{l}}
 *        {T_\mathrm{l}-T_\mathrm{s}} & \text{if} \;
 * T\in[T_\mathrm{s},T_\mathrm{l}]\\ C_\mathrm{p,l} & \text{if} \;
 * T>T_\mathrm{l} \end{cases} \f$
 *
 * where \f$C_\mathrm{p,s}\f$ and \f$C_\mathrm{p,l}\f$ are the solid and liquid
 * specific heat, respectively, and \f$h_\mathrm{l}\f$ is the latent enthalpy
 * (enthalpy related to the phase change).
 */
class PhaseChangeSpecificHeat : public SpecificHeatModel
{
public:
  /**
   * @brief Constructor of the phase change specific heat model.
   *
   * @param[in] p_phase_change_parameters Set of parameters of the phase change
   * model.
   */
  PhaseChangeSpecificHeat(
    const Parameters::PhaseChange p_phase_change_parameters)
    : param(p_phase_change_parameters)
  {
    this->model_depends_on[field::temperature]    = true;
    this->model_depends_on[field::temperature_p1] = true;
    this->model_depends_on[field::temperature_p2] = true;
  }

  /**
   * @brief Compute the specific heat of the isothermal ideal gas.
   *
   * @param[in] field_values Values of the various fields on which the property
   * may depend. In this case, the specific heat depends on temperature.
   * The map stores a single value per field.
   *
   * @return Value of the specific heat computed with the @p field_values.
   */
  double
  value(const std::map<field, double> &field_values) override
  {
    Assert(field_values.find(field::temperature) != field_values.end(),
           PhysicialPropertyModelFieldUndefined("PhaseChangeSpecificHeat",
                                                "temperature"));
    Assert(field_values.find(field::temperature_p1) != field_values.end(),
           PhysicialPropertyModelFieldUndefined("PhaseChangeSpecificHeat",
                                                "temperature_p1"));
    Assert(field_values.find(field::temperature_p2) != field_values.end(),
           PhysicialPropertyModelFieldUndefined("PhaseChangeSpecificHeat",
                                                "temperature_p2"));
    double temperature    = field_values.at(field::temperature);
    double temperature_p1 = field_values.at(field::temperature_p1);
    double temperature_p2 = field_values.at(field::temperature_p2);


    // Gather information required from the simulation control to have the time
    // history
    std::shared_ptr<SimulationControl> simulation_control =
      get_simulation_control();

    // Time stepping information
    auto method = simulation_control->get_assembly_method();
    // Vector for the BDF coefficients
    const Vector<double> &bdf_coefs =
      simulation_control->get_bdf_coefficients();

    switch (method)
      {
        case Parameters::SimulationControl::TimeSteppingMethod::bdf1:
          {
            // Clipping to ensure a minimal temperature difference
            if (temperature > temperature_p1)
              temperature = std::max(temperature, temperature_p1 + 1e-6);
            else
              temperature_p1 = std::max(temperature_p1, temperature + 1e-6);

            const double enthalpy_current = enthalpy(temperature);
            const double enthalpy_p1      = enthalpy(temperature_p1);
            const double dH =
              bdf_coefs[0] * enthalpy_current + bdf_coefs[1] * enthalpy_p1;
            const double dT =
              bdf_coefs[0] * temperature + bdf_coefs[1] * temperature_p1;
            return dH / dT;
          }

        case Parameters::SimulationControl::TimeSteppingMethod::bdf2:
          {
            // Clipping to ensure a minimal temperature difference
            if (temperature > temperature_p2)
              temperature = std::max(temperature, temperature_p2 + 1e-6);
            else
              temperature_p2 = std::max(temperature_p2, temperature + 1e-6);

            const double enthalpy_current = enthalpy(temperature);
            const double enthalpy_p1      = enthalpy(temperature_p1);
            const double enthalpy_p2      = enthalpy(temperature_p2);
            const double dH               = bdf_coefs[0] * enthalpy_current +
                              bdf_coefs[1] * enthalpy_p1 +
                              bdf_coefs[2] * enthalpy_p2;
            const double dT = bdf_coefs[0] * temperature +
                              bdf_coefs[1] * temperature_p1 +
                              bdf_coefs[2] * temperature_p2;
            return dH / dT;
          }

        default:
          throw(std::runtime_error(
            "Other time integration schemes are not supported by the phase change model"));
      }
  }

  /**
   * @brief Compute a vector of specific heat values for an isothermal ideal gas.
   *
   * @param[in] field_vectors Vectors of the fields on which the specific heat
   * may depend. In this case, the specific heat depends on temperature. The map
   * stores a vector of values per field.
   *
   * @param[out] property_vector Vectors of computed specific heat values.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    Assert(field_vectors.find(field::temperature) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined("PhaseChangeSpecificHeat",
                                                "temperature"));
    Assert(field_vectors.find(field::temperature_p1) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined("PhaseChangeSpecificHeat",
                                                "temperature_p1"));
    Assert(field_vectors.find(field::temperature_p2) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined("PhaseChangeSpecificHeat",
                                                "temperature_p2"));
    const std::vector<double> &temperature_vec =
      field_vectors.at(field::temperature);
    const std::vector<double> &p1_temperature_vec =
      field_vectors.at(field::temperature_p1);
    const std::vector<double> &p2_temperature_vec =
      field_vectors.at(field::temperature_p2);

    const unsigned int n_values = temperature_vec.size();

    Assert(n_values == p2_temperature_vec.size(),
           SizeOfFields(n_values, p2_temperature_vec.size()));
    Assert(n_values == p1_temperature_vec.size(),
           SizeOfFields(n_values, p1_temperature_vec.size()));
    Assert(n_values == property_vector.size(),
           SizeOfFields(n_values, property_vector.size()));

    // Gather information required from the simulation control to have the time
    // history
    std::shared_ptr<SimulationControl> simulation_control =
      get_simulation_control();
    std::vector<double> time_steps_vector =
      simulation_control->get_time_steps_vector();

    auto                  method = simulation_control->get_assembly_method();
    const Vector<double> &bdf_coefs =
      simulation_control->get_bdf_coefficients();

    for (unsigned int i = 0; i < n_values; ++i)
      {
        double temperature    = temperature_vec[i];
        double temperature_p1 = p1_temperature_vec[i];
        double temperature_p2 = p2_temperature_vec[i];

        switch (method)
          {
            case Parameters::SimulationControl::TimeSteppingMethod::bdf1:
              {
                // Clipping to ensure a minimal temperature difference
                if (temperature > temperature_p1)
                  temperature = std::max(temperature, temperature_p1 + 1e-6);
                else
                  temperature_p1 = std::max(temperature_p1, temperature + 1e-6);

                const double enthalpy_current = enthalpy(temperature);
                const double enthalpy_p1      = enthalpy(temperature_p1);
                const double dH =
                  bdf_coefs[0] * enthalpy_current + bdf_coefs[1] * enthalpy_p1;
                const double dT =
                  bdf_coefs[0] * temperature + bdf_coefs[1] * temperature_p1;

                // The division does not need a tolerance to avoid a division by
                // 0 since the temperature is clipped
                property_vector[i] = dH / dT;

                break;
              }

            case Parameters::SimulationControl::TimeSteppingMethod::bdf2:
              {
                // Clipping to ensure a minimal temperature difference
                if (temperature > temperature_p2)
                  temperature = std::max(temperature, temperature_p2 + 1e-6);
                else
                  temperature_p2 = std::max(temperature_p2, temperature + 1e-6);

                const double enthalpy_current = enthalpy(temperature);
                const double enthalpy_p1      = enthalpy(temperature_p1);
                const double enthalpy_p2      = enthalpy(temperature_p2);
                const double dH = bdf_coefs[0] * enthalpy_current +
                                  bdf_coefs[1] * enthalpy_p1 +
                                  bdf_coefs[2] * enthalpy_p2;
                const double dT = bdf_coefs[0] * temperature +
                                  bdf_coefs[1] * temperature_p1 +
                                  bdf_coefs[2] * temperature_p2;
                property_vector[i] = dH / (dT + 1e-8);
                break;
              }

            default:
              throw(std::runtime_error(
                "Other time integration schemes are not supported by the phase change model"));
          }
      }
  }

  /**
   * @brief Compute the jacobian (the partial derivative) of the specific heat
   * with respect to a specified field.
   *
   * @param[in] field_values Values of the various fields on which the specific
   * heat may depend. The map stores a single value per field.
   *
   * @param[in] id Indicator of the field with respect to which the jacobian
   * should be computed.
   *
   * @return Value of the derivative of the specific heat with respect to the
   * specified field.
   */
  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    if (id == field::temperature)
      {
        Assert(field_values.find(field::temperature) != field_values.end(),
               PhysicialPropertyModelFieldUndefined(
                 "EvaporationModelTemperature", "temperature"));
        return numerical_jacobian(field_values, field::temperature);
      }
    else
      return 0;
  };

  /**
   * @brief Compute the derivative of the specific heat with respect to a field
   * for an isothermal ideal gas.
   *
   * @param[in] field_vectors Vector of values of the fields used to evaluate
   * the specific heat. The map stores a vector of values per field.
   *
   * @param[in] id Identifier of the field with respect to which a derivative
   * should be computed.
   *
   * @param[out] jacobian Vector of computed derivative values of the specific
   * heat with respect to the field of the specified @p id.
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    Assert(field_vectors.find(field::temperature) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined("EvaporationModelTemperature",
                                                "temperature"));
    vector_numerical_jacobian(field_vectors, id, jacobian_vector);
  };

  /**
   * @brief Compute the enthalpy of a phase change material for a temperature
   * \f$T\f$.
   *
   * The enthalpy \f$(H)\f$ is defined:
   *
   * - In pure liquid \f$(T>T_\mathrm{l})\f$ as:\n
   * \f$H(T) = c_\mathrm{p,s} T_\mathrm{s} + 0.5(c_\mathrm{p,s}+c_\mathrm{p,l})
   * (T_\mathrm{l}-T_\mathrm{s}) + h_\mathrm{l} + c_\mathrm{p,l}
   * (T-T_\mathrm{l})\f$\n\n
   *
   * - In liquid-solid mixture \f$(T_\mathrm{s} \leq T \leq T_\mathrm{l})\f$
   * as:\n
   *\f$H(T) = c_\mathrm{p,s} T_\mathrm{s} + 0.5(c_\mathrm{p,s}+c_\mathrm{p,l})
   * (T_\mathrm{l}-T_\mathrm{s}) + \alpha_\mathrm{l} h_\mathrm{l}\f$\n
   * where \f$\alpha_\mathrm{l} = \frac{T-T_\mathrm{s}}{T_\mathrm{l}-
   * T_\mathrm{s}}\f$ is the liquid fraction.\n\n
   *
   * - In pure solid \f$(T<T_\mathrm{s})\f$ as:\n
   * \f$H(T) = c_\mathrm{p,s} T\f$
   *
   * @param[in] T temperature at which the enthalpy is computed.
   *
   * @return Value of the enthalpy evaluated at \f$T\f$.
   */

  inline double
  enthalpy(const double T)
  {
    if (T > param.T_liquidus)
      return (param.cp_s * param.T_solidus +
              0.5 * (param.cp_l + param.cp_s) *
                (param.T_liquidus - param.T_solidus) +
              param.latent_enthalpy + (T - param.T_liquidus) * param.cp_l);
    else if (T > param.T_solidus)
      {
        liquid_fraction = calculate_liquid_fraction(T, param);
        return (param.cp_s * param.T_solidus +
                0.5 * (param.cp_l + param.cp_s) * (T - param.T_solidus) +
                param.latent_enthalpy * liquid_fraction);
      }
    else
      return param.cp_s * T;
  }

private:
  const Parameters::PhaseChange param;
  double                        liquid_fraction;
};

#endif
