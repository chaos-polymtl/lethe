/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */
#ifndef lethe_evaporation_model_h
#define lethe_evaporation_model_h

#include <core/parameters.h>
// only for the field definition
#include <core/physical_property_model.h>

using namespace dealii;

/**
 * @brief EvaporationModel. Abstract class that allows to calculate the
 * evaporation mass, heat and momentum fluxes.
 */
class EvaporationModel
{
public:
  EvaporationModel()
  {
    model_depends_on[temperature] = false;
  }

  /**
   * @brief Instantiates and returns a pointer to a EvaporationModel
   * object by casting it to the proper child class
   *
   * @param evaporation_parameters Evaporation model parameters
   */
  static std::shared_ptr<EvaporationModel>
  model_cast(const Parameters::Evaporation &evaporation_parameters);


  /**
   * @brief Returns true if the EvaporationModel depends on a field, false if not.
   */
  inline bool
  depends_on(const field &id)
  {
    return model_depends_on[id];
  }

  /**
   * @brief mass_flux Calculates the value of the evaporation mass flux.
   * @param fields_value Value of the various field on which the flux may depend.
   * @return value of the flux calculated with the fields_value
   */
  virtual double
  mass_flux(const std::map<field, double> &fields_value) = 0;

  /**
   * @brief vector_mass_flux Calculates the values of the evaporation mass flux
   * for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param mass_flux_vector Vectors of the mass flux values
   */
  virtual void
  vector_mass_flux(const std::map<field, std::vector<double>> &field_vectors,
                   std::vector<double> &mass_flux_vector) = 0;

  /**
   * @brief heat_flux Calculates the value of the evaporation heat flux.
   * @param fields_value Value of the various field on which the flux may depend.
   * @return value of the heat flux calculated with the fields_value
   */
  virtual double
  heat_flux(const std::map<field, double> &fields_value) = 0;

  /**
   * @brief vector_heat_flux Calculates the values of the evaporation heat flux
   * for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param heat_flux_vector Vectors of the heat flux values
   */
  virtual void
  vector_heat_flux(const std::map<field, std::vector<double>> &field_vectors,
                   std::vector<double> &heat_flux_vector) = 0;

  /**
   * @brief heat_flux_jacobian Calculates the jacobian (the partial derivative)
   * of the evaporation heat flux with respect to a field
   * @param field_values Value of the various fields on which the flux may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the heat flux with respect to the field.
   */
  virtual double
  heat_flux_jacobian(const std::map<field, double> &field_values,
                     const field                    id) = 0;

  /**
   * @brief vector_heat_flux_jacobian Calculate the derivative of the evaporation
   * heat flux with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate
   * the flux
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated
   * @param jacobian Vector of the values of the derivative of the heat flux with
   * respect to the field id
   */
  virtual void
  vector_heat_flux_jacobian(
    const std::map<field, std::vector<double>> &field_vectors,
    const field                                 id,
    std::vector<double>                        &jacobian_vector) = 0;

  /**
   * @brief momentum_flux Calculates the value of the evaporation recoil pressure.
   * @param fields_value Value of the various field on which the recoil pressure
   * may depend.
   * @return value of the recoil pressure calculated with the fields_value
   */
  virtual double
  momentum_flux(const std::map<field, double> &fields_value) = 0;

  /**
   * @brief vector_momentum_flux Calculates the values of the evaporation recoil
   * pressure for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param momentum_flux_vector Vectors of the recoil pressure values
   */
  virtual void
  vector_momentum_flux(
    const std::map<field, std::vector<double>> &field_vectors,
    std::vector<double>                        &momentum_flux_vector) = 0;

  /**
   * @brief momentum_flux_jacobian Calculates the jacobian (the partial derivative)
   * of the evaporation recoil pressure with respect to a field
   * @param field_values Value of the various fields on which the recoil pressure
   * may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the recoil pressure with respect
   * to the field.
   */
  virtual double
  momentum_flux_jacobian(const std::map<field, double> &field_values,
                         const field                    id) = 0;

  /**
   * @brief vector_momentum_flux_jacobian Calculate the derivative of the
   * evaporation recoil pressure with respect to a field
   * @param field_vectors Vector for the values of the fields on which the recoil
   * pressure may depend.
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated
   * @param jacobian Vector of the value of the derivative of the recoil pressure
   * with respect to the field id
   */
  virtual void
  vector_momentum_flux_jacobian(
    const std::map<field, std::vector<double>> &field_vectors,
    const field                                 id,
    std::vector<double>                        &jacobian_vector) = 0;

protected:
  // Map to indicate on which variables the model depends on
  std::map<field, bool> model_depends_on;
};

/**
 * @brief EvaporationModelConstant. Class that allows to calculate the
 * evaporation heat and momentum fluxes when a constant mass flux is
 * considered.
 */
class EvaporationModelConstant : public EvaporationModel
{
public:
  EvaporationModelConstant(const Parameters::Evaporation &p_evaporation)
    : n_evaporation(p_evaporation.n_evaporation)
    , latent_heat_evaporation(p_evaporation.latent_heat_evaporation)
    , ambient_gas_density_inv(1.0 / p_evaporation.ambient_gas_density)
    , liquid_density_inv(1.0 / p_evaporation.liquid_density)
  {
    model_depends_on[temperature] = false;
  }


  /**
   * @brief mass_flux Calculates the value of the evaporation mass flux.
   * @param fields_value Value of the various field on which the flux may depend.
   * @return value of the flux calculated with the fields_value
   */
  double
  mass_flux(const std::map<field, double> & /*fields_value*/) override
  {
    return n_evaporation;
  }

  /**
   * @brief vector_mass_flux Calculates the values of the evaporation mass flux
   * for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param mass_flux_vector Vectors of the mass flux values
   */
  void
  vector_mass_flux(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    std::vector<double> &mass_flux_vector) override
  {
    std::fill(mass_flux_vector.begin(), mass_flux_vector.end(), n_evaporation);
  }

  /**
   * @brief heat_flux Calculates the value of the evaporation heat flux.
   * @param fields_value Value of the various field on which the flux may depend.
   * @return value of the heat flux calculated with the fields_value
   */
  double
  heat_flux(const std::map<field, double> & /*fields_value*/) override
  {
    return n_evaporation * latent_heat_evaporation;
  }

  /**
   * @brief vector_heat_flux Calculates the values of the evaporation heat flux
   * for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param heat_flux_vector Vectors of the heat flux values
   */
  void
  vector_heat_flux(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    std::vector<double> &heat_flux_vector) override
  {
    std::fill(heat_flux_vector.begin(),
              heat_flux_vector.end(),
              n_evaporation * latent_heat_evaporation);
  }

  /**
   * @brief heat_flux_jacobian Calculates the jacobian (the partial derivative)
   * of the evaporation heat flux with respect to a field
   * @param field_values Value of the various fields on which the flux may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the heat flux with respect to the field.
   */
  double
  heat_flux_jacobian(const std::map<field, double> & /*fields_value*/,
                     const field /*id*/) override
  {
    return 0.0;
  }

  /**
   * @brief vector_heat_flux_jacobian Calculate the derivative of the evaporation
   * heat flux with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate
   * the flux
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated
   * @param jacobian Vector of the values of the derivative of the heat flux with
   * respect to the field id
   */
  void
  vector_heat_flux_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }
  /**
   * @brief momentum_flux Calculates the value of the evaporation recoil pressure.
   * @param fields_value Value of the various field on which the recoil pressure
   * may depend.
   * @return value of the recoil pressure calculated with the fields_value
   */
  double
  momentum_flux(const std::map<field, double> & /*fields_value*/) override
  {
    return -n_evaporation * n_evaporation *
           (liquid_density_inv - ambient_gas_density_inv);
  }

  /**
   * @brief vector_momentum_flux Calculates the values of the evaporation recoil
   * pressure for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param momentum_flux_vector Vectors of the recoil pressure values
   */
  void
  vector_momentum_flux(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    std::vector<double> &momentum_flux_vector) override
  {
    const double momentum_flux_value =
      -n_evaporation * n_evaporation *
      (liquid_density_inv - ambient_gas_density_inv);

    std::fill(momentum_flux_vector.begin(),
              momentum_flux_vector.end(),
              momentum_flux_value);
  }

  /**
   * @brief momentum_flux_jacobian Calculates the jacobian (the partial derivative)
   * of the evaporation recoil pressure with respect to a field
   * @param field_values Value of the various fields on which the recoil pressure
   * may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the recoil pressure with respect
   * to the field.
   */
  double
  momentum_flux_jacobian(const std::map<field, double> & /*field_values*/,
                         const field /*id*/) override
  {
    return 0.0;
  }

  /**
   * @brief vector_momentum_flux_jacobian Calculate the derivative of the
   * evaporation recoil pressure with respect to a field
   * @param field_vectors Vector for the values of the fields on which the recoil
   * pressure may depend.
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated
   * @param jacobian Vector of the value of the derivative of the recoil pressure
   * with respect to the field id
   */
  void
  vector_momentum_flux_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }

private:
  const double n_evaporation;
  const double latent_heat_evaporation;
  const double ambient_gas_density_inv;
  const double liquid_density_inv;
};

/**
 * @brief EvaporationModelConstant. Class that allows to calculate the
 * evaporation heat and momentum fluxes when a constant mass flux is
 * considered.
 */
class EvaporationModelTemperature : public EvaporationModel
{
public:
  EvaporationModelTemperature(const Parameters::Evaporation &p_evaporation)
    : evaporation_coefficient(p_evaporation.evaporation_coefficient)
    , recoil_pressure_coefficient(p_evaporation.recoil_pressure_coefficient)
    , molar_mass(p_evaporation.molar_mass)
    , boiling_temperature_inv(1.0 / p_evaporation.boiling_temperature)
    , latent_heat_evaporation(p_evaporation.latent_heat_evaporation)
    , ambient_pressure(p_evaporation.ambient_pressure)
    , ambient_gas_density_inv(1.0 / p_evaporation.ambient_gas_density)
    , liquid_density_inv(1.0 / p_evaporation.liquid_density)
  {
    model_depends_on[temperature] = true;
    this->R_inv                   = 1.0 / 8.3145;
    this->L_vap_x_M_x_R_inv =
      latent_heat_evaporation * molar_mass * this->R_inv;
    this->M_x_2PI_inv_x_R_inv = molar_mass / (2.0 * M_PI) * this->R_inv;
  }

  /**
   * @brief saturation_pressure Calculates the value of the staturation pressure.
   * @param fields_value Value of the various field on which the flux may depend.
   * @return value of the flux calculated with the fields_value
   */
  double
  saturation_pressure(const std::map<field, double> &fields_value)
  {
    const double temperature_inv = 1.0 / fields_value.at(field::temperature);

    return ambient_pressure *
           std::exp(-L_vap_x_M_x_R_inv *
                    (temperature_inv - boiling_temperature_inv));
  }

  /**
   * @brief vector_mass_flux Calculates the values of the staturation pressure
   * for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param mass_flux_vector Vectors of the mass flux values
   */
  void
  vector_saturation_pressure(
    const std::map<field, std::vector<double>> &field_vectors,
    std::vector<double>                        &saturation_pressure_vector)
  {
    const std::vector<double> &temperature =
      field_vectors.at(field::temperature);

    for (unsigned int i = 0; i < saturation_pressure_vector.size(); ++i)
      {
        const double temperature_inv = 1.0 / temperature[i];

        saturation_pressure_vector[i] =
          ambient_pressure *
          std::exp(-L_vap_x_M_x_R_inv *
                   (temperature_inv - boiling_temperature_inv));
      }
  }

  /**
   * @brief mass_flux Calculates the value of the saturated vapor pressure.
   * @param fields_value Value of the various field on which the saturated
   * vapor pressure may depend.
   * @return value of the pressure calculated with the fields_value
   */
  double
  mass_flux(const std::map<field, double> &fields_value) override
  {
    const double saturation_pressure_value = saturation_pressure(fields_value);
    const double temperature_inv = 1.0 / fields_value.at(field::temperature);

    const double mass_flux_value =
      evaporation_coefficient * saturation_pressure_value *
      std::sqrt(M_x_2PI_inv_x_R_inv * temperature_inv);

    return mass_flux_value;
  }

  /**
   * @brief vector_mass_flux Calculates the values of the evaporation mass flux
   * for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param mass_flux_vector Vectors of the mass flux values
   */
  void
  vector_mass_flux(const std::map<field, std::vector<double>> &field_vectors,
                   std::vector<double> &mass_flux_vector) override
  {
    const std::vector<double> &temperature =
      field_vectors.at(field::temperature);

    const unsigned int n_pts = mass_flux_vector.size();

    std::vector<double> saturation_pressure_vector(n_pts);
    vector_saturation_pressure(field_vectors, saturation_pressure_vector);

    for (unsigned int i = 0; i < mass_flux_vector.size(); ++i)
      {
        const double temperature_inv = 1.0 / temperature[i];

        mass_flux_vector[i] = evaporation_coefficient *
                              saturation_pressure_vector[i] *
                              std::sqrt(M_x_2PI_inv_x_R_inv * temperature_inv);
      }
  }

  /**
   * @brief heat_flux Calculates the value of the evaporation heat flux.
   * @param fields_value Value of the various field on which the flux may depend.
   * @return value of the heat flux calculated with the fields_value
   */
  double
  heat_flux(const std::map<field, double> &fields_value) override
  {
    const double mass_flux_value = mass_flux(fields_value);

    return mass_flux_value * latent_heat_evaporation;
  }

  /**
   * @brief vector_heat_flux Calculates the values of the evaporation heat flux
   * for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param heat_flux_vector Vectors of the heat flux values
   */
  void
  vector_heat_flux(const std::map<field, std::vector<double>> &field_vectors,
                   std::vector<double> &heat_flux_vector) override
  {
    const unsigned int n_pts = heat_flux_vector.size();

    std::vector<double> mass_flux_vector(n_pts);
    vector_mass_flux(field_vectors, mass_flux_vector);

    for (unsigned int i = 0; i < n_pts; ++i)
      {
        heat_flux_vector[i] = mass_flux_vector[i] * latent_heat_evaporation;
      }
  }

  /**
   * @brief heat_flux_jacobian Calculates the jacobian (the partial derivative)
   * of the evaporation heat flux with respect to a field
   * @param field_values Value of the various fields on which the flux may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the heat flux with respect to the field.
   */
  double
  heat_flux_jacobian(const std::map<field, double> &fields_value,
                     const field                    id) override
  {
    const double temperature_inv = 1.0 / fields_value.at(field::temperature);

    if (id == field::temperature)
      return latent_heat_evaporation * ambient_pressure *
             evaporation_coefficient * temperature_inv *
             std::sqrt(M_x_2PI_inv_x_R_inv * temperature_inv) *
             std::exp(-L_vap_x_M_x_R_inv *
                      (temperature_inv - boiling_temperature_inv)) *
             (L_vap_x_M_x_R_inv * temperature_inv - 0.5);
    else
      return 0.0;
  }

  /**
   * @brief vector_heat_flux_jacobian Calculate the derivative of the evaporation
   * heat flux with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate
   * the flux
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated
   * @param jacobian Vector of the values of the derivative of the heat flux with
   * respect to the field id
   */
  void
  vector_heat_flux_jacobian(
    const std::map<field, std::vector<double>> &field_vectors,
    const field                                 id,
    std::vector<double>                        &jacobian_vector) override
  {
    if (id == field::temperature)
      {
        const std::vector<double> &temperature =
          field_vectors.at(field::temperature);
        for (unsigned int i = 0; i < jacobian_vector.size(); ++i)
          {
            const double temperature_inv = 1.0 / temperature[i];
            jacobian_vector[i] =
              latent_heat_evaporation * ambient_pressure *
              evaporation_coefficient * temperature_inv *
              std::sqrt(M_x_2PI_inv_x_R_inv * temperature_inv) *
              std::exp(-L_vap_x_M_x_R_inv *
                       (temperature_inv - boiling_temperature_inv)) *
              (L_vap_x_M_x_R_inv * temperature_inv - 0.5);
          }
      }
    else
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }
  /**
   * @brief momentum_flux Calculates the value of the evaporation recoil pressure.
   * @param fields_value Value of the various field on which the recoil pressure
   * may depend.
   * @return value of the recoil pressure calculated with the fields_value
   */
  double
  momentum_flux(const std::map<field, double> &fields_value) override
  {
    const double saturation_pressure_value = saturation_pressure(fields_value);
    const double mass_flux_value           = mass_flux(fields_value);
    const double temperature_inv = 1.0 / fields_value.at(field::temperature);

    const double vapor_saturation_density_value =
      saturation_pressure_value * R_inv * temperature_inv;

    // rho_vap = 0.31*rho_sat according to Anisimov and Khokhlov 1995
    const double vapor_density_inv =
      1.0 / (0.31 * vapor_saturation_density_value);

    return -mass_flux_value * mass_flux_value *
             (liquid_density_inv - vapor_density_inv) +
           recoil_pressure_coefficient * saturation_pressure_value;
  }

  /**
   * @brief vector_momentum_flux Calculates the values of the evaporation recoil
   * pressure for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param momentum_flux_vector Vectors of the recoil pressure values
   */
  void
  vector_momentum_flux(
    const std::map<field, std::vector<double>> &field_vectors,
    std::vector<double>                        &momentum_flux_vector) override
  {
    const unsigned int         n_pts = momentum_flux_vector.size();
    const std::vector<double> &temperature =
      field_vectors.at(field::temperature);

    std::vector<double> saturation_pressure_vector(n_pts);
    vector_saturation_pressure(field_vectors, saturation_pressure_vector);

    std::vector<double> mass_flux_vector(n_pts);
    vector_mass_flux(field_vectors, mass_flux_vector);

    for (unsigned int i = 0; i < momentum_flux_vector.size(); ++i)
      {
        const double temperature_inv = 1.0 / temperature[i];

        const double vapor_saturation_density_value =
          saturation_pressure_vector[i] * R_inv * temperature_inv;

        // rho_vap = 0.31*rho_sat according to Anisimov and Khokhlov 1995
        const double vapor_density_inv =
          1.0 / (0.31 * vapor_saturation_density_value);

        momentum_flux_vector[i] =
          -mass_flux_vector[i] * mass_flux_vector[i] *
            (liquid_density_inv - vapor_density_inv) +
          recoil_pressure_coefficient * saturation_pressure_vector[i];
      }
  }

  /**
   * @brief momentum_flux_jacobian Calculates the jacobian (the partial derivative)
   * of the evaporation recoil pressure with respect to a field
   * @param field_values Value of the various fields on which the recoil pressure
   * may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the recoil pressure with respect
   * to the field.
   */
  double
  momentum_flux_jacobian(const std::map<field, double> & /*field_values*/,
                         const field /*id*/) override
  {
    return 0.0;
  }

  /**
   * @brief vector_momentum_flux_jacobian Calculate the derivative of the
   * evaporation recoil pressure with respect to a field
   * @param field_vectors Vector for the values of the fields on which the recoil
   * pressure may depend.
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated
   * @param jacobian Vector of the value of the derivative of the recoil pressure
   * with respect to the field id
   */
  void
  vector_momentum_flux_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }

private:
  const double evaporation_coefficient;
  const double recoil_pressure_coefficient;

  const double molar_mass;
  const double boiling_temperature_inv;
  const double latent_heat_evaporation;
  const double ambient_pressure;
  const double ambient_gas_density_inv;
  const double liquid_density_inv;

  // 1/R
  double R_inv;

  // L_vap*M/R
  double L_vap_x_M_x_R_inv;

  // M/(2*PI*R)
  double M_x_2PI_inv_x_R_inv;
};

#endif
