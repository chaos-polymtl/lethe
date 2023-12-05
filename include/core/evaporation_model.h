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
#include <core/physical_property_model.h>

using namespace dealii;

/**
 * @brief EvaporationModel. Abstract class that allows to calculate the
 * evaporation mass, heat and momentum fluxes. This model computes all the terms
 * required for the assembly of the governing equations, i.e., Navier-Stokes
 * momentum and energy, when evaporation is considered. The former terms are
 * the momentum flux (i.e., pressure) acting on the evaporative front in the
 * normal direction, the heat flux, acting at the same location as a cooling
 * term, and the fluxes' jacobian, if required.
 *
 * The EvaporativeModel is the parent class, and its children are specific
 * types of evaporative models, e.g., constant, temperature dependent. The
 * selected model is instantiated in the assemblers' constructor
 * (NavierStokesVOFAssemblerEvaporation and/or
 * HeatTransferAssemblerVOFEvaporation), using the model_cast method according
 * to the evaporation model parameters.
 *
 * The assemblers call the required method of the EvaporationModel for either
 * the assembly of the rhs or matrix terms. Hence, no matter the selected
 * evaporation model, the same assembler is called, avoiding repetitions of
 * the assembler classes.
 *
 */
class EvaporationModel
{
public:
  EvaporationModel()
  {}

  /**
   * @brief Instantiates and returns a pointer to an EvaporationModel
   * object by casting it to the proper child class.
   *
   * @param evaporation_parameters Evaporation model parameters.
   */
  static std::shared_ptr<EvaporationModel>
  model_cast(const Parameters::Evaporation &evaporation_parameters);


  /**
   * @brief Returns true if the EvaporationModel depends on a field, false if not.
   * @param id Identifier of the field for which the function is called.
   */
  inline bool
  depends_on(const field &id)
  {
    return model_depends_on[id];
  }

  /**
   * @brief mass_flux Calculates the value of the evaporation mass flux.
   * @param field_values Map storing the values of the various fields on which the
   * mass flux may depend.
   * @return Value of the flux calculated with the field_values.
   */
  virtual double
  mass_flux(const std::map<field, double> &field_values) = 0;

  /**
   * @brief mass_flux Calculates the values of the evaporation mass flux
   * for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the mass flux may depend. Each vector contains the field values at
   * the points of interest.
   * @param mass_flux_vector Vector of mass flux values.
   */
  virtual void
  mass_flux(const std::map<field, std::vector<double>> &field_vectors,
                   std::vector<double> &mass_flux_vector) = 0;

  /**
   * @brief heat_flux Calculates the value of the evaporation heat flux.
   * @param field_values Map storing the values of the various fields on which
   * the heat flux may depend.
   * @return Value of the heat flux calculated with the field_values
   */
  virtual double
  heat_flux(const std::map<field, double> &field_values) = 0;

  /**
   * @brief heat_flux Calculates the values of the evaporation heat flux
   * for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the heat flux may depend. Each vector contains the field values at
   * the points of interest.
   * @param heat_flux_vector Vector of heat flux values.
   */
  virtual void
  heat_flux(const std::map<field, std::vector<double>> &field_vectors,
                   std::vector<double> &heat_flux_vector) = 0;

  /**
   * @brief heat_flux_jacobian Calculates the jacobian (the partial derivative)
   * of the evaporation heat flux with respect to a field.
   * @param field_values Map storing the values of the various fields on which the
   * heat flux evaluated at a point may depend.
   * @param id Identifier of the field with respect to which the jacobian
   * should be calculated.
   */
  virtual double
  heat_flux_jacobian(const std::map<field, double> &field_values,
                     const field                    id) = 0;

  /**
   * @brief heat_flux_jacobian Calculates the derivative of the evaporation
   * heat flux with respect to a field for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the heat flux may depend. Each vector contains the field values at
   * the points of interest.
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated.
   * @param jacobian_vector Vector of the values of the derivative of the heat
   * flux with respect to the field.
   */
  virtual void
  heat_flux_jacobian(
    const std::map<field, std::vector<double>> &field_vectors,
    const field                                 id,
    std::vector<double>                        &jacobian_vector) = 0;

  /**
   * @brief momentum_flux Calculates the value of the evaporation momentum flux.
   * @param field_values Map storing the values of the various fields on which the
   * momentum flux may depend.
   * @return Value of the momentum flux calculated with the field_values.
   */
  virtual double
  momentum_flux(const std::map<field, double> &field_values) = 0;

  /**
   * @brief momentum_flux Calculates the values of the evaporation momentum
   * flux for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the evaporation momentum flux may depend. Each vector contains the
   * field values at the points of interest.
   * @param momentum_flux_vector Vector of momentum flux values.
   */
  virtual void
  momentum_flux(
    const std::map<field, std::vector<double>> &field_vectors,
    std::vector<double>                        &momentum_flux_vector) = 0;

  /**
   * @brief momentum_flux_jacobian Calculates the jacobian (the partial derivative)
   * of the evaporation momentum flux with respect to a field.
   * @param field_values Map storing the values of the various fields on which the
   * momentum flux may depend.
   * @param id Identifier of the field with respect to which the jacobian
   * should be calculated.
   * @return Value of the partial derivative of the momentum flux with respect
   * to the field.
   */
  virtual double
  momentum_flux_jacobian(const std::map<field, double> &field_values,
                         const field                    id) = 0;

  /**
   * @brief momentum_flux_jacobian Calculates the derivative of the
   * evaporation momentum flux with respect to a field for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the momentum flux may depend. Each vector contains the field
   * values at the points of interest.
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated.
   * @param jacobian_vector Vector of the value of the derivative of the momentum
   * flux with respect to the field.
   */
  virtual void
  momentum_flux_jacobian(
    const std::map<field, std::vector<double>> &field_vectors,
    const field                                 id,
    std::vector<double>                        &jacobian_vector) = 0;

protected:
  // Map that indicates on which fields the model depends on
  std::map<field, bool> model_depends_on;
};

/**
 * @brief EvaporationModelConstant. Class that allows to calculate the
 * evaporation mass, heat, and momentum fluxes when a constant mass flux is
 * considered.
 *
 * @param p_evaporation evaporation model parameters.
 */
class EvaporationModelConstant : public EvaporationModel
{
public:
  EvaporationModelConstant(const Parameters::Evaporation &p_evaporation)
    : n_evaporation(p_evaporation.evaporation_mass_flux)
    , latent_heat_evaporation(p_evaporation.latent_heat_evaporation)
    , ambient_gas_density_inv(1.0 / p_evaporation.ambient_gas_density)
    , liquid_density_inv(1.0 / p_evaporation.liquid_density)
  {
    model_depends_on[temperature] = false;
  }


  /**
   * @brief mass_flux Calculates the value of the evaporation mass flux.
   * @param field_values Map storing the values of the various fields on which the
   * mass flux may depend.
   * @return Value of the flux calculated with the field_values.
   */
  double
  mass_flux(const std::map<field, double> & /*field_values*/) override
  {
    return n_evaporation;
  }

  /**
   * @brief mass_flux Calculates the values of the evaporation mass flux
   * for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the mass flux may depend. Each vector contains the field values at
   * the points of interest.
   * @param mass_flux_vector Vector of mass flux values.
   */
  void
  mass_flux(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    std::vector<double> &mass_flux_vector) override
  {
    std::fill(mass_flux_vector.begin(), mass_flux_vector.end(), n_evaporation);
  }

  /**
   * @brief heat_flux Calculates the value of the evaporation heat flux.
   * @param field_values Map storing the values of the various fields on which the
   * heat flux may depend.
   * @return Value of the heat flux calculated with the field_values.
   */
  double
  heat_flux(const std::map<field, double> & /*field_values*/) override
  {
    return n_evaporation * latent_heat_evaporation;
  }

  /**
   * @brief heat_flux Calculates the values of the evaporation heat flux
   * for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the heat flux may depend. Each vector contains the field values at
   * the points of interest.
   * @param heat_flux_vector Vectors of the heat flux values.
   */
  void
  heat_flux(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    std::vector<double> &heat_flux_vector) override
  {
    std::fill(heat_flux_vector.begin(),
              heat_flux_vector.end(),
              n_evaporation * latent_heat_evaporation);
  }

  /**
   * @brief heat_flux_jacobian Calculates the jacobian (the partial derivative)
   * of the evaporation heat flux with respect to a field.
   * @param field_values Map storing the values of the various fields on which the
   * heat flux evaluated at a point may depend.
   * @param id Identifier of the field with respect to which the jacobian
   * should be calculated.
   */
  double
  heat_flux_jacobian(const std::map<field, double> & /*field_values*/,
                     const field /*id*/) override
  {
    return 0.0;
  }

  /**
   * @brief heat_flux_jacobian Calculates the derivative of the evaporation
   * heat flux with respect to a field.
   * @param field_vectors Vector for the values of the fields on which the heat
   * flux may depend. Each vector contains the field values at the points of
   * interest.
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated.
   * @param jacobian_vector Vector of the values of the derivative of the heat
   * flux with respect to the field.
   */
  void
  heat_flux_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }

  /**
   * @brief momentum_flux Calculates the value of the evaporation momentum flux.
   * @param field_values Map storing the values of the various fields on which the
   * momentum flux may depend.
   * @return Value of the momentum flux calculated with the field_values
   */
  double
  momentum_flux(const std::map<field, double> & /*field_values*/) override
  {
    return -n_evaporation * n_evaporation *
           (liquid_density_inv - ambient_gas_density_inv);
  }

  /**
   * @brief momentum_flux Calculates the values of the evaporation momentum
   * flux for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the momentum flux may depend. Each vector contains the field
   * values at the points of interest.
   * @param momentum_flux_vector Vector of momentum flux values.
   */
  void
  momentum_flux(
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
   * of the evaporation momentum flux with respect to a field.
   * @param field_values Map storing the values of the various fields on which the
   * momentum flux may depend.
   * @param id Identifier of the field with respect to which the jacobian
   * should be calculated
   * @return Value of the partial derivative of the momentum flux with respect
   * to the field.
   */
  double
  momentum_flux_jacobian(const std::map<field, double> & /*field_values*/,
                         const field /*id*/) override
  {
    return 0.0;
  }

  /**
   * @brief momentum_flux_jacobian Calculates the derivative of the
   * evaporation momentum flux with respect to a field
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the momentum flux may depend. Each vector contains the field
   * values at the points of interest.
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated.
   * @param jacobian_vector Vector of the values of the derivative of the momentum
   * flux with respect to the field.
   */
  void
  momentum_flux_jacobian(
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
 * @brief EvaporationModelTemperature. Class that allows to calculate the
 * evaporation mass, heat, and momentum fluxes  when a temperature-dependent
 * mass flux is considered.
 *
 * @param p_evaporation evaporation model parameters.
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
  }

  /**
   * @brief saturation_pressure Calculates the value of the saturation pressure.
   * @param field_values Map storing the values of the various field on which the
   * saturated vapor pressure evaluated at a point may depend.
   * @return Value of the saturation pressure calculated with the field_values.
   */
  double
  saturation_pressure(const std::map<field, double> &field_values)
  {
    const double temperature_inv =
      1.0 / (field_values.at(field::temperature) + 1e-16);

    const double R_inv = 1.0 / 8.3145;

    const double L_vap_x_M_x_R_inv =
      latent_heat_evaporation * molar_mass * R_inv;

    return ambient_pressure *
           std::exp(-L_vap_x_M_x_R_inv *
                    (temperature_inv - boiling_temperature_inv));
  }

  /**
   * @brief mass_flux Calculates the values of the saturation pressure
   * for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the saturated vapor pressure may depend. Each vector contains the
   * field values at the points of interest.
   * @param saturation_pressure_vector Vector of the saturation pressure values.
   */
  void
  saturation_pressure(
    const std::map<field, std::vector<double>> &field_vectors,
    std::vector<double>                        &saturation_pressure_vector)
  {
    const std::vector<double> &temperature =
      field_vectors.at(field::temperature);

    const double R_inv = 1.0 / 8.3145;

    const double L_vap_x_M_x_R_inv =
      latent_heat_evaporation * molar_mass * R_inv;

    for (unsigned int i = 0; i < saturation_pressure_vector.size(); ++i)
      {
        const double temperature_inv = 1.0 / (temperature[i] + 1e-16);

        saturation_pressure_vector[i] =
          ambient_pressure *
          std::exp(-L_vap_x_M_x_R_inv *
                   (temperature_inv - boiling_temperature_inv));
      }
  }

  /**
   * @brief mass_flux Calculates the value of the evaporation mass flux.
   * @param field_values Map storing the values of the various fields on which the
   * evaporation mass flux may depend.
   * @return Value of the evaporation mass flux calculated with the field_values.
   */
  double
  mass_flux(const std::map<field, double> &field_values) override
  {
    const double saturation_pressure_value = saturation_pressure(field_values);
    const double temperature_inv =
      1.0 / (field_values.at(field::temperature) + 1e-16);

    const double R_inv               = 1.0 / 8.3145;
    const double M_x_2PI_inv_x_R_inv = molar_mass / (2.0 * M_PI) * R_inv;

    const double mass_flux_value =
      evaporation_coefficient * saturation_pressure_value *
      std::sqrt(M_x_2PI_inv_x_R_inv * temperature_inv);

    return mass_flux_value;
  }

  /**
   * @brief mass_flux Calculates the values of the evaporation mass flux
   * for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the evaporation mass flux may depend. Each vector contains the
   * field values at the points of interest.
   * @param mass_flux_vector Vector of mass flux values.
   */
  void
  mass_flux(const std::map<field, std::vector<double>> &field_vectors,
                   std::vector<double> &mass_flux_vector) override
  {
    const std::vector<double> &temperature =
      field_vectors.at(field::temperature);

    const double R_inv               = 1.0 / 8.3145;
    const double M_x_2PI_inv_x_R_inv = molar_mass / (2.0 * M_PI) * R_inv;

    const unsigned int n_pts = mass_flux_vector.size();

    std::vector<double> saturation_pressure_vector(n_pts);
    saturation_pressure(field_vectors, saturation_pressure_vector);

    for (unsigned int i = 0; i < n_pts; ++i)
      {
        const double temperature_inv = 1.0 / (temperature[i] + 1e-16);

        mass_flux_vector[i] = evaporation_coefficient *
                              saturation_pressure_vector[i] *
                              std::sqrt(M_x_2PI_inv_x_R_inv * temperature_inv);
      }
  }

  /**
   * @brief heat_flux Calculates the value of the evaporation heat flux.
   * @param field_values Map storing the values of the various fields on which the
   * heat flux may depend.
   * @return Value of the heat flux calculated with the field_values.
   */
  double
  heat_flux(const std::map<field, double> &field_values) override
  {
    const double mass_flux_value = mass_flux(field_values);

    return mass_flux_value * latent_heat_evaporation;
  }

  /**
   * @brief heat_flux Calculates the values of the evaporation heat flux
   * for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the heat flux may depend. Each vector contains the field values at
   * the points of interest.
   * @param heat_flux_vector Vector of heat flux values.
   */
  void
  heat_flux(const std::map<field, std::vector<double>> &field_vectors,
                   std::vector<double> &heat_flux_vector) override
  {
    const unsigned int n_pts = heat_flux_vector.size();

    std::vector<double> mass_flux_vector(n_pts);
    mass_flux(field_vectors, mass_flux_vector);

    for (unsigned int i = 0; i < n_pts; ++i)
      {
        heat_flux_vector[i] = mass_flux_vector[i] * latent_heat_evaporation;
      }
  }

  /**
   * @brief heat_flux_jacobian Calculates the jacobian (the partial derivative)
   * of the evaporation heat flux with respect to a field.
   * @param field_values Map storing the values of the various fields on which the
   * heat flux evaluated at a point may depend.
   * @param id Identifier of the field with respect to which the jacobian
   * should be calculated.
   * @return Value of the partial derivative of the heat flux with respect to
   * the field.
   */
  double
  heat_flux_jacobian(const std::map<field, double> &field_values,
                     const field                    id) override
  {
    const double temperature_inv =
      1.0 / (field_values.at(field::temperature) + 1e-16);

    const double R_inv = 1.0 / 8.3145;
    const double L_vap_x_M_x_R_inv =
      latent_heat_evaporation * molar_mass * R_inv;
    const double M_x_2PI_inv_x_R_inv = molar_mass / (2.0 * M_PI) * R_inv;

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
   * @brief heat_flux_jacobian Calculates the derivative of the evaporation
   * heat flux with respect to a field.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which heat flux evaluated at a point may depend. Each vector contains
   * the field values at the points of interest.
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated.
   * @param jacobian_vector Vector of the values of the derivative of the heat flux
   * with respect to the field.
   */
  void
  heat_flux_jacobian(
    const std::map<field, std::vector<double>> &field_vectors,
    const field                                 id,
    std::vector<double>                        &jacobian_vector) override
  {
    const double R_inv = 1.0 / 8.3145;
    const double L_vap_x_M_x_R_inv =
      latent_heat_evaporation * molar_mass * R_inv;
    const double M_x_2PI_inv_x_R_inv = molar_mass / (2.0 * M_PI) * R_inv;

    if (id == field::temperature)
      {
        const std::vector<double> &temperature =
          field_vectors.at(field::temperature);
        for (unsigned int i = 0; i < jacobian_vector.size(); ++i)
          {
            const double temperature_inv = 1.0 / (temperature[i] + 1e-16);
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
   * @brief momentum_flux Calculates the value of the evaporation momentum flux.
   * @param field_values Map storing the values of the various fields on which the
   * momentum flux may depend.
   * @return Value of the momentum flux calculated with the field_values
   */
  double
  momentum_flux(const std::map<field, double> &field_values) override
  {
    const double saturation_pressure_value = saturation_pressure(field_values);
    const double mass_flux_value           = mass_flux(field_values);
    const double temperature_inv =
      1.0 / (field_values.at(field::temperature) + 1e-16);

    const double R_inv = 1.0 / 8.3145;

    const double vapor_saturation_density_value =
      molar_mass * saturation_pressure_value * R_inv * temperature_inv;

    // rho_vap = 0.31*rho_sat according to Anisimov and Khokhlov 1995
    const double vapor_density_inv =
      1.0 / (0.31 * vapor_saturation_density_value);

    const double pressure =
      -mass_flux_value * mass_flux_value *
        (liquid_density_inv - vapor_density_inv) +
      recoil_pressure_coefficient * saturation_pressure_value;

    return std::max(pressure - ambient_pressure, 0.0);
  }

  /**
   * @brief momentum_flux Calculates the values of the evaporation momentum
   * flux for multiple points.
   * @param field_vectors Map storing the vectors of the values for each fields
   * on which the momentum flux may depend. Each vector contains the field
   * values at the points of interest.
   * @param momentum_flux_vector Vectors of the momentum flux values.
   */
  void
  momentum_flux(
    const std::map<field, std::vector<double>> &field_vectors,
    std::vector<double>                        &momentum_flux_vector) override
  {
    const unsigned int         n_pts = momentum_flux_vector.size();
    const std::vector<double> &temperature =
      field_vectors.at(field::temperature);

    const double R_inv = 1.0 / 8.3145;

    std::vector<double> saturation_pressure_vector(n_pts);
    saturation_pressure(field_vectors, saturation_pressure_vector);

    std::vector<double> mass_flux_vector(n_pts);
    mass_flux(field_vectors, mass_flux_vector);

    for (unsigned int i = 0; i < n_pts; ++i)
      {
        const double temperature_inv = 1.0 / (temperature[i] + 1e-16);

        const double vapor_saturation_density_value =
          molar_mass * saturation_pressure_vector[i] * R_inv * temperature_inv;

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
   * of the evaporation momentum flux with respect to a field.
   * @param field_values Map storing the values of the various fields on which the
   * momentum flux may depend.
   * @param id Identifier of the field with respect to which the jacobian
   * should be calculated
   * @return Value of the partial derivative of the momentum flux with respect
   * to the field.
   */
  double
  momentum_flux_jacobian(const std::map<field, double> & /*field_values*/,
                         const field /*id*/) override
  {
    return 0.0;
  }

  /**
   * @brief momentum_flux_jacobian Calculates the derivative of the
   * evaporation momentum flux with respect to a field.
   * @param field_vectors Map storing the vectors of the values for each fields
   * flux may depend. Each vector contains the field values at the points of
   * interest.
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated.
   * @param jacobian_vector Vector of the values of the derivative of the momentum
   * flux with respect to the field.
   */
  void
  momentum_flux_jacobian(
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
};

#endif
