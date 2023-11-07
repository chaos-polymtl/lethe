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
               std::vector<double>                        &mass_flux_vector) = 0;
               
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
               std::vector<double>                        &heat_flux_vector) = 0;

  /**
   * @brief heat_flux_jacobian Calculates the jacobian (the partial derivative) 
   * of the evaporation heat flux with respect to a field
   * @param field_values Value of the various fields on which the flux may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the heat flux with respect to the field.
   */
  virtual double
  heat_flux_jacobian(const std::map<field, double> &field_values, const field id) = 0;

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
  vector_heat_flux_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) = 0;

  /**
   * @brief recoil_pressure Calculates the value of the evaporation recoil pressure.
   * @param fields_value Value of the various field on which the recoil pressure 
   * may depend.
   * @return value of the recoil pressure calculated with the fields_value
   */
  virtual double
  recoil_pressure(const std::map<field, double> &fields_value) = 0;
  
  /**
   * @brief vector_recoil_pressure Calculates the values of the evaporation recoil 
   * pressure for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param recoil_pressure_vector Vectors of the recoil pressure values
   */
  virtual void
  vector_recoil_pressure(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double>                        &recoil_pressure_vector) = 0;

  /**
   * @brief recoil_pressure_jacobian Calculates the jacobian (the partial derivative) 
   * of the evaporation recoil pressure with respect to a field
   * @param field_values Value of the various fields on which the recoil pressure 
   * may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the recoil pressure with respect 
   * to the field.
   */
  virtual double
  recoil_pressure_jacobian(const std::map<field, double> &field_values, const field id) = 0;

  /**
   * @brief vector_recoil_pressure_jacobian Calculate the derivative of the  
   * evaporation recoil pressure with respect to a field
   * @param field_vectors Vector for the values of the fields on which the recoil 
   * pressure may depend.
   * @param id Identifier of the field with respect to which a derivative should 
   * be calculated
   * @param jacobian Vector of the value of the derivative of the recoil pressure
   * with respect to the field id
   */
  virtual void
  vector_recoil_pressure_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) = 0;

protected:
  // Map to indicate on which variables the model depends on
  std::map<field, bool> model_depends_on;
};

/**
 * @brief EvaporationModel. Abstract class that allows to calculate the
 * evaporation mass, heat and momentum fluxes.
 */
class EvaporationModelConstant : public EvaporationModel
{
public:
  EvaporationModelConstant(const Parameters::Evaporation &p_evaporation): n_evaporation(p_evaporation.n_evaporation), latent_heat_evaporation(p_evaporation.latent_heat_evaporation), ambient_gas_density_inv(1.0/p_evaporation.ambient_gas_density), liquid_density_inv(1.0/p_evaporation.liquid_density)
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
  vector_mass_flux(const std::map<field, std::vector<double>> & /*fields_value*/,
               std::vector<double>                        &mass_flux_vector) override
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
    return n_evaporation*latent_heat_evaporation;
  }
  
  /**
   * @brief vector_heat_flux Calculates the values of the evaporation heat flux 
   * for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param heat_flux_vector Vectors of the heat flux values
   */
  void
  vector_heat_flux(const std::map<field, std::vector<double>> &/*field_vectors*/,
               std::vector<double>                        &heat_flux_vector) override
  {
    std::fill(heat_flux_vector.begin(), heat_flux_vector.end(),  n_evaporation*latent_heat_evaporation);
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
  heat_flux_jacobian(const std::map<field, double> & /*fields_value*/, const field /*id*/) override
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
  vector_heat_flux_jacobian(const std::map<field, std::vector<double>> &/*field_vectors*/,
                  const field                                 /*id*/,
                  std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }
  /**
   * @brief recoil_pressure Calculates the value of the evaporation recoil pressure.
   * @param fields_value Value of the various field on which the recoil pressure 
   * may depend.
   * @return value of the recoil pressure calculated with the fields_value
   */
  double
  recoil_pressure(const std::map<field, double> & /*fields_value*/) override
  {
    return -n_evaporation*n_evaporation*(liquid_density_inv - ambient_gas_density_inv);
  }
  
  /**
   * @brief vector_recoil_pressure Calculates the values of the evaporation recoil 
   * pressure for multiple points
   * @param field_vectors Value of the various fields on which the flux may depend.
   * @param recoil_pressure_vector Vectors of the recoil pressure values
   */
  void
  vector_recoil_pressure(const std::map<field, std::vector<double>> &/*field_vectors*/,
               std::vector<double>                        &recoil_pressure_vector) override
  { 
    
    const double recoil_pressure_value = -n_evaporation*n_evaporation*(liquid_density_inv - ambient_gas_density_inv);
    
    std::fill(recoil_pressure_vector.begin(), recoil_pressure_vector.end(), recoil_pressure_value);
  }

  /**
   * @brief recoil_pressure_jacobian Calculates the jacobian (the partial derivative) 
   * of the evaporation recoil pressure with respect to a field
   * @param field_values Value of the various fields on which the recoil pressure 
   * may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the recoil pressure with respect 
   * to the field.
   */
  double
  recoil_pressure_jacobian(const std::map<field, double> &/*field_values*/, const field /*id*/) override
  {
    return 0.0;
  }

  /**
   * @brief vector_recoil_pressure_jacobian Calculate the derivative of the  
   * evaporation recoil pressure with respect to a field
   * @param field_vectors Vector for the values of the fields on which the recoil 
   * pressure may depend.
   * @param id Identifier of the field with respect to which a derivative should 
   * be calculated
   * @param jacobian Vector of the value of the derivative of the recoil pressure
   * with respect to the field id
   */
  void
  vector_recoil_pressure_jacobian(const std::map<field, std::vector<double>> & /*field_vectors*/,
                  const field                                 /*id*/,
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

#endif
