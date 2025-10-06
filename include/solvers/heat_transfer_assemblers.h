// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_heat_transfer_assemblers_h
#define lethe_heat_transfer_assemblers_h

#include <core/evaporation_model.h>
#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/heat_transfer_scratch_data.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/physics_assemblers.h>

/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for heat transfer
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 *
 * @ingroup assemblers
 */
template <int dim>
using HeatTransferAssemblerBase =
  PhysicsAssemblerBase<HeatTransferScratchData<dim>, StabilizedMethodsCopyData>;


/**
 * @brief Class that assembles the core of the heat transfer equation.
 * This class assembles the weak form of:
 * \f$ - k * \nabla^2 T + \rho * cp * \mathbf{u} * \nabla T - f - \tau :
 * \nabla \mathbf{u} =0 \f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerCore : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerCore(
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;


  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the heat transfer solver.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerBDF : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerBDF(
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

  const bool GGLS = true;

  const std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief Class that assembles the Robin boundary condition for the heat
 * transfer solver.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 * @param p_boundary_conditions_ht HTBoundaryConditions object that hold
 * boundary condition information for the Heat-Transfer solver
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerRobinBC : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerRobinBC(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const BoundaryConditions::HTBoundaryConditions<dim>
      &p_boundary_conditions_ht)
    : simulation_control(simulation_control)
    , boundary_conditions_ht(p_boundary_conditions_ht)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

  const std::shared_ptr<SimulationControl>             simulation_control;
  const BoundaryConditions::HTBoundaryConditions<dim> &boundary_conditions_ht;
};


/**
 * @brief Class that assembles the viscous dissipation for the heat transfer
 * solver.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerViscousDissipation
  : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerViscousDissipation(
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief Class that assembles the viscous dissipation for the heat transfer
 * solver, for the specific case of VOF simulations. The only difference
 * compared to the regular one is that the viscous dissipation can be applied in
 * one of the fluids rather than both, through the viscous_dissipative_fluid
 * parameter.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 * @param p_viscous_dissipative_fluid A FluidIndicator enum element indicating
 * the selected viscous dissipative fluid(s).
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerViscousDissipationVOF
  : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerViscousDissipationVOF(
    const std::shared_ptr<SimulationControl> &simulation_control,
    Parameters::FluidIndicator                p_viscous_dissipative_fluid)
    : simulation_control(simulation_control)
    , viscous_dissipative_fluid(p_viscous_dissipative_fluid)
  {}

  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

protected:
  const std::shared_ptr<SimulationControl> simulation_control;
  Parameters::FluidIndicator               viscous_dissipative_fluid;
};


/**
 * @brief Class that assembles the Discontinuity-Capturing Directional
 * Dissipation stabilization term for the heat transfer
 * solver. For more information see Tezduyar, T. E. (2003). Computation of
 * moving boundaries and interfaces and stabilization parameters. International
 * Journal for Numerical Methods in Fluids, 43(5), 555-575. The implementation
 * is based on equations (70) and (79), which are adapted for the heat transfer
 * solver.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 *
 * @ingroup assemblers
 */

template <int dim>
class HeatTransferAssemblerDCDDstabilization
  : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerDCDDstabilization(
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief Class that assembles the laser heating as a volumetric source for
 * the heat transfer solver. Exponentially decaying model is used to simulate
 * the laser heat source: "Liu, S., Zhu, H., Peng, G., Yin, J. and Zeng, X.,
 * 2018. Microstructure prediction of selective laser melting AlSi10Mg
 * using finite element analysis. Materials & Design, 142, pp.319-328."
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 * @param p_laser_parameters Shared pointer of the laser parameters
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerLaserExponentialDecay
  : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerLaserExponentialDecay(
    const std::shared_ptr<SimulationControl> &simulation_control,
    std::shared_ptr<Parameters::Laser<dim>>   p_laser_parameters)
    : simulation_control(simulation_control)
    , laser_parameters(p_laser_parameters)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

protected:
  const std::shared_ptr<SimulationControl> simulation_control;
  std::shared_ptr<Parameters::Laser<dim>>  laser_parameters;
};


/**
 * @brief Class that assembles the laser heating as a surface flux for the
 * heat transfer solver when VOF is enabled. The laser heat flux is
 * applied at the VOF interface (where the phase gradient is non-null).
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 * @param p_laser_parameters Shared pointer of the laser parameters
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerLaserGaussianHeatFluxVOFInterface
  : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerLaserGaussianHeatFluxVOFInterface(
    const std::shared_ptr<SimulationControl> &simulation_control,
    std::shared_ptr<Parameters::Laser<dim>>   p_laser_parameters)
    : simulation_control(simulation_control)
    , laser_parameters(p_laser_parameters)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

protected:
  const std::shared_ptr<SimulationControl> simulation_control;
  std::shared_ptr<Parameters::Laser<dim>>  laser_parameters;
};


/**
 * @brief Class that assembles the laser heating as a uniform surface flux for
 * the heat transfer solver when VOF is enabled. The laser heat flux is
 * applied at the VOF interface (where the phase gradient is non-null).
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 * @param p_laser_parameters Shared pointer of the laser parameters
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerLaserUniformHeatFluxVOFInterface
  : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerLaserUniformHeatFluxVOFInterface(
    const std::shared_ptr<SimulationControl> &simulation_control,
    std::shared_ptr<Parameters::Laser<dim>>   p_laser_parameters)
    : simulation_control(simulation_control)
    , laser_parameters(p_laser_parameters)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

protected:
  const std::shared_ptr<SimulationControl> simulation_control;
  std::shared_ptr<Parameters::Laser<dim>>  laser_parameters;
};


/**
 * @brief Class that assembles the laser heating as a volumetric source for
 * the heat transfer solver when VOF is enabled. Exponentially decaying model is
 * used to simulate the laser heat source: "Liu, S., Zhu, H., Peng, G.,
 * Yin, J. and Zeng, X., 2018. Microstructure prediction of selective
 * laser melting AlSi10Mg using finite element analysis. Materials &
 * Design, 142, pp.319-328." The laser heat source is only applied in the metal
 * (when phase value is non-null) using the phase value alpha as a multiplying
 * factor on the laser heat source.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 * @param p_laser_parameters Shared pointer of the laser parameters
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerLaserExponentialDecayVOF
  : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerLaserExponentialDecayVOF(
    const std::shared_ptr<SimulationControl> &simulation_control,
    std::shared_ptr<Parameters::Laser<dim>>   p_laser_parameters)
    : simulation_control(simulation_control)
    , laser_parameters(p_laser_parameters)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

protected:
  const std::shared_ptr<SimulationControl> simulation_control;
  std::shared_ptr<Parameters::Laser<dim>>  laser_parameters;
};


/**
 * @brief Class that assembles the radiation sink for the heat
 * transfer solver at the free surface (air/metal interface) when VOF and the
 * laser are active. The phase gradient of the VOF solver is used to transform
 * radiative boundary condition at free surface into a volumetric sink at the
 * air/metal interface: "Tao Yu, Jidong Zhao,Semi-coupled resolved CFD–DEM
 * simulation of powder-based selective laser melting for additive
 * manufacturing, Computer Methods in Applied Mechanics and Engineering, Volume
 * 377, 2021, 113707, ISSN 0045-7825,
 * https://doi.org/10.1016/j.cma.2021.113707."
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 * @param p_laser_parameters Shared pointer of the laser parameters
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerFreeSurfaceRadiationVOF
  : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerFreeSurfaceRadiationVOF(
    const std::shared_ptr<SimulationControl> &simulation_control,
    std::shared_ptr<Parameters::Laser<dim>>   p_laser_parameters)
    : simulation_control(simulation_control)
    , laser_parameters(p_laser_parameters)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

protected:
  const std::shared_ptr<SimulationControl> simulation_control;
  std::shared_ptr<Parameters::Laser<dim>>  laser_parameters;
};


/**
 * @brief Class that assembles the evaporation sink for the heat
 * transfer solver at the free surface (air/metal interface) when VOF is
 * enabled.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation.
 * @param p_evaporation Struct that holds all evaporation model
 * parameters.
 * @ingroup assemblers.
 */
template <int dim>
class HeatTransferAssemblerVOFEvaporation
  : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerVOFEvaporation(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::Evaporation            &p_evaporation)
    : simulation_control(simulation_control)
  {
    this->evaporation_model = EvaporationModel::model_cast(p_evaporation);
  }

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

private:
  const std::shared_ptr<SimulationControl> simulation_control;

  // Evaporation model
  std::shared_ptr<EvaporationModel> evaporation_model;
};

#endif
