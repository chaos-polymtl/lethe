// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vans_assemblers_h
#define lethe_vans_assemblers_h

#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/navier_stokes_assemblers.h>
#include <solvers/navier_stokes_scratch_data.h>

#include <fem-dem/cfd_dem_simulation_parameters.h>

#include <deal.II/particles/particle_handler.h>


/**
 * @file vans_assemblers.h
 * @brief Assemblers for the Volume-Averaged Navier-Stokes (VANS) equations.
 *
 * This file contains the assemblers used to build the matrix and right-hand
 * side contributions for the VANS equations in CFD-DEM simulations. It
 * includes the core equation assemblers (Model A and Model B), various drag
 * models (Di Felice, Rong, Dallavalle, Koch-Hill, Beetstra, Gidaspow), lift
 * force models (Saffman-Mei, Magnus), torque models (viscous, vortical),
 * buoyancy, pressure and shear force assemblers, and the fluid-particle
 * interaction assemblers.
 */

/**
 * @brief Interface for all assemblers of particle-fluid interactions in the
 * VANS equations.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class ParticleFluidAssemblerBase
{
public:
  /**
   * @brief Calculate the particle-fluid interaction for all particles in a
   * cell.
   *
   * The forces calculated on the particles are stored within the
   * @p fem_force field. The overall \f$ \beta \f$ coefficient for the cell,
   * which is defined as the drag coefficient divided by the cell volume, is
   * also calculated within this function.
   *
   * @param[in,out] scratch_data Scratch data containing the Navier-Stokes
   * information. The scratch data must have been re-initialized before calling
   * this function.
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) = 0;

  /**
   * @brief Default destructor.
   */
  virtual ~ParticleFluidAssemblerBase() = default;
};


/**
 * @brief Assembler for the core of Model B of the Volume-Averaged
 * Navier-Stokes (VANS) equations.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerCoreModelB : public NavierStokesAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerCoreModelB object.
   *
   * @param[in] simulation_control Shared pointer to the simulation control
   * object used to retrieve time-stepping information.
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerCoreModelB(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::CFDDEM                 &cfd_dem)
    : simulation_control(simulation_control)
    , cfd_dem(cfd_dem)
  {}

  /**
   * @brief Assemble the matrix.
   *
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Assemble the right-hand side.
   *
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /// Flag indicating whether SUPG stabilization is used.
  const bool SUPG = true;

  /// Shared pointer to the simulation control object.
  const std::shared_ptr<SimulationControl> simulation_control;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Assembler for the core of Model A of the Volume-Averaged
 * Navier-Stokes (VANS) equations.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerCoreModelA : public NavierStokesAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerCoreModelA object.
   *
   * @param[in] simulation_control Shared pointer to the simulation control
   * object used to retrieve time-stepping information.
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerCoreModelA(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::CFDDEM                 &cfd_dem)
    : simulation_control(simulation_control)
    , cfd_dem(cfd_dem)
  {}

  /**
   * @brief Assemble the matrix.
   *
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Assemble the right-hand side.
   *
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /// Flag indicating whether SUPG stabilization is used.
  const bool SUPG = true;

  /// Shared pointer to the simulation control object.
  const std::shared_ptr<SimulationControl> simulation_control;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Assembler for the transient term arising from BDF time integration
 * for the VANS equations.
 *
 * For example, if a BDF1 scheme is chosen, the following is assembled:
 * \f[
 * \frac{(\rho \mathbf{u})^{t+\Delta t} - (\rho \mathbf{u})^{t}}{\Delta t}
 * \f]
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerBDF : public NavierStokesAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerBDF object.
   *
   * @param[in] simulation_control Shared pointer to the simulation control
   * object used to retrieve time-stepping information.
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerBDF(const std::shared_ptr<SimulationControl> &simulation_control,
                   const Parameters::CFDDEM                 &cfd_dem)
    : simulation_control(simulation_control)
    , cfd_dem(cfd_dem)
  {}

  /**
   * @brief Assemble the matrix.
   *
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Assemble the right-hand side.
   *
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /// Shared pointer to the simulation control object.
  const std::shared_ptr<SimulationControl> simulation_control;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};


/**
 * @brief Assembler for the drag force using the Di Felice model.
 *
 * The drag coefficient is computed as:
 * \f[
 * C_d = \left(0.63 + \frac{4.8}{\sqrt{\text{Re}_p}}\right)^2
 *       \varepsilon^{-(3.7 - 0.65
 *       \exp(-(1.5 - \log_{10} \text{Re}_p)^2 / 2))}
 * \f]
 *
 * and the momentum exchange coefficient for each particle is:
 * \f[
 * \beta = \frac{1}{2} C_d \frac{\pi d_p^2}{4} \|\mathbf{u}_r\|
 * \f]
 *
 * where \f$ \varepsilon \f$ is the cell void fraction, \f$ \text{Re}_p \f$ is
 * the particle Reynolds number, \f$ d_p \f$ is the particle diameter, and
 * \f$ \mathbf{u}_r \f$ is the relative velocity between the fluid and the
 * particle.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerDiFelice : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerDiFelice object.
   *
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerDiFelice(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief Calculate the particle-fluid interactions resulting from drag, using
   * the Di Felice model.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Assembler for the drag force using the Rong model.
 *
 * The drag coefficient is computed as:
 * \f[
 * C_d = \left(0.63 + \frac{4.8}{\sqrt{\text{Re}_p}}\right)^2
 *       \varepsilon^{-(2.65(\varepsilon + 1) -
 *       (5.3 - 3.5\varepsilon)\varepsilon^2
 *       \exp(-(1.5 - \log_{10}\text{Re}_p)^2 / 2))}
 * \f]
 *
 * and the momentum exchange coefficient for each particle is:
 * \f[
 * \beta = \frac{1}{2} C_d \frac{\pi d_p^2}{4} \|\mathbf{u}_r\|
 * \f]
 *
 * where \f$ \varepsilon \f$ is the cell void fraction, \f$ \text{Re}_p \f$ is
 * the particle Reynolds number, \f$ d_p \f$ is the particle diameter, and
 * \f$ \mathbf{u}_r \f$ is the relative velocity between the fluid and the
 * particle.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerRong : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerRong object.
   *
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerRong(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
  * @brief Calculate the particle-fluid interactions resulting from drag, using
   * the Rong model.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Assembler for the drag force using the Dallavalle model.
 *
 * The drag coefficient is computed as:
 * \f[
 * C_d = \left(0.63 + \frac{4.8}{\sqrt{\text{Re}_p}}\right)^2
 * \f]
 *
 * and the momentum exchange coefficient for each particle is:
 * \f[
 * \beta = \frac{1}{2} C_d \frac{\pi d_p^2}{4} \|\mathbf{u}_r\|
 * \f]
 *
 * where \f$ \text{Re}_p \f$ is the particle Reynolds number, \f$ d_p \f$ is
 * the particle diameter, and \f$ \mathbf{u}_r \f$ is the relative velocity
 * between the fluid and the particle.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerDallavalle : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerDallavalle object.
   *
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerDallavalle(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief Calculate the particle-fluid interactions resulting from drag, using
   * the Dallavalle model.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Assembler for the drag force using the Koch and Hill drag model.
 *
 * The momentum exchange coefficient is computed as:
 * \f[
 * \beta = \frac{18 \mu \varepsilon^2 (1 - \varepsilon)}{d_p^2}
 *         \left(F_0 + \frac{1}{2} F_3 \varepsilon \text{Re}_p\right)
 *         \frac{V_p}{1 - \varepsilon}
 * \f]
 *
 * where \f$ F_0 \f$ and \f$ F_3 \f$ are defined as:
 *
 * For \f$ (1 - \varepsilon) < 0.4 \f$:
 * \f[
 * F_0 = \frac{1 + 3\sqrt{(1-\varepsilon)/2}
 *        + \frac{135}{64}(1-\varepsilon)\ln(1-\varepsilon)
 *        + 16.14(1-\varepsilon)}
 *       {1 + 0.681(1-\varepsilon)
 *        - 8.48(1-\varepsilon)^2
 *        + 8.14(1-\varepsilon)^3}
 * \f]
 *
 * For \f$ (1 - \varepsilon) \geq 0.4 \f$:
 * \f[
 * F_0 = \frac{10(1-\varepsilon)}{\varepsilon^3}
 * \f]
 *
 * and:
 * \f[
 * F_3 = 0.0673 + 0.212(1-\varepsilon) + \frac{0.0232}{\varepsilon^5}
 * \f]
 *
 * where \f$ \varepsilon \f$ is the cell void fraction, \f$ \mu \f$ is the
 * dynamic viscosity, \f$ d_p \f$ is the particle diameter, \f$ V_p \f$ is the
 * particle volume, and \f$ \text{Re}_p \f$ is the particle Reynolds number.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerKochHill : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerKochHill object.
   *
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerKochHill(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief Calculate the particle-fluid interactions resulting from drag, using
   * the Koch-Hill model.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Assembler for the drag force using the Beetstra model.
 *
 * The normalized drag force is computed as:
 * \f[
 * F_{\text{norm}} = \frac{10(1-\varepsilon)}{\varepsilon^2}
 *   + \varepsilon^2 \left(1 + 1.5\sqrt{1-\varepsilon}\right)
 *   + \frac{0.413\,\text{Re}_p}{24\varepsilon^2}
 *     \frac{\frac{1}{\varepsilon}
 *           + 3(1-\varepsilon)\varepsilon
 *           + 8.4\,\text{Re}_p^{-0.343}}
 *          {1 + 10^{3(1-\varepsilon)}
 *           \text{Re}_p^{-0.5(1 + 4(1-\varepsilon))}}
 * \f]
 *
 * The drag coefficient is then:
 * \f[
 * C_d = \frac{24\,F_{\text{norm}}}{\text{Re}_p}
 * \f]
 *
 * where \f$ \varepsilon \f$ is the cell void fraction and \f$ \text{Re}_p \f$
 * is the particle Reynolds number.
 *
 * Reference: Beetstra, R., Martin Anton van der Hoef, and J. A. M. Kuipers.
 * "Drag force of intermediate Reynolds number flow past mono‐and bidisperse
 * arrays of spheres." AIChE journal 53.2 (2007): 489-501.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerBeetstra : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerBeetstra object.
   *
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerBeetstra(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief Calculate the particle-fluid interactions resulting from drag, using
   * the Beetstra model..
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Assembler for the drag force using the Gidaspow model.
 *
 * The momentum transfer coefficient is calculated following Marchelli et al.
 * (2020). For \f$ \varepsilon > 0.8 \f$:
 * \f[
 * \beta = 18 \varepsilon^{-3.65}
 *         \left(1 + 0.15\,\text{Re}_p^{0.687}\right)
 *         \frac{m_p \mu}{d_p^2 \rho_p}
 * \f]
 *
 * For \f$ \varepsilon \leq 0.8 \f$ (assuming spherical particles):
 * \f[
 * \beta = \left(\frac{150(1-\varepsilon)}{\varepsilon^2}
 *         + \frac{1.75\,\text{Re}_p}{\varepsilon^2}\right)
 *         \frac{m_p \mu}{d_p^2 \rho_p}
 * \f]
 *
 * where \f$ \varepsilon \f$ is the cell void fraction, \f$ \text{Re}_p \f$ is
 * the particle Reynolds number, \f$ m_p \f$ is the particle mass, \f$ \mu \f$
 * is the dynamic viscosity, \f$ d_p \f$ is the particle diameter, and
 * \f$ \rho_p \f$ is the particle density.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerGidaspow : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerGidaspow object.
   *
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerGidaspow(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief Calculate the particle-fluid interactions resulting from drag, using
   * the Gidaspow model.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Assembler for the lift force using the Saffman-Mei model.
 *
 * This implementation follows the formulation in the book "Multiphase Flows
 * with Droplets and Particles" by Crowe et al. (2011) and the brief
 * communication article "An approximate expression for the shear lift force
 * on a spherical particle at finite Reynolds number" by Mei (1992).
 *
 * The lift force is computed as:
 * \f[
 * \mathbf{F}_L = 1.61 \, C_s \, d_p^2 \, \rho_f
 *   \frac{\sqrt{\nu}}{\sqrt{\|\boldsymbol{\omega}_f\|}}
 *   \left(\mathbf{u}_r \times \boldsymbol{\omega}_f\right)
 * \f]
 *
 * where the parameter \f$ \alpha \f$ is defined as:
 * \f[
 * \alpha = \frac{d_p}{2 \|\mathbf{u}_r\|} \|\boldsymbol{\omega}_f\|
 * \f]
 *
 * and the Saffman-Mei correction coefficient \f$ C_s \f$ is:
 *
 * For \f$ \text{Re}_p \leq 40 \f$:
 * \f[
 * C_s = (1 - 0.3314\sqrt{\alpha}) \exp(-0.1\,\text{Re}_p)
 *       + 0.3314\sqrt{\alpha}
 * \f]
 *
 * For \f$ \text{Re}_p > 40 \f$:
 * \f[
 * C_s = 0.0524\sqrt{\alpha \, \text{Re}_p}
 * \f]
 *
 * where \f$ d_p \f$ is the particle diameter, \f$ \rho_f \f$ is the fluid
 * density, \f$ \nu \f$ is the kinematic viscosity, \f$ \boldsymbol{\omega}_f
 * \f$ is the fluid vorticity, \f$ \mathbf{u}_r \f$ is the relative velocity
 * between the fluid and the particle, and \f$ \text{Re}_p \f$ is the particle
 * Reynolds number.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerSaffmanMei : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Default constructor.
   */
  VANSAssemblerSaffmanMei()
  {}

  /**
   * @brief Calculate the Saffman-Mei lift force on particles.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;
};

/**
 * @brief Assembler for the lift force using the Magnus model.
 *
 * This implementation follows the formulation in the book "Multiphase Flows
 * with Droplets and Particles" by Crowe et al. (2011).
 *
 * The Magnus lift force is computed as:
 * \f[
 * \mathbf{F}_L = \frac{\pi}{8} d_p^2 \, C_m \, \rho_f \,
 *   \|\mathbf{u}_r\| \left(\hat{\boldsymbol{\omega}}_p \times
 *   \mathbf{u}_r\right)
 * \f]
 *
 * where \f$ \hat{\boldsymbol{\omega}}_p =
 * \boldsymbol{\omega}_p / \|\boldsymbol{\omega}_p\| \f$ is the unit rotational
 * vector. The spin parameter is defined as:
 * \f[
 * \Omega_s = \frac{d_p \|\boldsymbol{\omega}_p\|}{2 \|\mathbf{u}_r\|}
 * \f]
 *
 * The Magnus lift coefficient \f$ C_m \f$ is computed as:
 *
 * For \f$ 1 < \Omega_s < 6 \f$ and \f$ 10 < \text{Re}_p < 140 \f$
 * (Oesterle and Dinh, 1998):
 * \f[
 * C_m = 0.45 + (2\Omega_s - 0.45)
 *       \exp\!\left(-0.075 \, \Omega_s^{0.4} \, \text{Re}_p^{0.7}\right)
 * \f]
 *
 * Otherwise:
 * \f[
 * C_m = 2 \Omega_s
 * \f]
 *
 * where \f$ d_p \f$ is the particle diameter, \f$ \rho_f \f$ is the fluid
 * density, \f$ \mathbf{u}_r \f$ is the relative velocity between the fluid and
 * the particle, \f$ \boldsymbol{\omega}_p \f$ is the particle angular velocity,
 * and \f$ \text{Re}_p \f$ is the particle Reynolds number.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerMagnus : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Default constructor.
   */
  VANSAssemblerMagnus()
  {}

  /**
   * @brief Calculate the Magnus lift force on particles.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;
};

/**
 * @brief Assembler for the viscous dissipative torque due to particle rotation
 * as defined by Derksen (2004).
 *
 * The viscous torque due to particle rotation is:
 * \f[
 * \mathbf{M}_{\text{viscous}} = -\frac{\pi}{2} d_p^3 \mu
 *   \boldsymbol{\omega}_p
 * \f]
 *
 * where \f$ d_p \f$ is the particle diameter, \f$ \mu = \nu \rho_f \f$ is
 * the dynamic viscosity, and \f$ \boldsymbol{\omega}_p \f$ is the particle
 * angular velocity. This torque opposes the particle rotation.
 *
 * @note The complete model described by Derksen is composed of
 * VANSAssemblerViscousTorque + VANSAssemblerVorticalTorque.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerViscousTorque : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Default constructor.
   */
  VANSAssemblerViscousTorque()
  {}

  /**
   * @brief Calculate the viscous torque dissipation on particles.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;
};

/**
 * @brief Assembler for the viscous dissipative torque due to fluid vorticity
 * as defined by Derksen (2004).
 *
 * The viscous torque due to the fluid vorticity is:
 * \f[
 * \mathbf{M}_{\text{vorticity}} = \frac{\pi}{2} d_p^3 \mu
 *   \boldsymbol{\omega}_f
 * \f]
 *
 * where \f$ d_p \f$ is the particle diameter, \f$ \mu = \nu \rho_f \f$ is
 * the dynamic viscosity, and \f$ \boldsymbol{\omega}_f \f$ is the fluid
 * vorticity. This torque drives the particle rotation towards the local fluid
 * vorticity.
 *
 * @note The complete model described by Derksen is composed of
 * VANSAssemblerViscousTorque + VANSAssemblerVorticalTorque.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerVorticalTorque : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Default constructor.
   */
  VANSAssemblerVorticalTorque()
  {}

  /**
   * @brief Calculate the vortical torque on particles.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;
};

/**
 * @brief Assembler for the buoyancy force on particles in the VANS equations.
 *
 * The buoyancy force applied to each particle is computed as:
 * \f[
 * \mathbf{F}_b = -\rho_f \mathbf{g} \frac{\pi}{6} d_p^3
 * \f]
 *
 * which is equivalent to \f$ -\rho_f \mathbf{g} V_p \f$, where \f$ V_p =
 * \frac{4}{3}\pi\left(\frac{d_p}{2}\right)^3 \f$ is the particle volume,
 * \f$ \rho_f \f$ is the fluid density, \f$ \mathbf{g} \f$ is the
 * gravitational acceleration, and \f$ d_p \f$ is the particle diameter.
 *
 * @note The buoyancy force is stored as a one-way coupling force on the
 * particle and does not generate a reaction force on the fluid.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerBuoyancy : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerBuoyancy object.
   *
   * @param[in] p_gravity The gravity acceleration vector applied to the
   * particles.
   */
  VANSAssemblerBuoyancy(const Tensor<1, 3> &p_gravity)
    : gravity(p_gravity)
  {}

  /**
   * @brief Calculate the buoyancy force on particles.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  /// Gravity acceleration applied to the particles.
  const Tensor<1, 3> gravity;
};

/**
 * @brief Assembler for the pressure force on particles.
 *
 * The undisturbed pressure force applied to each particle is computed as:
 * \f[
 * \mathbf{F}_p = -\rho_f V_p \nabla p
 * \f]
 *
 * where \f$ \rho_f \f$ is the fluid density, \f$ V_p = \frac{\pi}{6} d_p^3
 * \f$ is the particle volume, \f$ d_p \f$ is the particle diameter, and
 * \f$ \nabla p \f$ is the fluid pressure gradient evaluated at the particle
 * location.
 *
 * @note When using VANS Model A, the pressure force is applied only to the
 * particle (one-way coupling). When using Model B, the pressure force is
 * also applied back on the fluid as a reaction force (two-way coupling).
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerPressureForce : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerPressureForce object.
   *
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerPressureForce(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief Calculate the pressure force on particles.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Assembler for the shear force on particles.
 *
 * The undisturbed shear (viscous stress) force applied to each particle is
 * computed as:
 * \f[
 * \mathbf{F}_s = -V_p \mu \nabla^2 \mathbf{u}
 * \f]
 *
 * where \f$ V_p = \frac{\pi}{6} d_p^3 \f$ is the particle volume, \f$ \mu =
 * \nu \rho_f \f$ is the dynamic viscosity, \f$ d_p \f$ is the particle
 * diameter, and \f$ \nabla^2 \mathbf{u} \f$ is the Laplacian of the fluid
 * velocity evaluated at the particle location.
 *
 * @note When using VANS Model A, the shear force is applied only to the
 * particle (one-way coupling). When using Model B, the shear force is also
 * applied back on the fluid as a reaction force (two-way coupling).
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerShearForce : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerShearForce object.
   *
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerShearForce(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief Calculate the shear force on particles.
   *
   * @param[in,out] scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Assembler for the fluid-particle interactions (FPI) in the VANS
 * equations, such as the drag force.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerFPI : public NavierStokesAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerFPI object.
   *
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerFPI(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief Assemble the matrix.
   *
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Assemble the right-hand side.
   *
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Assembler for the fluid-particle interactions (FPI) in the VANS
 * equations when the forces are projected from the particles to the fluid grid.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerFPIProjection : public NavierStokesAssemblerBase<dim>
{
public:
  /**
   * @brief Construct a VANSAssemblerFPIProjection object.
   *
   * @param[in] cfd_dem CFD-DEM simulation parameters.
   */
  VANSAssemblerFPIProjection(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief Assemble the matrix.
   *
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Assemble the right-hand side.
   *
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /// CFD-DEM simulation parameters.
  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Calculate the grad-div stabilization constant \f$ \gamma \f$ for the
 * VANS equations.
 *
 * The stabilization constant is computed as:
 * \f[
 * \gamma = \nu + c^* \|\mathbf{u}\|
 * \f]
 *
 * where \f$ \nu \f$ is the kinematic viscosity, \f$ c^* \f$ is a scaling
 * constant with units of length, and \f$ \|\mathbf{u}\| \f$ is the magnitude
 * of the velocity at the quadrature point.
 *
 * @param[in] velocity Magnitude of the velocity at the quadrature point.
 * @param[in] kinematic_viscosity Kinematic viscosity of the fluid.
 * @param[in] h The element size (unused).
 * @param[in] c_star Scaling constant with units of length.
 *
 * @return The grad-div stabilization constant \f$ \gamma \f$.
 */
inline double
calculate_gamma(double velocity,
                double kinematic_viscosity,
                double /*h*/,
                double c_star)
{
  return kinematic_viscosity + c_star * velocity;
}
#endif
