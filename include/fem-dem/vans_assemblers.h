// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/navier_stokes_assemblers.h>
#include <solvers/navier_stokes_scratch_data.h>

#include <fem-dem/cfd_dem_simulation_parameters.h>

#include <deal.II/particles/particle_handler.h>

#ifndef lethe_vans_assemblers_h
#  define lethe_vans_assemblers_h

/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for the particle_fluid interactions of the
 * VANS equations
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class ParticleFluidAssemblerBase
{
public:
  /**
   * @brief calculate_particle_fluid_interactions calculated the solid_fluid
   * interactions
   * @param scratch_data Scratch data containing the Navier-Stokes information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) = 0;
};


/**
 * @brief Class that assembles the core of Model B of the Volume Averaged Navier-Stokes equations.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerCoreModelB : public NavierStokesAssemblerBase<dim>
{
public:
  VANSAssemblerCoreModelB(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::CFDDEM                 &cfd_dem)
    : simulation_control(simulation_control)
    , cfd_dem(cfd_dem)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)Particles::ParticleHandler
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  const bool SUPG = true;

  const std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::CFDDEM                 cfd_dem;
};

/**
 * @brief Class that assembles the core of Model A of the Volume Averaged Navier-Stokes equations.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class VANSAssemblerCoreModelA : public NavierStokesAssemblerBase<dim>
{
public:
  VANSAssemblerCoreModelA(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::CFDDEM                 &cfd_dem)
    : simulation_control(simulation_control)
    , cfd_dem(cfd_dem)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)Particles::ParticleHandler
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  const bool SUPG = true;

  const std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::CFDDEM                 cfd_dem;
};

/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the Navier-Stokes equation with
 * free surface using VOF modeling.. For example, if a BDF1 scheme is
 * chosen, the following is assembled
 * \f$\frac{(\rho \mathbf{u})^{t+\Delta t}-(\rho \mathbf{u})^{t}}{\Delta t}\f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerBDF : public NavierStokesAssemblerBase<dim>
{
public:
  VANSAssemblerBDF(const std::shared_ptr<SimulationControl> &simulation_control,
                   const Parameters::CFDDEM                 &cfd_dem)
    : simulation_control(simulation_control)
    , cfd_dem(cfd_dem)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::CFDDEM                 cfd_dem;
};


/**
 * @brief Class that assembles the drag force using DiFelice model for the
 * VANS equations where the drag coefficient c_d = pow((0.63 + 4.8 / sqrt(re)),
 2) * pow(cell_void_fraction,
                -(3.7 - 0.65 * exp(-pow((1.5 - log10(re)), 2) / 2)))
 *  and the momentum exchange coefficient
 *  beta =(0.5 * c_d * M_PI *
         pow(particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp],
 2) / 4) * relative_velocity.norm()
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerDiFelice : public ParticleFluidAssemblerBase<dim>
{
public:
  VANSAssemblerDiFelice(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief calculate_particle_fluid_interactions calculated the solid_fluid interactions
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Class that assembles the drag force using Rong model for the
 * VANS equations where the drag coefficient c_d =
        pow((0.63 + 4.8 / sqrt(re)), 2) *
        pow(cell_void_fraction,
            -(2.65 * (cell_void_fraction + 1) -
              (5.3 - (3.5 * cell_void_fraction)) * pow(cell_void_fraction, 2) *
                exp(-pow(1.5 - log10(re), 2) / 2)))
 * and the momentum exchange coefficient
 *  beta =(0.5 * c_d * M_PI *
         pow(particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp],
 2) / 4) * relative_velocity.norm()
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerRong : public ParticleFluidAssemblerBase<dim>
{
public:
  VANSAssemblerRong(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief calculate_particle_fluid_interactions calculated the solid_fluid interactions
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Class that assembles the drag force using Dallavalle model for the
 * VANS equations where the drag coefficient c_d =
        pow((0.63 + 4.8 / sqrt(re)), 2)
 * and the momentum exchange coefficient
 *  beta =(0.5 * c_d * M_PI *
         pow(particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp],
 2) / 4) * relative_velocity.norm()
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerDallavalle : public ParticleFluidAssemblerBase<dim>
{
public:
  VANSAssemblerDallavalle(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief calculate_particle_fluid_interactions calculted the solid_fluid interactions
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Class that assembles the drag force using the Koch and Hill drag model for the
 * VANS equations where the momentum exchange coefficient
 *  beta =   ((18 * mu * cell_void_fraction^2 *
          (1 - cell_void_fraction)) / pow(dp, 2)) *
        (f0 + 0.5 * f3 * cell_void_fraction * re) *
        Vp /(1 - cell_void_fraction) where f0 and f3 are functions given by:
      if ((1 - cell_void_fraction) < 0.4)
        {
          f0 = (1 + 3 * sqrt((1 - cell_void_fraction) / 2) +
                (135.0 / 64) * (1 - cell_void_fraction) *
                  log(1 - cell_void_fraction) +
                16.14 * (1 - cell_void_fraction)) /
               (1 + 0.681 * (1 - cell_void_fraction) -
                8.48 * (1 - cell_void_fraction)^2 +
                8.14 * (1 - cell_void_fraction)^3);
        }
      else if ((1 - cell_void_fraction) >= 0.4)
        {
          f0 = 10 * (1 - cell_void_fraction) / cell_void_fraction^3;
        }

      f3 = 0.0673 + 0.212 * (1 - cell_void_fraction) +
           0.0232 / pow(cell_void_fraction, 5);
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerKochHill : public ParticleFluidAssemblerBase<dim>
{
public:
  VANSAssemblerKochHill(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief calculate_particle_fluid_interactions  calculates the solid-fluid interaction of the Koch-Hill drag model.
   * @param scratch_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Class that assembles the drag force using Beetstra model for the
 * VANS equations where the
 * normalized drag force = 10 * (1 - cell_void_fraction) /
          (pow(cell_void_fraction, 2)) + pow(cell_void_fraction, 2) * (1 + 1.5 *
 pow((1 - cell_void_fraction), 0.5)) + 0.413 * re / (24 *
 pow(cell_void_fraction, 2)) *
          ((1 / cell_void_fraction) + 3 * (1 - cell_void_fraction) *
 cell_void_fraction
          + 8.4 * pow(re, -0.343)) / (1 + pow(10, 3 * (1 - cell_void_fraction))
 * pow(re,
          -(1 + 4 * (1 - cell_void_fraction)) * 0.5));
 *the drag coefficient = normalized_drag_force * 24 / Re_p
 *The reference for this formulation of the Beestra model is given in the
 article *Complete liquid-solid momentum coupling for unresolved CFD-DEM
 simulations:
 https://www.sciencedirect.com/science/article/pii/S0301932220305346
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerBeetstra : public ParticleFluidAssemblerBase<dim>
{
public:
  VANSAssemblerBeetstra(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief calculate_particle_fluid_interactions calculates the solid_fluid interactions
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Class that assembles the drag force using Gidaspow model for the
 * VANS equations where the momentum_transfer_coefficient is calculated
   according to the following Marchelli et al. (2020):
      if (cell_void_fraction > 0.8)
        {
          momentum_transfer_coefficient =
            (18 * pow(cell_void_fraction, -3.65) *
             (1 + 0.15 * pow(Re_p[particle_number], 0.687))) *
            (particle_properties[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::mass]
 * viscosity /
             (pow(particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp],
 2) * particle_density));
        }
      else
        {
          // Assuming the sphericity of particles = 1
          momentum_transfer_coefficient =
            (150 * (1 - cell_void_fraction) / pow(cell_void_fraction, 2) +
             1.75 * Re_p[particle_number] / pow(cell_void_fraction, 2)) *
            (particle_properties[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::mass]
 * viscosity /
             (pow(particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp],
 2) * particle_density));
        }


 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerGidaspow : public ParticleFluidAssemblerBase<dim>
{
public:
  VANSAssemblerGidaspow(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief calculate_particle_fluid_interactions calculted the solid_fluid interactions
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  const Parameters::CFDDEM cfd_dem;
};



/**
 * @brief Class that assembles the Lift force using Saffman-Mei model
 *
 * This implementation follows the formulation in the book "Multiphase Flows
 * with Droplets and Particles" by Crowe et al. (2011) and the brief
 * communication article "An approximate expression for the shear lift force
 * on a spherical particle at finite reynolds number" by Mei (1992)
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerSaffmanMei : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Constructor. Does not do anything.
   */
  VANSAssemblerSaffmanMei()
  {}

  /**
   * @brief calculate_particle_fluid_interactions calculates the saffman force.
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;
};


/**
 * @brief Class that assembles the Lift force using Magnus model
 *
 * This implementation follows the formulation in the book "Multiphase Flows
 * with Droplets and Particles" by Crowe et al. (2011).
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerMagnus : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Constructor. Does not do anything.
   */
  VANSAssemblerMagnus()
  {}

  /**
   * @brief calculate_particle_fluid_interactions calculates the magnus force.
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;
};

/**
 * @brief Class that assembles the Viscous dissipative torque due to particle's rotation
 * as defined
 * by Derksen (2004).
 * M_viscous_rotation = pi * pow(dp, 3) * mu * (- omega_p)
 *
 * The complete model described by Derksen is composed of
 * VANSAssemblerViscousTorque + VANSAssemblerVorticalTorque
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerViscousTorque : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Constructor. Does not do anything.
   */
  VANSAssemblerViscousTorque()
  {}

  /**
   * @brief calculate_particle_fluid_interactions calculates the viscous torque dissipation
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;
};

/**
 * @brief Class that assembles the Viscous dissipative torque due to fluid vorticity as defined
 * by Derksen (2004).
 * M_viscous_vorticity = pi * pow(dp, 3) * mu * (0.5 * vorticity)
 *
 * The complete model described by Derksen is composed of
 * VANSAssemblerViscousTorque + VANSAssemblerVorticalTorque
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerVorticalTorque : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Constructor. Does not do anything.
   */
  VANSAssemblerVorticalTorque()
  {}

  /**
   * @brief calculate_particle_fluid_interactions calculates the viscous torque dissipation
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;
};

/**
 * @brief Class that assembles the Buoyancy force  for the
 * VANS equations whe F_b =  -g *
        density * (4.0 / 3) * M_PI *
        pow((dp /2.0),3)
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerBuoyancy : public ParticleFluidAssemblerBase<dim>
{
public:
  /**
   * @brief Constructor. Keeps an internal copy of the gravity.
   *
   * @param gravity the gravity applied to the particles.
   */

  VANSAssemblerBuoyancy(const Tensor<1, 3> &p_gravity)
    : gravity(p_gravity)

  {}

  /**
   * @brief calculate_particle_fluid_interactions calculates the buoyancy force.
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  /// Gravity acceleration applied to the particles
  const Tensor<1, 3> gravity;
};


/**
 * @brief Class that assembles the particle pressure force that will then
 * be added to particle_fluid_interactions.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerPressureForce : public ParticleFluidAssemblerBase<dim>
{
public:
  VANSAssemblerPressureForce(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief calculate_particle_fluid_interactions calculates the pressure force
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Class that assembles the particle shear force that will then
 * be added to particle_fluid_interactions
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerShearForce : public ParticleFluidAssemblerBase<dim>
{
public:
  VANSAssemblerShearForce(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)
  {}

  /**
   * @brief calculate_particle_fluid_interactions calculates the shear force
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Class that assembles the fluid_particle interactions (FPI) for the
 * VANS equations such as the drag force.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class VANSAssemblerFPI : public NavierStokesAssemblerBase<dim>
{
public:
  VANSAssemblerFPI(const Parameters::CFDDEM &cfd_dem)
    : cfd_dem(cfd_dem)

  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  const Parameters::CFDDEM cfd_dem;
};

/**
 * @brief Calculate gamma grad-div stabilization constant for the VANS equations
 *
 * @param[in] velocity Magnitude of the velocity at the quadrature point
 * @param[in] kinematic_viscosity
 * @param[in] h The element size
 * @param[in] c_star Scaling constant with units of length
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
