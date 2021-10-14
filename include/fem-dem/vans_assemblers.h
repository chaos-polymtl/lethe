/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 * ---------------------------------------------------------------------*/


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
 * of the assemblers for the particle_fluid interqactions of the
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
   * @brief calculate_particle_fluid_interactions calculted the solid_fluid interactions
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
 * @brief Class that assembles the core of the Volume Averaged Navier-Stokes equations.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSVansAssemblerCore : public NavierStokesAssemblerBase<dim>
{
public:
  GLSVansAssemblerCore(std::shared_ptr<SimulationControl> simulation_control,
                       Parameters::PhysicalProperties     physical_properties,
                       Parameters::CFDDEM                 cfd_dem)
    : simulation_control(simulation_control)
    , physical_properties(physical_properties)
    , cfd_dem(cfd_dem)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)Particles::ParticleHandler
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  const bool SUPG = true;

  std::shared_ptr<SimulationControl> simulation_control;
  Parameters::PhysicalProperties     physical_properties;
  Parameters::CFDDEM                 cfd_dem;
};

/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the Navier-Stokes equation with
 * free surface using VOF modeling.. For example, if a BDF1 scheme is
 * chosen, the following is assembled
 * $$\frac{(\rho \mathbf{u})^{t+\Delta t}-(\rho \mathbf{u})^{t}{\Delta t}
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class GLSVansAssemblerBDF : public NavierStokesAssemblerBase<dim>
{
public:
  GLSVansAssemblerBDF(std::shared_ptr<SimulationControl> simulation_control,
                      Parameters::CFDDEM                 cfd_dem)
    : simulation_control(simulation_control)
    , cfd_dem(cfd_dem)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
  Parameters::CFDDEM                 cfd_dem;
};

/**
 * @brief Class that assembles the drag force using DiFelice model for the
 * VANS equations where the drag coefficient c_d = pow((0.63 + 4.8 / sqrt(re)),
 2) * pow(cell_void_fraction,
                -(3.7 - 0.65 * exp(-pow((1.5 - log10(re)), 2) / 2)))
 *  and the momentum exchange coefficient
 *  beta =(0.5 * c_d * M_PI *
         pow(particle_properties[DEM::PropertiesIndex::dp], 2) / 4) *
        relative_velocity.norm()
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class GLSVansAssemblerDiFelice : public ParticleFluidAssemblerBase<dim>
{
public:
  GLSVansAssemblerDiFelice(
    std::shared_ptr<SimulationControl> simulation_control,
    Parameters::PhysicalProperties     physical_properties)

    : simulation_control(simulation_control)
    , physical_properties(physical_properties)

  {}

  /**
   * @brief calculate_particle_fluid_interactions calculted the solid_fluid interactions
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
  Parameters::PhysicalProperties     physical_properties;
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
         pow(particle_properties[DEM::PropertiesIndex::dp], 2) / 4) *
        relative_velocity.norm()
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class GLSVansAssemblerRong : public ParticleFluidAssemblerBase<dim>
{
public:
  GLSVansAssemblerRong(std::shared_ptr<SimulationControl> simulation_control,
                       Parameters::PhysicalProperties     physical_properties)

    : simulation_control(simulation_control)
    , physical_properties(physical_properties)

  {}

  /**
   * @brief calculate_particle_fluid_interactions calculted the solid_fluid interactions
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
  Parameters::PhysicalProperties     physical_properties;
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
class GLSVansAssemblerBuoyancy : public ParticleFluidAssemblerBase<dim>
{
public:
  GLSVansAssemblerBuoyancy(
    std::shared_ptr<SimulationControl> simulation_control,
    Parameters::PhysicalProperties     physical_properties,
    Parameters::Lagrangian::LagrangianPhysicalProperties<dim>
      lagrangian_physical_properties)

    : simulation_control(simulation_control)
    , physical_properties(physical_properties)
    , lagrangian_physical_properties(lagrangian_physical_properties)

  {}

  /**
   * @brief calculate_particle_fluid_interactions calculted the solid_fluid interactions
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  calculate_particle_fluid_interactions(
    NavierStokesScratchData<dim> &scratch_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
  Parameters::PhysicalProperties     physical_properties;
  Parameters::Lagrangian::LagrangianPhysicalProperties<dim>
    lagrangian_physical_properties;
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
class GLSVansAssemblerFPI : public NavierStokesAssemblerBase<dim>
{
public:
  GLSVansAssemblerFPI(std::shared_ptr<SimulationControl> simulation_control,
                      Parameters::PhysicalProperties     physical_properties,
                      Parameters::CFDDEM                 cfd_dem)

    : simulation_control(simulation_control)
    , physical_properties(physical_properties)
    , cfd_dem(cfd_dem)

  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
  Parameters::PhysicalProperties     physical_properties;
  Parameters::CFDDEM                 cfd_dem;
};


#endif
