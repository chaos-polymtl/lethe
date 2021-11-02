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
#include <solvers/navier_stokes_scratch_data.h>

#ifndef lethe_navier_stokes_assemblers_h
#  define lethe_navier_stokes_assemblers_h

/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for the Navier-Stokes Equations
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class NavierStokesAssemblerBase
{
public:
  /**
   * @brief assemble_matrix Interface for the call to matrix assembly
   * @param scratch_data Scratch data containing the Navier-Stokes information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) = 0;


  /**
   * @brief assemble_matrix Interface for the call to rhs
   * @param scratch_data Scratch data containing the Navier-Stokes information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) = 0;
};

/**
 * @brief Class that assembles the core of the Navier-Stokes equation.
 * This class assembles the weak form of:
 * $$\mathbf{u} \cdot \nabla \mathbf{u} - \nabla p - \nu \nabla^2 \mathbf{u}
 * =0 $$ with an SUPG and PSPG stabilziation
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */


template <int dim>
class GLSNavierStokesAssemblerCore : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesAssemblerCore(
    std::shared_ptr<SimulationControl> simulation_control,
    Parameters::PhysicalProperties     physical_properties)
    : simulation_control(simulation_control)
    , physical_properties(physical_properties)
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

  /**
   * Enables SUPG stabilization for the Navier-Stokes formulation.
   * We have not found any scenarios where it is relevant not to use SUPG
   * stabilization yet.
   */
  const bool SUPG = true;

  std::shared_ptr<SimulationControl> simulation_control;
  Parameters::PhysicalProperties     physical_properties;
};


/**
 * @brief Class that assembles the coriolis and the centrifugal
 * if a simulation is carried out in a rotating frame of reference.
 * This class assembles the following term:
 * $$2\mathbf{\omega} \times \mathbf{u} + \mathbf{\omega}\times (\mathbf{\omega}
 * \times \mathbf{r})$$ Where $\mathbf{\omega}$ is the rotation vector of the
 * frame of reference and $\mathbf{r}$ is the position vector (e.g. a vector
 * between a point on the rotation axis and the gauss point) By default, it is
 * assumed that the rotation vector passes through the point (0,0) in 2D or
 * (0,0,0) in 3D
 *
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class GLSNavierStokesAssemblerSRF : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesAssemblerSRF(Parameters::VelocitySource velocity_sources)
    : velocity_sources(velocity_sources)
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

  Parameters::VelocitySource velocity_sources;
};


/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the Navier Stokes equations. For example, if a BDF1 scheme is
 * chosen, the following is assembled
 * $$\frac{\mathbf{u}^{t+\Delta t}-\mathbf{u}^{t}{\Delta t}
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSNavierStokesAssemblerBDF : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesAssemblerBDF(
    std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
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
};

/**
 * @brief Class that assembles the transient time arising from SDIRK time
 * integration for the Navier Stokes equatios.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSNavierStokesAssemblerSDIRK : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesAssemblerSDIRK(
    std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
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
};


/**
 * @brief Class that assembles the core of the Navier-Stokes equation.
 * This class assembles the weak form of:
 * $$\mathbf{u} \cdot \nabla \mathbf{u} - \nabla p - \nu \nabla^2 \mathbf{u}
 * =0 $$ with a grad-div stabilization
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */


template <int dim>
class GDNavierStokesAssemblerCore : public NavierStokesAssemblerBase<dim>
{
public:
  GDNavierStokesAssemblerCore(
    std::shared_ptr<SimulationControl> simulation_control,
    Parameters::PhysicalProperties     physical_properties,
    const double                       gamma)
    : simulation_control(simulation_control)
    , physical_properties(physical_properties)
    , gamma(gamma)
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
  double                             gamma;
};


/**
 * @brief Class that assembles a Poisson problem for all velocity components and pressure variables.
 * This class assembles the weak form of: d^2 U/dx^2=0 and  d^2 P/dx^2=0
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class LaplaceAssembly : public NavierStokesAssemblerBase<dim>
{
public:
  LaplaceAssembly(std::shared_ptr<SimulationControl> simulation_control,
                  Parameters::PhysicalProperties     physical_properties)
    : simulation_control(simulation_control)
    , physical_properties(physical_properties)
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
};


/**
 * @brief Class that assembles the Buoyancy using Boussinesq approximation
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class BuoyancyAssembly : public NavierStokesAssemblerBase<dim>
{
public:
  BuoyancyAssembly(std::shared_ptr<SimulationControl> simulation_control,
                   Parameters::PhysicalProperties     physical_properties,
                   Parameters::VelocitySource         velocity_sources)
    : simulation_control(simulation_control)
    , physical_properties(physical_properties)
    , velocity_sources(velocity_sources)
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
  Parameters::VelocitySource         velocity_sources;
};


#endif
