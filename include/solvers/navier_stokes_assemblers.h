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


#include <core/boundary_conditions.h>
#include <core/rheological_model.h>
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
   * @param copy_data Stores the local_rhs and local_matrix that will be
   * written into the global_rhs and global_matrix
   */

  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) = 0;


  /**
   * @brief assemble_matrix Interface for the call to rhs
   * @param scratch_data Scratch data containing the Navier-Stokes information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Stores the local_rhs and local_matrix that will be
   * written into the global_rhs and global_matrix
   */

  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) = 0;
};

/**
 * @brief Class that assembles the core of the Navier-Stokes equation.
 * This class assembles the weak form of:
 * $$\mathbf{u} \cdot \nabla \mathbf{u} - \nabla p - \nu \nabla^2 \mathbf{u}
 * =0 $$ with an SUPG and PSPG stabilization
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */


template <int dim>
class PSPGSUPGNavierStokesAssemblerCore : public NavierStokesAssemblerBase<dim>
{
public:
  PSPGSUPGNavierStokesAssemblerCore(
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
 * =0 $$ with a full GLS stabilization including the laplacian of the test
 * function.
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
 * @brief Class that assembles the core of the Navier-Stokes equation
 * using a Rheological model to predict non Newtonian behaviors
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */


template <int dim>
class GLSNavierStokesAssemblerNonNewtonianCore
  : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesAssemblerNonNewtonianCore(
    std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief Calculates an approximation of the gradient of the viscosity
   * @param velocity_gradient The velocity gradient tensor on the quadrature point
     @param velocity_hessians The velocity hessian tensor on the quadrature point
     @param non_newtonian_viscosity The viscosity at which the gradient is calculated
     @param d_gamma_dot Th difference in the shear rate magnitude to approximate the
     viscosity variation with a slight change in the shear_rate magnitude
   */
  inline Tensor<1, dim>
  get_viscosity_gradient(const Tensor<2, dim> &velocity_gradient,
                         const Tensor<3, dim> &velocity_hessians,
                         const double          shear_rate_magnitude,
                         const double          grad_viscosity_shear_rate) const
  {
    // Calculates an approximation of the shear rate magnitude gradient using
    // the derived form, since it does not change with rheological models
    Tensor<1, dim> grad_shear_rate;
    for (unsigned int d = 0; d < dim; ++d)
      {
        if (dim == 2)
          {
            for (unsigned int k = 0; k < dim; ++k)
              {
                grad_shear_rate[d] +=
                  2 * (velocity_gradient[k][k] * velocity_hessians[k][d][k]) /
                  shear_rate_magnitude;
              }
            grad_shear_rate[d] +=
              (velocity_gradient[0][1] + velocity_gradient[1][0]) *
              (velocity_hessians[0][d][1] + velocity_hessians[1][d][0]) /
              shear_rate_magnitude;
          }
        else
          {
            for (unsigned int k = 0; k < dim; ++k)
              {
                grad_shear_rate[d] +=
                  2 * (velocity_gradient[k][k] * velocity_hessians[k][d][k]) /
                    shear_rate_magnitude +
                  (velocity_gradient[(k + 1) % dim][(k + 2) % dim] +
                   velocity_gradient[(k + 2) % dim][(k + 1) % dim]) *
                    (velocity_hessians[(k + 1) % dim][d][(k + 2) % dim] +
                     velocity_hessians[(k + 2) % dim][d][(k + 1) % dim]) /
                    shear_rate_magnitude;
              }
          }
      }
    return grad_shear_rate * grad_viscosity_shear_rate;
  };

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
 * integration for the Navier Stokes equations.
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
    const double                       gamma)
    : simulation_control(simulation_control)
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
  double                             gamma;
};

template <int dim>
class GDNavierStokesAssemblerNonNewtonianCore
  : public NavierStokesAssemblerBase<dim>
{
public:
  GDNavierStokesAssemblerNonNewtonianCore(
    std::shared_ptr<SimulationControl> simulation_control,
    const double                       gamma)
    : simulation_control(simulation_control)
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
  LaplaceAssembly(std::shared_ptr<SimulationControl> simulation_control)
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
 * @brief Class that assembles a Neumann boundary condition.
 * This class assembles the weak form of: (p-mu*grad_u)*n at the boundary
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 * @param pressure_boundary_condition The boundary condition objects use to store the function.
 * @ingroup assemblers
 */

template <int dim>
class PressureBoundaryCondition : public NavierStokesAssemblerBase<dim>
{
public:
  PressureBoundaryCondition(
    std::shared_ptr<SimulationControl> simulation_control,
    const BoundaryConditions::NSBoundaryConditions<dim>
      &pressure_boundary_conditions_input)
    : simulation_control(simulation_control)
    , pressure_boundary_conditions(pressure_boundary_conditions_input)
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
  const BoundaryConditions::NSBoundaryConditions<dim>
    &pressure_boundary_conditions;
};

/**
 * @brief Class that assembles the weak formulation of a Dirichlet boundary condition using the Nitsche method.
 * This class assembles the weak form of: (u_ib-u)-(u,grad v)
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 * @param boundary_condition The boundary condition objects us to store the function.
 * @ingroup assemblers
 */

template <int dim>
class WeakDirichletBoundaryCondition : public NavierStokesAssemblerBase<dim>
{
public:
  WeakDirichletBoundaryCondition(
    std::shared_ptr<SimulationControl> simulation_control,
    const BoundaryConditions::NSBoundaryConditions<dim>
      &boundary_conditions_input)
    : simulation_control(simulation_control)
    , boundary_conditions(boundary_conditions_input)
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


  std::shared_ptr<SimulationControl>                   simulation_control;
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions;
};

/**
 * @brief Class that assembles the special case of partial slip condition using the weak formulation of a Dirichlet boundary condition using the Nitsche method.
 * This class assembles the weak form of: beta n(nu-nv) +
 * mu/boundary_layer_thickness t(tu-tv)
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 * @param boundary_condition The boundary condition objects us to store the function.
 * @ingroup assemblers
 */

template <int dim>
class PartialSlipDirichletBoundaryCondition
  : public NavierStokesAssemblerBase<dim>
{
public:
  PartialSlipDirichletBoundaryCondition(
    std::shared_ptr<SimulationControl> simulation_control,
    const BoundaryConditions::NSBoundaryConditions<dim>
      &boundary_conditions_input)
    : simulation_control(simulation_control)
    , boundary_conditions(boundary_conditions_input)
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


  std::shared_ptr<SimulationControl>                   simulation_control;
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions;
};

/**
 * @brief Class that assembles the weak formulation of an outlet boundary condition.
 * This class assembles the weak form of (nu grad(u) * grad(v) + pI - (beta *
 * u)_ * n See the paper by Arndt, Braack and Lube
 * https://www.mathsim.eu/~darndt/files/ENUMATH_2015.pdf
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 * @param boundary_condition The boundary condition objects us to store the function.
 * @ingroup assemblers
 */

template <int dim>
class OutletBoundaryCondition : public NavierStokesAssemblerBase<dim>
{
public:
  OutletBoundaryCondition(std::shared_ptr<SimulationControl> simulation_control,
                          const BoundaryConditions::NSBoundaryConditions<dim>
                            &boundary_conditions_input)
    : simulation_control(simulation_control)
    , boundary_conditions(boundary_conditions_input)
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


  std::shared_ptr<SimulationControl>                   simulation_control;
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions;
};



/**
 * @brief Class that assembles a buoyancy forcing term using the Boussinesq
 * approximation. For more information, read Chapter 10 of Transport phenomena
 * by Bird et al., or "Boussinesq approximation (buoyancy)" page on Wikipedia.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class BuoyancyAssembly : public NavierStokesAssemblerBase<dim>
{
public:
  BuoyancyAssembly(std::shared_ptr<SimulationControl> simulation_control)
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
   * @brief assemble_rhs Assembles the weak form of: $$-\mathbf{g} \times \alpha \times (T - T_0)$$
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;


  std::shared_ptr<SimulationControl> simulation_control;
};

#endif
