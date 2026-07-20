// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_navier_stokes_assemblers_h
#define lethe_navier_stokes_assemblers_h

#include <core/boundary_conditions.h>
#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/navier_stokes_scratch_data.h>
#include <solvers/physics_assemblers.h>

/*
 * Exceptions used to capture incoherent setup of assemblers
 */

DeclExceptionMsg(
  PhaseChangePermeabilityModelRequiresTemperature,
  "Using a phase change permeability model requires running a multiphysics simulation with the heat transfer solver enabled.");

DeclExceptionMsg(
  PhaseChangePermeabilityModelDoesNotSupportCHN,
  "The phase change permeability models do not currently have a Cahn-Hilliard implementation.");

DeclExceptionMsg(
  SingleFluidPhaseChangeDarcyModelDoesNotSupportDensityMultiplication,
  "Single-fluid simulations with Darcy phase change model does not support the multiplication by density since the kinematic pressure is used in the momentum equation.");


/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for the Navier-Stokes Equations
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
using NavierStokesAssemblerBase =
  PhysicsAssemblerBase<NavierStokesScratchData<dim>,
                       StabilizedMethodsTensorCopyData<dim>>;


/**
 * @brief Class that assembles the core of the Navier-Stokes equation.
 * According to the following weak form:
 * \f$\mathbf{u} \cdot \nabla \mathbf{u} - \nabla p - \nu \nabla^2 \mathbf{u}
 * =0 \f$ with an SUPG and PSPG stabilization
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
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
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
};


/**
 * @brief Class that assembles the core of the Navier-Stokes equation.
 * According to the following weak form:
 * \f$\mathbf{u} \cdot \nabla \mathbf{u} - \nabla p - \nu \nabla^2 \mathbf{u}
 * =0 \f$ with a full GLS stabilization including the laplacian of the test
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
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
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
};


/**
 * @brief Class that assembles the core of the Navier-Stokes equation with the
 * residual-based variational multiscale (RBVMS) stabilization of Bazilevs,
 * Calo, Cottrell, Hughes, Reali & Scovazzi, "Variational multiscale
 * residual-based turbulence modeling for large eddy simulation of
 * incompressible flows", CMAME 197 (2007) 173-201.
 *
 * On top of the SUPG/PSPG terms it adds the metric-based stabilization
 * parameters tau_M (eq. 64) and tau_C (eq. 65) built from the element metric
 * tensor G (eq. 66) and vector g (eq. 69), the LSIC/grad-div term, and the two
 * fine-scale terms that define RBVMS: the cross-stress term (u·(∇w)^T, tau_M
 * r_M) and the Reynolds-stress term -(∇w, tau_M r_M ⊗ tau_M r_M) (eq. 72). This
 * is the matrix-based counterpart of the matrix-free rbvms implementation.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class RBVMSNavierStokesAssemblerCore : public NavierStokesAssemblerBase<dim>
{
public:
  RBVMSNavierStokesAssemblerCore(
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
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
};


/**
 * @brief Class that assembles the coriolis and the centrifugal
 * if a simulation is carried out in a rotating frame of reference
 * according to the following term:
 * \f$2\mathbf{\omega} \times \mathbf{u} + \mathbf{\omega}\times
 * (\mathbf{\omega} \times \mathbf{r})\f$ Where $\mathbf{\omega}$ is the
 * rotation vector of the frame of reference and $\mathbf{r}$ is the position
 * vector (e.g. a vector between a point on the rotation axis and the gauss
 * point) By default, it is assumed that the rotation vector passes through the
 * point (0,0) in 2D or (0,0,0) in 3D
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

  Parameters::VelocitySource velocity_sources;
};


/**
 * @brief Class that assembles the core of the Navier-Stokes equation
 * using a Rheological model to predict non-Newtonian behaviors
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
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief Calculates an approximation of the gradient of the kinematic viscosity
   * @param velocity_gradient The velocity gradient tensor on the quadrature point
     @param velocity_hessians The velocity hessian tensor on the quadrature point
     @param non_newtonian_viscosity The viscosity at which the gradient is calculated
     @param d_gamma_dot The difference in the shear rate magnitude to approximate the
     viscosity variation with a slight change in the shear_rate magnitude
   */
  inline Tensor<1, dim>
  get_kinematic_viscosity_gradient(
    const Tensor<2, dim> &velocity_gradient,
    const Tensor<3, dim> &velocity_hessians,
    const double          shear_rate_magnitude,
    const double          grad_kinematic_viscosity_shear_rate) const
  {
    // Calculates an approximation of the shear rate magnitude gradient using
    // the derived form, since it does not change with rheological models
    Tensor<1, dim> grad_shear_rate;
    for (int d = 0; d < dim; ++d)
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
    return grad_shear_rate * grad_kinematic_viscosity_shear_rate;
  };

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

  /**
   * Enables SUPG stabilization for the Navier-Stokes formulation.
   * We have not found any scenarios where it is relevant not to use SUPG
   * stabilization yet.
   */
  const bool SUPG = true;

  const std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the Navier Stokes equations. For example, if a BDF1 scheme is
 * chosen, the following is assembled
 * \f$\frac{\mathbf{u}^{t+\Delta t}-\mathbf{u}^{t}}{\Delta t}\f$
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
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
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
};

/**
 * @brief Class that assembles the transient term arising from SDIRK time
 * integration for the Navier Stokes equations.
 * For a given stage \f[ i \f], the following expression is assembled and
 * treated implicitly: \f[ \frac{\mathbf{u}_i^* - \mathbf{u}_n}{h a_{ii}} -
 * \sum_{j=0}^{i-1} \frac{a_{ij}}{a_{ii}} \mathbf{k}_j \f] where:
 * - \f[ \mathbf{u}_i^* \f] is the current stage solution,
 * - \f[ \mathbf{u}_n \f] is the solution at the previous time step,
 * - \f[ \mathbf{k}_j \f] coefficients are the solution increments at the
 * previous stages,
 * - \f[ h \f] is the time step size,
 * - \f[ a_{ij} \), \( a_{ii} \f] are coefficients from the SDIRK Butcher
 * tableau. These coefficients are given data by the method and constant for the
 * entire simulation.
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
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
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
};


/**
 * @brief Class that assembles the core of the Navier-Stokes equation.
 * According to the following weak form:
 * \f$\mathbf{u} \cdot \nabla \mathbf{u} - \nabla p - \nu \nabla^2 \mathbf{u}
 * =0 \f$ with a grad-div stabilization
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class BlockNavierStokesAssemblerCore : public NavierStokesAssemblerBase<dim>
{
public:
  BlockNavierStokesAssemblerCore(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const double                              gamma)
    : simulation_control(simulation_control)
    , gamma(gamma)
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
  const double                             gamma;
};


template <int dim>
class BlockNavierStokesAssemblerNonNewtonianCore
  : public NavierStokesAssemblerBase<dim>
{
public:
  BlockNavierStokesAssemblerNonNewtonianCore(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const double                              gamma)
    : simulation_control(simulation_control)
    , gamma(gamma)
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
  const double                             gamma;
};


/**
 * @brief Class that assembles a Poisson problem for all velocity components and pressure variables.
 * According to the following weak form: d^2 U/dx^2=0 and  d^2 P/dx^2=0
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class LaplaceAssembly : public NavierStokesAssemblerBase<dim>
{
public:
  LaplaceAssembly(const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
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
};


/**
 * @brief Class that assembles a Neumann boundary condition.
 * According to the following weak form: (p-mu*grad_u)*n at the boundary
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
    const std::shared_ptr<SimulationControl> &simulation_control,
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
  const BoundaryConditions::NSBoundaryConditions<dim>
    &pressure_boundary_conditions;
};

/**
 * @brief Class that assembles a Neumann traction boundary condition.
 * According to the following weak form: \f$(p-\mu\nabla \mathbf{u})\cdot
 * \mathbf{n}\f$ at the boundary.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @param neumann_traction_boundary_condition The boundary condition objects use to store the function.
 * @ingroup assemblers
 */
template <int dim>
class NeumannTractionBoundaryCondition : public NavierStokesAssemblerBase<dim>
{
public:
  NeumannTractionBoundaryCondition(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const BoundaryConditions::NSBoundaryConditions<dim>
      &neumann_traction_boundary_condition_input)
    : simulation_control(simulation_control)
    , neumann_traction_boundary_condition(
        neumann_traction_boundary_condition_input)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   * @note This method is just for consistency and Neumann boundary traction doesn't
   *       require matrix assembly
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
  const BoundaryConditions::NSBoundaryConditions<dim>
    &neumann_traction_boundary_condition;
};



/**
 * @brief Class that assembles the weak formulation of a Dirichlet boundary condition using the Nitsche method.
 * According to the following weak form: (u_ib-u)-(u,grad v)
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
    const std::shared_ptr<SimulationControl> &simulation_control,
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

  const std::shared_ptr<SimulationControl>             simulation_control;
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions;
};


/**
 * @brief Class that assembles the special case of partial slip condition using the weak formulation of a Dirichlet boundary condition using the Nitsche method.
 * According to the following weak form: beta n(nu-nv) +
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
    const std::shared_ptr<SimulationControl> &simulation_control,
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

  const std::shared_ptr<SimulationControl>             simulation_control;
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions;
};


/**
 * @brief Class that assembles the weak formulation of an outlet boundary condition
 * according to the following equation: (nu grad(u) * grad(v) + pI - (beta *
 * u)_ * n. See the paper by Arndt, Braack and Lube
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
  OutletBoundaryCondition(
    const std::shared_ptr<SimulationControl> &simulation_control,
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

  const std::shared_ptr<SimulationControl>             simulation_control;
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions;
};


/**
 * @brief Class that assembles a thermal buoyancy forcing term using the Boussinesq
 * approximation. For more information, read Chapter 10 of Transport phenomena
 * by Bird et al., or "Boussinesq approximation (buoyancy)" page on Wikipedia.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class ThermalBuoyancyAssembly : public NavierStokesAssemblerBase<dim>
{
public:
  ThermalBuoyancyAssembly(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const double                              reference_temperature)
    : simulation_control(simulation_control)
    , reference_temperature(reference_temperature)
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
   * @brief assemble_rhs Assembles the weak form of: \f$-\mathbf{g} \times \alpha \times (T - T_0)\f$
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;


private:
  const std::shared_ptr<SimulationControl> simulation_control;
  const double                             reference_temperature;
};


/**
 * @brief Class that assembles a phase change Darcy forcing term. This term adds
 * \f$-\beta_D  \mathbf{u} \f$ to the right hand-side of the Navier-Stokes
 * equations to prohibit the motion of a material. In the phase change model,
 * the value of the \f$ \beta_D \f$ coefficient depends on the temperature
 * field. Generally, this is used to prohibit fluid motion in the solid phase
 * within phase change problem. This generally leads to a better conditioning of
 * the linear system than increasing the viscosity of the solid phase.
 *
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class PhaseChangeDarcyAssembly : public NavierStokesAssemblerBase<dim>
{
public:
  PhaseChangeDarcyAssembly(
    const Parameters::PhaseChange phase_change_parameters)
    : phase_change_parameters(phase_change_parameters)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix of: \f$-\beta_D  \mathbf{u} \f$
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;


  /**
   * @brief assemble_rhs Assembles the weak form of: \f$-\beta_D  \mathbf{u} \f$
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;


private:
  /*
   * Phase change parameters are kept within the assembler and are used to
   * calculate, on the fly, the inverse permeability (\f$ \beta_D \f$).
   */
  const Parameters::PhaseChange phase_change_parameters;
};

/**
 * @brief Assembler of the non-linear Carman-Kozeny permeability model for
 * simulations with phase change.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 */
template <int dim>
class PhaseChangeCarmanKozenyAssembler : public NavierStokesAssemblerBase<dim>
{
public:
  /**
   * @brief Assembler for the Carman-Kozeny permeability source term used for
   * phase change modelling.
   *
   * \f[
   * \boldsymbol{F}_\mathrm{Carman-Kozeny} = \frac{-\nu}{A_\mathrm{perm}}
   * \left[ \frac{(1-\alpha_\mathrm{l})^2}{(\alpha_\mathrm{l})^3 +
   * \delta}\right] \boldsymbol{u}
   * \f]
   *
   * with \f$\nu\f$ the kinematic viscosity, \f$A_\mathrm{perm}\f$ the
   * permeability area, \f$\alpha_\mathrm{l}\f$ the liquid fraction,
   * \f$\delta\f$ a tolerance to avoid division by zero in the solid, and
   * \f$\boldsymbol{u}\f$ the velocity.
   *
   * @param[in] carman_kozeny_permeability_area Permeability area of the solid
   * (pseudo-porous media).
   * @param[in] carman_kozeny_tolerance Tolerance to avoid division by zero in
   * the solid.
   * @param[in] liquidus_temperature Liquidus temperature of the fluid used for
   * computing the liquid fraction.
   * @param[in] solidus_temperature Solidus temperature of the fluid used for
   * computing the liquid fraction.
   */
  PhaseChangeCarmanKozenyAssembler(const double carman_kozeny_permeability_area,
                                   const double carman_kozeny_tolerance,
                                   const double liquidus_temperature,
                                   const double solidus_temperature)
    : carman_kozeny_permeability_area_inv(1. / carman_kozeny_permeability_area)
    , carman_kozeny_tolerance(carman_kozeny_tolerance)
    , liquidus_temperature(liquidus_temperature)
    , solidus_temperature(solidus_temperature)
  {}

  /**
   * @brief Assembles the matrix of:
   * \f[
   * \boldsymbol{F}_\mathrm{Carman-Kozeny} = \frac{-\nu}{A_\mathrm{perm}}
   * \left[ \frac{(1-\alpha_\mathrm{l})^2}{(\alpha_\mathrm{l})^3 +
   * \delta}\right] \boldsymbol{u} \f]
   *
   * with \f$\nu\f$ the kinematic viscosity, \f$A_\mathrm{perm}\f$ the
   * permeability area, \f$\alpha_\mathrm{l}\f$ the liquid fraction,
   * \f$\delta\f$ a tolerance to avoid division by zero in the solid, and
   * \f$\boldsymbol{u}\f$ the velocity.
   *
   * @param[in] scratch_data Scratch data containing the information required
   * for system assembly.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   *
   * @param[in,out] copy_data Destination where the local_matrix is copied to.
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Assemble the right-hand side (rhs) of:
   * \f[
   * \boldsymbol{F}_\mathrm{Carman-Kozeny} = \frac{-\nu}{A_\mathrm{perm}}
   * \left[ \frac{(1-\alpha_\mathrm{l})^2}{(\alpha_\mathrm{l})^3 +
   * \delta}\right] \boldsymbol{u} \f]
   *
   * with \f$\nu\f$ the kinematic viscosity, \f$A_\mathrm{perm}\f$ the
   * permeability area, \f$\alpha_\mathrm{l}\f$ the liquid fraction,
   * \f$\delta\f$ a tolerance to avoid division by zero in the solid, and
   * \f$\boldsymbol{u}\f$ the velocity.
   *
   * @param[in] scratch_data Scratch data containing the information required
   * for system assembly.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for rhs assembly.
   *
   * @param[in,out] copy_data Destination where the local_rhs is copied to.
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

private:
  /// One over the permeability area of the pseudo-porous bed (solid phase).
  const double carman_kozeny_permeability_area_inv;

  /// Tolerance in the Carman-Kozeny source term that avoids division by zero.
  const double carman_kozeny_tolerance;

  /// Liquidus temperature
  const double liquidus_temperature;

  /// Solidus temperature
  const double solidus_temperature;
};


#endif
