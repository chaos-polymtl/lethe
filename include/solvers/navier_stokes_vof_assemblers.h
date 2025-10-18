// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_navier_stokes_vof_assemblers_h
#define lethe_navier_stokes_vof_assemblers_h

#include <core/evaporation_model.h>
#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/navier_stokes_assemblers.h>
#include <solvers/navier_stokes_scratch_data.h>
#include <solvers/simulation_parameters.h>


/**
 * @brief Assembles the core of the Navier-Stokes equation with
 * free surface using VOF modeling.
 * According to the following weak form:
 * \f$ \nabla \cdot \mathbf{u} + \rho \mathbf{u} \cdot \nabla
 * \mathbf{u} + \nabla p - \mu \nabla\cdot (\nabla \mathbf{u} +
 * (\nabla \mathbf{u})^T) = 0 \f$ with an SUPG and PSPG stabilization
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSNavierStokesVOFAssemblerCore : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesVOFAssemblerCore(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const SimulationParameters<dim>          &nsparam)
    : simulation_control(simulation_control)
    , vof_parameters(nsparam.multiphysics.vof_parameters)
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

  const bool SUPG = true;

  const std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::VOF                    vof_parameters;
};


/**
 * @brief Assembles the transient time arising from BDF time
 * integration for the Navier-Stokes equation with
 * free surface using VOF modeling. For example, if a BDF1 scheme is
 * chosen, the following is assembled
 * \f$\frac{(\rho \mathbf{u})^{t+\Delta t}-(\rho \mathbf{u})^{t}}{\Delta t} \f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSNavierStokesVOFAssemblerBDF : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesVOFAssemblerBDF(
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
 * @brief Assembles a phase change Darcy forcing term for a VOF two-fluid
 * simulation. The term  \f$-\phi \beta_D  \mathbf{u} \f$ is added to the right
 * hand-side of the Navier-Stokes equations to prohibit the motion of a
 * material. In the phase change model, the value of the \f$ \beta_D \f$
 * coefficient depends on the temperature field and the material (fluid)
 * properties. Generally, this is used to impose the stasis in the solid
 * phase. This generally leads to a better conditioning of the linear
 * system than increasing the viscosity of the solid phase.
 *
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class PhaseChangeDarcyVOFAssembler : public NavierStokesAssemblerBase<dim>
{
public:
  PhaseChangeDarcyVOFAssembler(
    const std::vector<Parameters::PhaseChange> &phase_change_parameters_vector)
    : phase_change_parameters_vector(phase_change_parameters_vector)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix of: \f$-\beta_D  \mathbf{u}\f$.
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
  const std::vector<Parameters::PhaseChange> phase_change_parameters_vector;
};


/**
 * @brief Assembles the Surface Tension Force (STF) for the
 * Navier-Stokes equations according to the following equation:
 *
 * \f$\mathbf{F_{CSV}}=\sigma k \nabla \phi \frac{2 \rho}{\rho_0 + \rho_1} \f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSNavierStokesVOFAssemblerSTF : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesVOFAssemblerSTF(
    const std::shared_ptr<SimulationControl> &p_simulation_control,
    const SimulationParameters<dim>          &nsparam)
    : simulation_control(p_simulation_control)
    , STF_parameters(nsparam.multiphysics.vof_parameters.surface_tension_force)
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

  // Surface tension force (STF)
  const Parameters::VOF_SurfaceTensionForce STF_parameters;
};


/**
 * @brief Assembles the Marangoni effect for the
 * Navier-Stokes equations according to the following equation:
 *
 * \f$\mathbf{F_{Ma}}= \frac{\partial \sigma}{\partial T} \left[ \nabla T
 * - \frac{\nabla \phi}{| \nabla \phi |} \left( \frac{ \nabla \phi }
 * {| \nabla \phi |} \cdot \nabla T \right) \right] | \nabla \phi |
 * \frac{2 \rho}{\rho_0 + \rho_1} \f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSNavierStokesVOFAssemblerMarangoni
  : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesVOFAssemblerMarangoni(
    const std::shared_ptr<SimulationControl> &p_simulation_control,
    Parameters::VOF_SurfaceTensionForce       p_STF_properties)
    : simulation_control(p_simulation_control)
    , STF_properties(p_STF_properties)
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

  // Surface tension force (STF)
  const Parameters::VOF_SurfaceTensionForce STF_properties;
};


/**
 * @brief Assembles the momentum source due to evaporation for the
 * Navier-Stokes equations.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @param simulation_control Shared pointer of the SimulationControl object
 * controlling the current simulation
 * @param p_evaporation Struct that holds all evaporation model
 * parameters
 *
 * @ingroup assemblers
 */
template <int dim>
class NavierStokesVOFAssemblerEvaporation
  : public NavierStokesAssemblerBase<dim>
{
public:
  NavierStokesVOFAssemblerEvaporation(
    const std::shared_ptr<SimulationControl> &p_simulation_control,
    const Parameters::Evaporation            &p_evaporation)
    : simulation_control(p_simulation_control)
  {
    this->evaporation_model = EvaporationModel::model_cast(p_evaporation);
  }

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

private:
  const std::shared_ptr<SimulationControl> simulation_control;

  // Evaporation model
  std::shared_ptr<EvaporationModel> evaporation_model;
};


/**
 * @brief Assembles the core of the Navier-Stokes equation
 * using a Rheological model to predict non-Newtonian behaviors
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSNavierStokesVOFAssemblerNonNewtonianCore
  : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesVOFAssemblerNonNewtonianCore(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const SimulationParameters<dim>          &nsparam)
    : simulation_control(simulation_control)
    , vof_parameters(nsparam.multiphysics.vof_parameters)
  {}

  /**
   * @brief Calculates an approximation of the gradient of the viscosity
   * @param velocity_gradient The velocity gradient tensor on the quadrature point
     @param velocity_hessians The velocity hessian tensor on the quadrture point
     @param non_newtonian_viscosity The viscosity at which the gradient is calculated
     @param d_gamma_dot Th difference in the shear rate magnitude to approximate the
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
  const Parameters::VOF                    vof_parameters;
};
#endif
