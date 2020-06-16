/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Carole-Anne Daunais, Valérie Bibeau, Polytechnique Montreal, 2020-
 */

#include "solvers/gls_nitsche_navier_stokes.h"

#include <deal.II/particles/data_out.h>

#include "core/bdf.h"
#include "core/grids.h"
#include "core/manifolds.h"
#include "core/sdirk.h"
#include "core/time_integration_utilities.h"

// Constructor for class GLSNitscheNavierStokesSolver
template <int dim, int spacedim>
GLSNitscheNavierStokesSolver<dim,spacedim>::GLSNitscheNavierStokesSolver(
  NavierStokesSolverParameters<spacedim> &p_nsparam,
  const unsigned int                      p_degreeVelocity,
  const unsigned int                      p_degreePressure)
  : GLSNavierStokesSolver<spacedim>(
      p_nsparam,
      p_degreeVelocity,
      p_degreePressure)
{}

template <int dim, int spacedim>
GLSNitscheNavierStokesSolver<dim,spacedim>::~GLSNitscheNavierStokesSolver()
{

}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim,spacedim>::assemble_nitsche_restriction()
{
  Particles::ParticleHandler solid_ph = solid.get_solid_particle_handler();
}

template <int dim, int spacedim>
virtual void
GLSNitscheNavierStokesSolver<dim,spacedim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_system");

  if (this->nsparam.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::none)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
        assemble_nitsche_restriction();
    }

  else if (this->nsparam.velocitySource.type ==
           Parameters::VelocitySource::VelocitySourceType::srf)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
        assemble_nitsche_restriction();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
        assemble_nitsche_restriction();
    }
}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class GLSNitscheNavierStokesSolver<2>;
template class GLSNitscheNavierStokesSolver<2,3>;
template class GLSNitscheNavierStokesSolver<3>;
