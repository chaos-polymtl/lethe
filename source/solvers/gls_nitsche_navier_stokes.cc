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
 * Author: Bruno Blais, Carole-Anne Daunais, Valérie Bibeau, Polytechnique
 Montreal, 2020-
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
GLSNitscheNavierStokesSolver<dim, spacedim>::GLSNitscheNavierStokesSolver(
  NavierStokesSolverParameters<spacedim> &p_nsparam,
  const unsigned int                      p_degreeVelocity,
  const unsigned int                      p_degreePressure)
  : GLSNavierStokesSolver<spacedim>(p_nsparam,
                                    p_degreeVelocity,
                                    p_degreePressure)
  , solid(this->nsparam.nitsche, this->triangulation, p_degreeVelocity)
{}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::assemble_nitsche_restriction()
{
  std::shared_ptr<Particles::ParticleHandler<spacedim>> solid_ph =
    solid.get_solid_particle_handler();

  TimerOutput::Scope t(this->computing_timer, "Assemble Nitsche terms");

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(spacedim);

  const unsigned int dofs_per_cell = this->fe.dofs_per_cell;

  std::vector<types::global_dof_index> fluid_dof_indices(dofs_per_cell);

  FullMatrix<double>     local_matrix(dofs_per_cell, dofs_per_cell);
  dealii::Vector<double> local_rhs(dofs_per_cell);

  // Penalization terms
  const auto penalty_parameter =
    1.0 / GridTools::minimal_cell_diameter(*this->triangulation);
  double beta = this->nsparam.nitsche.beta;

  // Loop over all local particles
  auto particle = solid_ph->begin();
  while (particle != solid_ph->end())
    {
      local_matrix = 0;
      local_rhs    = 0;

      const auto &cell = particle->get_surrounding_cell(*this->triangulation);
      const auto &dh_cell =
        typename DoFHandler<spacedim>::cell_iterator(*cell, &this->dof_handler);
      dh_cell->get_dof_indices(fluid_dof_indices);

      const auto pic = solid_ph->particles_in_cell(cell);
      Assert(pic.begin() == particle, ExcInternalError());
      for (const auto &p : pic)
        {
          const auto &ref_q = p.get_reference_location();
          // const auto &real_q  = p.get_location();
          const auto &JxW = p.get_properties()[0];

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const auto comp_i = this->fe.system_to_component_index(i).first;
              if (comp_i < spacedim)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      const auto comp_j =
                        this->fe.system_to_component_index(j).first;
                      if (comp_i == comp_j)
                        local_matrix(i, j) += penalty_parameter * beta *
                                              this->fe.shape_value(i, ref_q) *
                                              this->fe.shape_value(j, ref_q) *
                                              JxW;
                    }
                  local_rhs(i) += 0;
                }
            }
        }
      const AffineConstraints<double> &constraints_used =
        this->zero_constraints;
      constraints_used.distribute_local_to_global(local_matrix,
                                                  local_rhs,
                                                  fluid_dof_indices,
                                                  this->system_matrix,
                                                  this->system_rhs);
      particle = pic.end();
    }
  this->system_matrix.compress(VectorOperation::add);
  this->system_rhs.compress(VectorOperation::add);
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_system");

  if (this->nsparam.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::none)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::bdf1,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::bdf2,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::bdf3,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::steady,
          Parameters::VelocitySource::VelocitySourceType::none>();
    }

  else if (this->nsparam.velocitySource.type ==
           Parameters::VelocitySource::VelocitySourceType::srf)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::bdf1,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::bdf2,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::bdf3,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        GLSNavierStokesSolver<spacedim>::template assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::steady,
          Parameters::VelocitySource::VelocitySourceType::srf>();
    }
  assemble_nitsche_restriction();
}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class GLSNitscheNavierStokesSolver<2>;
template class GLSNitscheNavierStokesSolver<2, 3>;
template class GLSNitscheNavierStokesSolver<3>;
