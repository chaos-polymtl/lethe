// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file qmom_source_terms.h
 * @brief Point-wise phase-space step of the QMOM physics.
 *
 * The QMOM physics is advanced with operator splitting: the physical-space
 * transport of all the fields is solved with finite elements by Lethe, and
 * the phase-space evolution (nucleation, growth and coalescence) is applied
 * independently at every support point of the finite element space. This file
 * defines the contract of that point-wise step, which is delegated to
 * QMOMKit.
 */

#ifndef lethe_qmom_source_terms_h
#define lethe_qmom_source_terms_h

#include <core/parameters.h>
#include <core/parameters_multiphysics.h>

#include <span>

/**
 * @brief Advance the local state of the QMOM physics over a time step under
 * the phase-space sources (nucleation, growth and coalescence), without
 * transport. The state is updated in place.
 *
 * The state holds, in this order, the ethanol concentration, the free lipids
 * concentration and the raw moments M_0 to M_{2N-1} of the particle size
 * distribution, where N is the number of quadrature nodes. This layout is the
 * component layout of the QMOM FESystem.
 *
 * The sources are not implemented yet: the state is left unchanged.
 *
 * @param[in,out] state Ethanol concentration, free lipids concentration and
 * raw moments at one support point.
 * @param[in] time_step Duration of the phase-space step.
 * @param[in] qmom_parameters Parameters of the QMOM physics (enabled sources,
 * operating conditions and constants).
 * @param[in] phase_properties Properties of the aqueous and organic phases.
 */
void
apply_qmom_sources(std::span<double>                      state,
                   const double                           time_step,
                   const Parameters::QMOM                &qmom_parameters,
                   const Parameters::QMOMPhaseProperties &phase_properties);

#endif
