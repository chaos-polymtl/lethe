// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_time_harmonic_maxwell_assemblers_h
#define lethe_time_harmonic_maxwell_assemblers_h

#include <core/boundary_conditions.h>
#include <core/parameters_multiphysics.h>

#include <solvers/copy_data.h>
#include <solvers/physics_assemblers.h>
#include <solvers/time_harmonic_maxwell_scratch_data.h>

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>

#include <complex>
#include <utility>
#include <vector>

/**
 * @brief A pure virtual class that serves as an interface for all the cell
 * assemblers of the time-harmonic Maxwell equations. The assemblers fill the
 * uncondensed local DPG system stored in the DPGCopyData.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
using TimeHarmonicMaxwellAssemblerBase =
  PhysicsAssemblerBase<TimeHarmonicMaxwellScratchData<dim>, DPGCopyData>;

/**
 * @brief A pure virtual class that serves as an interface for all the face
 * assemblers of the time-harmonic Maxwell equations. In the DPG method, every
 * face of a cell contributes to the local system through the skeleton trial
 * space, so these assemblers are called for each face of each cell once the
 * scratch has been reinitialized on that face.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
using TimeHarmonicMaxwellFaceAssemblerBase =
  PhysicsFaceAssemblerBase<TimeHarmonicMaxwellScratchData<dim>, DPGCopyData>;

/**
 * @brief Check if a boundary condition type of the time-harmonic Maxwell
 * equations is a Robin boundary condition, that is applied weakly through the
 * assembly instead of through the constraints.
 *
 * @param[in] type Boundary condition type.
 *
 * @return True for the silver_muller, impedance_boundary and waveguide_port
 * boundary conditions.
 */
inline bool
is_robin_boundary_type(const BoundaryConditions::BoundaryType type)
{
  return (type == BoundaryConditions::BoundaryType::silver_muller) ||
         (type == BoundaryConditions::BoundaryType::impedance_boundary) ||
         (type == BoundaryConditions::BoundaryType::waveguide_port);
}

/**
 * @brief Compute the incident electromagnetic fields of a rectangular
 * waveguide port from the input parameters. It is only implemented for 3D
 * problems.
 *
 * @tparam dim Spatial dimension.
 *
 * @param[in] time_harmonic_maxwell_parameters Parameters of the time-harmonic
 * Maxwell physics, which contain the waveguide port definitions.
 *
 * @param[in] p Position where the incident fields are computed.
 *
 * @param[in] normal Unit normal vector of the waveguide port face.
 *
 * @param[in] effective_electric_permittivity Effective electric permittivity
 * at the point p.
 *
 * @param[in] effective_magnetic_permeability Effective magnetic permeability
 * at the point p.
 *
 * @param[in] boundary_id_index Index of the waveguide port condition within
 * the waveguide port parameters. The default value is 0, which can be used
 * when there is only one waveguide port defined in the input file.
 *
 * @return A pair containing the incident electric and magnetic fields at the
 * point p.
 */
template <int dim>
std::pair<Tensor<1, dim, std::complex<double>>,
          Tensor<1, dim, std::complex<double>>>
compute_waveguide_port_incident_fields(
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters,
  const Point<dim>                           &p,
  const Tensor<1, dim>                       &normal,
  const std::complex<double>                 &effective_electric_permittivity,
  const std::complex<double>                 &effective_magnetic_permeability,
  const types::boundary_id                    boundary_id_index = 0);

/**
 * @brief Compute the excitation and the surface admittance of the Robin
 * boundary condition of a rectangular waveguide port from the input
 * parameters. It is only implemented for 3D problems.
 *
 * @tparam dim Spatial dimension.
 *
 * @param[in] time_harmonic_maxwell_parameters Parameters of the time-harmonic
 * Maxwell physics, which contain the waveguide port definitions.
 *
 * @param[in] waveguide_ports_electric_amplitudes Electric field amplitudes of
 * all the waveguide ports. The excitation is normalized by their maximum.
 *
 * @param[in] p Position where the excitation is computed.
 *
 * @param[in] normal Unit normal vector of the waveguide port face.
 *
 * @param[in] effective_electric_permittivity Effective electric permittivity
 * at the point p.
 *
 * @param[in] effective_magnetic_permeability Effective magnetic permeability
 * at the point p.
 *
 * @param[in] boundary_id_index Index of the waveguide port condition within
 * the waveguide port parameters. The default value is 0, which can be used
 * when there is only one waveguide port defined in the input file.
 *
 * @return A pair containing the excitation of the Robin boundary condition and
 * the surface admittance at the point p. The excitation is scaled by the ratio
 * of the port amplitude to the maximum amplitude of all the ports.
 */
template <int dim>
std::pair<Tensor<1, dim, std::complex<double>>, std::complex<double>>
compute_waveguide_port_excitation(
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters,
  const std::vector<double>  &waveguide_ports_electric_amplitudes,
  const Point<dim>           &p,
  const Tensor<1, dim>       &normal,
  const std::complex<double> &effective_electric_permittivity,
  const std::complex<double> &effective_magnetic_permeability,
  const types::boundary_id    boundary_id_index = 0);


/**
 * @brief Class that assembles the cell terms of the ultraweak DPG formulation
 * of the time-harmonic Maxwell equations. With the test functions
 * \f$\mathbf{F}\f$ (electric) and \f$\mathbf{I}\f$ (magnetic), it assembles
 * - the graph norm of the test space in the Gram matrix \f$G\f$:
 *   \f$(\mathbf{F}, \mathbf{F}) + (\nabla \times \mathbf{F}, \nabla \times
 *   \mathbf{F}) + (i\omega\varepsilon \mathbf{F}, i\omega\varepsilon
 *   \mathbf{F})\f$ and the equivalent terms coupling \f$\mathbf{F}\f$ and
 *   \f$\mathbf{I}\f$ and for \f$\mathbf{I}\f$ alone;
 * - the interior bilinear form in the matrix \f$B\f$:
 *   \f$(i\omega\varepsilon \mathbf{E}, \mathbf{F}) + (\mathbf{H}, \nabla
 *   \times \mathbf{F}) + (\mathbf{E}, \nabla \times \mathbf{I}) -
 *   (i\omega\mu \mathbf{H}, \mathbf{I})\f$;
 * - the imposed current density in the load vector \f$l\f$:
 *   \f$(\mathbf{J}, \mathbf{F})\f$.
 *
 * The current density is a phasor: the physical current is
 * \f$\mathbf{j}(\mathbf{x},t) = \Re(\mathbf{J}(\mathbf{x}) e^{-i\omega t})\f$,
 * so its imaginary part sets its phase relative to the other sources. The
 * sign of the current density follows from the \f$e^{-i\omega t}\f$ time
 * convention implied by the bilinear form, for which the Ampère law reads
 * \f$\nabla \times \mathbf{H} + i\omega\varepsilon\mathbf{E} = \mathbf{J}\f$.
 * The current density is imposed in the non-dimensional system as given.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class TimeHarmonicMaxwellAssemblerCore
  : public TimeHarmonicMaxwellAssemblerBase<dim>
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] time_harmonic_maxwell_parameters Parameters of the
   * time-harmonic Maxwell physics, used to get the angular frequency.
   */
  TimeHarmonicMaxwellAssemblerCore(const Parameters::TimeHarmonicMaxwell<dim>
                                     &time_harmonic_maxwell_parameters);

  /**
   * @brief Assemble the Gram matrix and the interior bilinear form.
   *
   * @param[in] scratch_data (see base class)
   *
   * @param[in,out] copy_data (see base class)
   */
  void
  assemble_matrix(const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
                  DPGCopyData &copy_data) override;

  /**
   * @brief Assemble the imposed current density in the load vector.
   *
   * @param[in] scratch_data (see base class)
   *
   * @param[in,out] copy_data (see base class)
   */
  void
  assemble_rhs(const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
               DPGCopyData                               &copy_data) override;

private:
  /**
   * Angular frequency of the electromagnetic fields.
   */
  const double omega;
};


/**
 * @brief Class that assembles the skeleton terms of the ultraweak DPG
 * formulation of the time-harmonic Maxwell equations. On every face, it
 * assembles \f$\langle \mathbf{n} \times \hat{\mathbf{E}}, \mathbf{I}
 * \rangle\f$ in the matrix \f$\hat{B}\f$ and, except on the faces where a
 * Robin boundary condition is applied, it assembles
 * \f$\langle \mathbf{n} \times \hat{\mathbf{H}}, \mathbf{F} \rangle\f$. On
 * Robin faces, the latter term is replaced by the Robin boundary condition
 * (see TimeHarmonicMaxwellAssemblerRobinBC).
 *
 * These terms are bilinear in the traces, so the boundary data is not
 * assembled here: the Dirichlet-type conditions (pec, pmc, electric field and
 * magnetic field) are imposed on the traces through the constraints, and the
 * Robin excitations are assembled by TimeHarmonicMaxwellAssemblerRobinBC.
 *
 * The load vector only receives the imposed surface current density
 * \f$\mathbf{K}\f$ on the interior faces. Across a current sheet, with
 * \f$\hat{\mathbf{n}}\f$ pointing from cell 1 to cell 2, the tangential
 * magnetic field jumps as \f$\hat{\mathbf{n}} \times (\mathbf{H}_2 -
 * \mathbf{H}_1) = \mathbf{K}\f$. Since the trace \f$\hat{\mathbf{H}}\f$ is
 * single-valued, it is defined as the average \f$\hat{\mathbf{H}} =
 * (\mathbf{H}_1 + \mathbf{H}_2)/2\f$, which gives
 * \f$\mathbf{n}_i \times \mathbf{H}_i = \mathbf{n}_i \times
 * \hat{\mathbf{H}} - \mathbf{K}/2\f$ for the outward normal of both cells.
 * Each cell adjacent to the face therefore receives
 * \f$\frac{1}{2}\langle \mathbf{K}_t, \mathbf{F} \rangle\f$ in its load
 * vector, where \f$\mathbf{K}_t\f$ is the tangential part of
 * \f$\mathbf{K}\f$. This contribution does not depend on the orientation of
 * the face.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class TimeHarmonicMaxwellAssemblerSkeleton
  : public TimeHarmonicMaxwellFaceAssemblerBase<dim>
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] boundary_conditions Boundary conditions of the time-harmonic
   * Maxwell physics, used to identify the Robin faces.
   */
  TimeHarmonicMaxwellAssemblerSkeleton(
    const BoundaryConditions::TimeHarmonicMaxwellBoundaryConditions<dim>
      &boundary_conditions)
    : boundary_conditions(boundary_conditions)
  {}

  /**
   * @brief Assemble the skeleton terms of the bilinear form.
   *
   * @param[in] scratch_data (see base class)
   *
   * @param[in,out] copy_data (see base class)
   */
  void
  assemble_matrix(const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
                  DPGCopyData &copy_data) override;

  /**
   * @brief Assemble the imposed surface current density in the load vector on
   * the interior faces.
   *
   * @param[in] scratch_data (see base class)
   *
   * @param[in,out] copy_data (see base class)
   */
  void
  assemble_rhs(const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
               DPGCopyData                               &copy_data) override;

private:
  /**
   * Boundary conditions of the time-harmonic Maxwell physics.
   */
  const BoundaryConditions::TimeHarmonicMaxwellBoundaryConditions<dim>
    &boundary_conditions;
};


/**
 * @brief Class that assembles the Robin boundary conditions of the
 * time-harmonic Maxwell equations (silver_muller, impedance_boundary and
 * waveguide_port). On these faces, the magnetic trace is expressed from the
 * electric trace with the surface admittance \f$Y\f$ and the excitation
 * \f$\mathbf{g}\f$, which gives the terms
 * \f$-\langle Y \hat{\mathbf{E}}, \mathbf{F} \rangle\f$ in \f$\hat{B}\f$ and
 * \f$-\langle \mathbf{g}, \mathbf{F} \rangle\f$ in \f$l\f$. Since the energy
 * norm that is minimized is also modified on these faces, it also assembles
 * the boundary terms of the Gram matrix \f$G\f$.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class TimeHarmonicMaxwellAssemblerRobinBC
  : public TimeHarmonicMaxwellFaceAssemblerBase<dim>
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] time_harmonic_maxwell_parameters Parameters of the
   * time-harmonic Maxwell physics, which contain the waveguide port
   * definitions.
   *
   * @param[in] boundary_conditions Boundary conditions of the time-harmonic
   * Maxwell physics.
   *
   * @param[in] waveguide_ports_electric_amplitudes Electric field amplitudes of
   * all the waveguide ports.
   */
  TimeHarmonicMaxwellAssemblerRobinBC(
    const Parameters::TimeHarmonicMaxwell<dim>
      &time_harmonic_maxwell_parameters,
    const BoundaryConditions::TimeHarmonicMaxwellBoundaryConditions<dim>
                              &boundary_conditions,
    const std::vector<double> &waveguide_ports_electric_amplitudes)
    : time_harmonic_maxwell_parameters(time_harmonic_maxwell_parameters)
    , boundary_conditions(boundary_conditions)
    , waveguide_ports_electric_amplitudes(waveguide_ports_electric_amplitudes)
  {}

  /**
   * @brief Assemble the Robin boundary terms of the Gram matrix and of the
   * skeleton bilinear form.
   *
   * @param[in] scratch_data (see base class)
   *
   * @param[in,out] copy_data (see base class)
   */
  void
  assemble_matrix(const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
                  DPGCopyData &copy_data) override;

  /**
   * @brief Assemble the excitation of the Robin boundary conditions in the
   * load vector.
   *
   * @param[in] scratch_data (see base class)
   *
   * @param[in,out] copy_data (see base class)
   */
  void
  assemble_rhs(const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
               DPGCopyData                               &copy_data) override;

private:
  /**
   * @brief Compute the excitation and the surface admittance of the Robin
   * boundary condition at a quadrature point of the current face.
   *
   * @param[in] scratch_data Scratch data reinitialized on the current face.
   *
   * @param[in] q Index of the face quadrature point.
   *
   * @param[in] bc_type Type of the Robin boundary condition of the face.
   *
   * @param[in] waveguide_port_index Index of the waveguide port condition of
   * the face. Only used if bc_type is waveguide_port.
   *
   * @return A pair containing the excitation and the surface admittance.
   */
  std::pair<Tensor<1, dim, std::complex<double>>, std::complex<double>>
  compute_robin_boundary_data(
    const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
    const unsigned int                         q,
    const BoundaryConditions::BoundaryType     bc_type,
    const unsigned int                         waveguide_port_index) const;

  /**
   * @brief Get the index of the waveguide port condition of the current face
   * within the waveguide port parameters.
   *
   * @param[in] scratch_data Scratch data reinitialized on the current face.
   *
   * @return Index of the waveguide port condition.
   */
  unsigned int
  get_waveguide_port_index(
    const TimeHarmonicMaxwellScratchData<dim> &scratch_data) const;

  /**
   * Parameters of the time-harmonic Maxwell physics.
   */
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters;

  /**
   * Boundary conditions of the time-harmonic Maxwell physics.
   */
  const BoundaryConditions::TimeHarmonicMaxwellBoundaryConditions<dim>
    &boundary_conditions;

  /**
   * Electric field amplitudes of all the waveguide ports.
   */
  const std::vector<double> waveguide_ports_electric_amplitudes;
};

#endif
