// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_cfd_dem.h>
#include <core/utilities.h>

namespace Parameters
{
  namespace
  {
    /// Deprecated strings of the void fraction "quadrature rule" parameter
    const DeprecatedEnumNames<VoidFractionQuadratureRule>
      deprecated_quadrature_rule_names = {
        {"gauss-lobatto", VoidFractionQuadratureRule::gauss_lobatto}};

    /// Deprecated strings of the "drag coupling" parameter
    const DeprecatedEnumNames<DragCoupling> deprecated_drag_coupling_names = {
      {"implicit", DragCoupling::fully_implicit},
      {"semi-implicit", DragCoupling::semi_implicit},
      {"explicit", DragCoupling::fully_explicit}};

    /// Deprecated strings of the "dem iteration control" parameter
    const DeprecatedEnumNames<SubSimulationControlDEM::DEMSubIterationLogic>
      deprecated_dem_iteration_control_names = {
        {"number of iterations",
         SubSimulationControlDEM::DEMSubIterationLogic::
           fixed_number_of_iterations},
        {"fraction of rayleigh time",
         SubSimulationControlDEM::DEMSubIterationLogic::
           fixed_fraction_of_rayleigh_time_step}};
  } // namespace

  template <int dim>
  void
  VoidFractionParameters<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("void fraction");
    prm.declare_entry(
      "mode",
      enum_to_string(mode),
      Patterns::Selection(enum_to_selection<Parameters::VoidFractionMode>()),
      "Choose the method for the calculation of the void fraction."
      "Choices are <function|pcm|qcm|spm>.");
    prm.enter_subsection("function");
    void_fraction.declare_parameters(prm);
    prm.leave_subsection();
    prm.declare_entry("read dem",
                      Patterns::Tools::Convert<bool>::to_string(read_dem),
                      Patterns::Bool(),
                      "Define particles using a DEM simulation results file.");
    prm.declare_entry("dem file name",
                      dem_file_name,
                      Patterns::FileName(),
                      "File output dem prefix");
    prm.declare_entry("l2 smoothing length",
                      Patterns::Tools::Convert<double>::to_string(
                        l2_smoothing_length),
                      Patterns::Double(),
                      "The smoothing length for void fraction L2 projection");
    prm.declare_entry(
      "particle refinement factor",
      Patterns::Tools::Convert<unsigned int>::to_string(
        particle_refinement_factor),
      Patterns::Double(),
      "The refinement factor used to calculate the number of pseudo-particles in the satellite point method");
    prm.declare_entry(
      "qcm smoothing length",
      Patterns::Tools::Convert<double>::to_string(qcm_smoothing_length),
      Patterns::Double(),
      "The smoothing length of the QCM filter. With the spherical filter, half of this value is the averaging-sphere radius; with the gaussian filter, half of this value is the standard deviation sigma.");
    prm.declare_entry(
      "qcm sphere equal cell volume",
      Patterns::Tools::Convert<bool>::to_string(qcm_sphere_equal_cell_volume),
      Patterns::Bool(),
      "Specify whether the virtual sphere has the same volume as the mesh element");
    prm.declare_entry(
      "qcm filter type",
      enum_to_string(qcm_filter_type),
      Patterns::Selection(enum_to_selection<Parameters::QCMFilterType>()),
      "Filter kernel used by the QCM to weigh particle contributions. With 'spherical' (default), half of 'qcm smoothing length' is the averaging-sphere radius. With 'gaussian', half of 'qcm smoothing length' is the standard deviation sigma of the Gaussian; sigma should be small compared to the QCM neighbor-cell stencil reach to avoid silent truncation bias.");
    prm.declare_entry(
      "quadrature rule",
      enum_to_string(quadrature_rule),
      Patterns::Selection(enum_to_selection<VoidFractionQuadratureRule>(
        deprecated_quadrature_rule_names)),
      "Choose which quadrature rule to follow when distributing quadrature points for the QCM void fraction scheme");
    prm.declare_entry(
      "n quadrature points",
      Patterns::Tools::Convert<unsigned int>::to_string(n_quadrature_points),
      Patterns::Integer(),
      "Number of quadrature points per cell used in the QCM void fraction scheme");
    prm.declare_entry(
      "project particle velocity",
      Patterns::Tools::Convert<bool>::to_string(project_particle_velocity),
      Patterns::Bool(),
      "Specify whether the particle velocity is projected using QCM");

    prm.leave_subsection();
  }

  template <int dim>
  void
  VoidFractionParameters<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("void fraction");
    mode = string_to_enum<Parameters::VoidFractionMode>(prm.get("mode"));
    prm.enter_subsection("function");
    void_fraction.parse_parameters(prm);
    prm.leave_subsection();

    read_dem                   = prm.get_bool("read dem");
    dem_file_name              = prm.get("dem file name");
    l2_smoothing_length        = prm.get_double("l2 smoothing length");
    particle_refinement_factor = prm.get_integer("particle refinement factor");
    qcm_smoothing_length       = prm.get_double("qcm smoothing length");
    qcm_sphere_equal_cell_volume = prm.get_bool("qcm sphere equal cell volume");

    qcm_filter_type =
      string_to_enum<Parameters::QCMFilterType>(prm.get("qcm filter type"));

    quadrature_rule = string_to_enum<Parameters::VoidFractionQuadratureRule>(
      prm.get("quadrature rule"),
      deprecated_quadrature_rule_names,
      "quadrature rule");

    n_quadrature_points = prm.get_integer("n quadrature points");

    project_particle_velocity = prm.get_bool("project particle velocity");

    prm.leave_subsection();
  }

  void
  CFDDEM::declare_parameters(ParameterHandler &prm)
  {
    const CFDDEM defaults;
    prm.enter_subsection("cfd-dem");
    prm.declare_entry("grad div",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.grad_div),
                      Patterns::Bool(),
                      "Choose whether or not to apply grad_div stabilization");
    prm.declare_entry("void fraction time derivative",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.void_fraction_time_derivative),
                      Patterns::Bool(),
                      "Choose whether or not to implement d(epsilon)/dt ");
    prm.declare_entry(
      "interpolated void fraction",
      Patterns::Tools::Convert<bool>::to_string(
        defaults.interpolated_void_fraction),
      Patterns::Bool(),
      "Choose whether the void fraction is the one of the cell or the one interpolated at the particle position.");
    prm.declare_entry("drag force",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.drag_force),
                      Patterns::Bool(),
                      "Choose whether or not to apply drag force");
    prm.declare_entry("buoyancy force",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.buoyancy_force),
                      Patterns::Bool(),
                      "Choose whether or not to apply buoyancy force");
    prm.declare_entry("shear force",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.shear_force),
                      Patterns::Bool(),
                      "Choose whether or not to apply shear force");
    prm.declare_entry("pressure force",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.pressure_force),
                      Patterns::Bool(),
                      "Choose whether or not to apply pressure force");
    prm.declare_entry("saffman lift force",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.saffman_lift_force),
                      Patterns::Bool(),
                      "Choose whether or not to apply Saffman-Mei lift force");
    prm.declare_entry("magnus lift force",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.magnus_lift_force),
                      Patterns::Bool(),
                      "Choose whether or not to apply Magnus lift force");
    prm.declare_entry(
      "rotational viscous torque",
      Patterns::Tools::Convert<bool>::to_string(
        defaults.rotational_viscous_torque),
      Patterns::Bool(),
      "Choose whether or not to apply rotational viscous torque on particles");
    prm.declare_entry(
      "vortical viscous torque",
      Patterns::Tools::Convert<bool>::to_string(
        defaults.vortical_viscous_torque),
      Patterns::Bool(),
      "Choose whether or not to apply vortical viscous torque on particles");
    prm.declare_entry("drag model",
                      enum_to_string(defaults.drag_model),
                      Patterns::Selection(
                        enum_to_selection<Parameters::DragModel>()),
                      "The drag model used to determine the drag coefficient");
    prm.declare_entry(
      "dem iteration control",
      enum_to_string(defaults.dem_iteration_control),
      Patterns::Selection(
        enum_to_selection<SubSimulationControlDEM::DEMSubIterationLogic>(
          deprecated_dem_iteration_control_names)),
      "The strategy used to control the DEM iterations in CFD-DEM simulations");
    prm.declare_entry("coupling frequency",
                      Patterns::Tools::Convert<unsigned int>::to_string(
                        defaults.coupling_frequency),
                      Patterns::Integer(1),
                      "dem-cfd coupling frequency");
    prm.declare_entry(
      "fraction rayleigh time",
      Patterns::Tools::Convert<double>::to_string(
        defaults.fraction_of_rayleigh_time),
      Patterns::Double(0., 1.),
      "Fraction of Rayleigh time used to control the DEM iterations.");
    prm.declare_entry("vans model",
                      enum_to_string(defaults.vans_model),
                      Patterns::Selection(
                        enum_to_selection<Parameters::VANSModel>()),
                      "The volume averaged Navier Stokes model to be solved.");
    prm.declare_entry(
      "grad-div length scale",
      Patterns::Tools::Convert<double>::to_string(defaults.cstar),
      Patterns::Double(),
      "Constant cs for the calculation of the grad-div stabilization (gamma = kinematic_viscosity + cs * velocity)");
    prm.declare_entry(
      "implicit stabilization",
      Patterns::Tools::Convert<bool>::to_string(
        defaults.implicit_stabilization),
      Patterns::Bool(),
      "Choose whether or not to use implicit or explicit stabilization");

    prm.declare_entry(
      "particle statistics",
      Patterns::Tools::Convert<bool>::to_string(defaults.particle_statistics),
      Patterns::Bool(),
      "Outputs statistics about the particles such as their total kinetic energy, angular momentum, etc.");

    prm.declare_entry(
      "drag coupling",
      enum_to_string(defaults.drag_coupling),
      Patterns::Selection(
        enum_to_selection<DragCoupling>(deprecated_drag_coupling_names)),
      "Formulation for the drag force. Choices are fully_implicit|semi_implicit|fully_explicit. The default value is semi_implicit, which represents the legacy coupling method.");

    prm.declare_entry(
      "project particle forces",
      Patterns::Tools::Convert<bool>::to_string(
        defaults.project_particle_forces),
      Patterns::Bool(),
      "In the VANS solver, specify whether the two-way coupling forces, including the drag, are calculated by projecting the forces acting on the particles onto the fluid grid using the QCM filter.");

    prm.leave_subsection();
  }

  void
  CFDDEM::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("cfd-dem");
    grad_div = prm.get_bool("grad div");
    void_fraction_time_derivative =
      prm.get_bool("void fraction time derivative");
    interpolated_void_fraction = prm.get_bool("interpolated void fraction");
    drag_force                 = prm.get_bool("drag force");
    buoyancy_force             = prm.get_bool("buoyancy force");
    shear_force                = prm.get_bool("shear force");
    pressure_force             = prm.get_bool("pressure force");
    saffman_lift_force         = prm.get_bool("saffman lift force");
    magnus_lift_force          = prm.get_bool("magnus lift force");
    rotational_viscous_torque  = prm.get_bool("rotational viscous torque");
    vortical_viscous_torque    = prm.get_bool("vortical viscous torque");
    coupling_frequency         = prm.get_integer("coupling frequency");
    fraction_of_rayleigh_time  = prm.get_double("fraction rayleigh time");
    cstar                      = prm.get_double("grad-div length scale");
    implicit_stabilization     = prm.get_bool("implicit stabilization");
    particle_statistics        = prm.get_bool("particle statistics");
    project_particle_forces    = prm.get_bool("project particle forces");

    dem_iteration_control =
      string_to_enum<SubSimulationControlDEM::DEMSubIterationLogic>(
        prm.get("dem iteration control"),
        deprecated_dem_iteration_control_names,
        "dem iteration control");

    drag_model = string_to_enum<Parameters::DragModel>(prm.get("drag model"));

    drag_coupling =
      string_to_enum<Parameters::DragCoupling>(prm.get("drag coupling"),
                                               deprecated_drag_coupling_names,
                                               "drag coupling");

    vans_model = string_to_enum<Parameters::VANSModel>(prm.get("vans model"));
    prm.leave_subsection();
  }
} // namespace Parameters
// Pre-compile the 2D and 3D
template class Parameters::VoidFractionParameters<2>;
template class Parameters::VoidFractionParameters<3>;
