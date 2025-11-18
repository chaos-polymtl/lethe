// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_cfd_dem.h>

namespace Parameters
{
  template <int dim>
  void
  VoidFractionParameters<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("void fraction");
    prm.declare_entry(
      "mode",
      "function",
      Patterns::Selection("function|pcm|qcm|spm"),
      "Choose the method for the calculation of the void fraction");
    prm.enter_subsection("function");
    void_fraction.declare_parameters(prm);
    prm.leave_subsection();
    prm.declare_entry("read dem",
                      "false",
                      Patterns::Bool(),
                      "Define particles using a DEM simulation results file.");
    prm.declare_entry("dem file name",
                      "dem",
                      Patterns::FileName(),
                      "File output dem prefix");
    prm.declare_entry("l2 smoothing length",
                      "0.001",
                      Patterns::Double(),
                      "The smoothing length for void fraction L2 projection");
    prm.declare_entry(
      "particle refinement factor",
      "0",
      Patterns::Double(),
      "The refinement factor used to calculate the number of pseudo-particles in the satellite point method");
    prm.declare_entry(
      "qcm sphere diameter",
      "0",
      Patterns::Double(),
      "The diameter of the reference sphere for QCM void fraction scheme");
    prm.declare_entry(
      "qcm sphere equal cell volume",
      "false",
      Patterns::Bool(),
      "Specify whether the virtual sphere has the same volume as the mesh element");
    prm.declare_entry(
      "quadrature rule",
      "gauss",
      Patterns::Selection("gauss|gauss-lobatto"),
      "Choose which quadrature rule to follow when distributing quadrature points for the QCM void fraction scheme");
    prm.declare_entry(
      "n quadrature points",
      "0",
      Patterns::Integer(),
      "Number of quadrature points per cell used in the QCM void fraction scheme");
    prm.declare_entry(
      "project particle velocity",
      "false",
      Patterns::Bool(),
      "Specify whether the particle velocity is projected using QCM");
    prm.declare_entry(
      "project particle forces",
      "false",
      Patterns::Bool(),
      "Specify whether in the VANS equations, the two-way coupling forces, including the drag, are calculated by projecting the forces acting on the particles to the fluid grid.");

    prm.leave_subsection();
  }

  template <int dim>
  void
  VoidFractionParameters<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("void fraction");
    const std::string op = prm.get("mode");
    if (op == "function")
      mode = Parameters::VoidFractionMode::function;
    else if (op == "pcm")
      mode = Parameters::VoidFractionMode::pcm;
    else if (op == "qcm")
      mode = Parameters::VoidFractionMode::qcm;
    else if (op == "spm")
      mode = Parameters::VoidFractionMode::spm;
    else
      throw(std::runtime_error("Invalid void fraction calculation scheme"));
    prm.enter_subsection("function");
    void_fraction.parse_parameters(prm);
    prm.leave_subsection();

    read_dem                   = prm.get_bool("read dem");
    dem_file_name              = prm.get("dem file name");
    l2_smoothing_length        = prm.get_double("l2 smoothing length");
    particle_refinement_factor = prm.get_integer("particle refinement factor");
    qcm_sphere_diameter        = prm.get_double("qcm sphere diameter");
    qcm_sphere_equal_cell_volume = prm.get_bool("qcm sphere equal cell volume");
    const std::string quadrature_rule_op = prm.get("quadrature rule");

    if (quadrature_rule_op == "gauss")
      quadrature_rule = Parameters::VoidFractionQuadratureRule::gauss;
    else if (quadrature_rule_op == "gauss-lobatto")
      quadrature_rule = Parameters::VoidFractionQuadratureRule::gauss_lobatto;
    else
      throw(std::runtime_error(
        "Invalid quadrature rule for the void fraction calculation scheme. Options are 'gauss' or 'gauss-lobatto'"));

    n_quadrature_points = prm.get_integer("n quadrature points");

    project_particle_velocity = prm.get_bool("project particle velocity");

    project_particle_forces = prm.get_bool("project particle forces");

    // if (project_particle_forces && !project_particle_velocity)
    //   throw(std::runtime_error(
    //     "'project particle forces' requires 'project particle velocity' to be true as well."));

    prm.leave_subsection();
  }

  void
  CFDDEM::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("cfd-dem");
    prm.declare_entry("grad div",
                      "true",
                      Patterns::Bool(),
                      "Choose whether or not to apply grad_div stabilization");
    prm.declare_entry("void fraction time derivative",
                      "true",
                      Patterns::Bool(),
                      "Choose whether or not to implement d(epsilon)/dt ");
    prm.declare_entry(
      "interpolated void fraction",
      "true",
      Patterns::Bool(),
      "Choose whether the void fraction is the one of the cell or the one interpolated at the particle position.");
    prm.declare_entry("drag force",
                      "true",
                      Patterns::Bool(),
                      "Choose whether or not to apply drag force");
    prm.declare_entry("buoyancy force",
                      "true",
                      Patterns::Bool(),
                      "Choose whether or not to apply buoyancy force");
    prm.declare_entry("shear force",
                      "true",
                      Patterns::Bool(),
                      "Choose whether or not to apply shear force");
    prm.declare_entry("pressure force",
                      "true",
                      Patterns::Bool(),
                      "Choose whether or not to apply pressure force");
    prm.declare_entry("saffman lift force",
                      "false",
                      Patterns::Bool(),
                      "Choose whether or not to apply Saffman-Mei lift force");
    prm.declare_entry("magnus lift force",
                      "false",
                      Patterns::Bool(),
                      "Choose whether or not to apply Magnus lift force");
    prm.declare_entry(
      "rotational viscous torque",
      "false",
      Patterns::Bool(),
      "Choose whether or not to apply rotational viscous torque on particles");
    prm.declare_entry(
      "vortical viscous torque",
      "false",
      Patterns::Bool(),
      "Choose whether or not to apply vortical viscous torque on particles");
    prm.declare_entry("drag model",
                      "difelice",
                      Patterns::Selection(
                        "difelice|rong|dallavalle|kochhill|beetstra|gidaspow"),
                      "The drag model used to determine the drag coefficient");
    prm.declare_entry("coupling frequency",
                      "100",
                      Patterns::Integer(),
                      "dem-cfd coupling frequency");
    prm.declare_entry("vans model",
                      "modelA",
                      Patterns::Selection("modelA|modelB"),
                      "The volume averaged Navier Stokes model to be solved.");
    prm.declare_entry(
      "grad-div length scale",
      "1",
      Patterns::Double(),
      "Constant cs for the calculation of the grad-div stabilization (gamma = kinematic_viscosity + cs * velocity)");
    prm.declare_entry(
      "implicit stabilization",
      "true",
      Patterns::Bool(),
      "Choose whether or not to use implicit or explicit stabilization");

    prm.declare_entry(
      "particle statistics",
      "true",
      Patterns::Bool(),
      "Outputs statistics about the particles such as their total kinetic energy, angular momentum, etc.");

    prm.declare_entry(
      "drag coupling",
      "semi-implicit",
      Patterns::Selection("implicit|semi-implicit|explicit"),
      "Formulation for the drag force. Choices are implicit|semi-implicit|explicit. The default value is semi-implicit, which represents the legacy coupling method.");

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
    cstar                      = prm.get_double("grad-div length scale");
    implicit_stabilization     = prm.get_bool("implicit stabilization");
    particle_statistics        = prm.get_bool("particle statistics");

    const std::string op = prm.get("drag model");
    if (op == "difelice")
      drag_model = Parameters::DragModel::difelice;
    else if (op == "rong")
      drag_model = Parameters::DragModel::rong;
    else if (op == "dallavalle")
      drag_model = Parameters::DragModel::dallavalle;
    else if (op == "kochhill")
      drag_model = Parameters::DragModel::kochhill;
    else if (op == "beetstra")
      drag_model = Parameters::DragModel::beetstra;
    else if (op == "gidaspow")
      drag_model = Parameters::DragModel::gidaspow;
    else
      AssertThrow(false, ExcMessage("Invalid drag model"));

    const std::string drag_coupling_str = prm.get("drag coupling");
    if (drag_coupling_str == "implicit")
      drag_coupling = Parameters::DragCoupling::fully_implicit;
    else if (drag_coupling_str == "explicit")
      drag_coupling = Parameters::DragCoupling::fully_explicit;
    else if (drag_coupling_str == "semi-implicit")
      drag_coupling = Parameters::DragCoupling::semi_implicit;
    else
      AssertThrow(false, ExcMessage("Drag coupling formulation"));

    const std::string op1 = prm.get("vans model");
    if (op1 == "modelA")
      vans_model = Parameters::VANSModel::modelA;
    else if (op1 == "modelB")
      vans_model = Parameters::VANSModel::modelB;
    else
      throw(std::runtime_error(
        "Invalid vans model. Valid choices are modelA and modelB."));
    prm.leave_subsection();
  }
} // namespace Parameters
// Pre-compile the 2D and 3D
template class Parameters::VoidFractionParameters<2>;
template class Parameters::VoidFractionParameters<3>;
