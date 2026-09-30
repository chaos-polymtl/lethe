// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_cfd_dem.h>

namespace Parameters
{
  namespace
  {
    std::string
    to_string(const VoidFractionMode mode)
    {
      switch (mode)
        {
          case VoidFractionMode::function:
            return "function";
          case VoidFractionMode::pcm:
            return "pcm";
          case VoidFractionMode::qcm:
            return "qcm";
          case VoidFractionMode::spm:
            return "spm";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const QCMFilterType type)
    {
      switch (type)
        {
          case QCMFilterType::spherical:
            return "spherical";
          case QCMFilterType::gaussian:
            return "gaussian";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const VoidFractionQuadratureRule rule)
    {
      switch (rule)
        {
          case VoidFractionQuadratureRule::gauss:
            return "gauss";
          case VoidFractionQuadratureRule::gauss_lobatto:
            return "gauss-lobatto";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const DragModel model)
    {
      switch (model)
        {
          case DragModel::difelice:
            return "difelice";
          case DragModel::rong:
            return "rong";
          case DragModel::dallavalle:
            return "dallavalle";
          case DragModel::kochhill:
            return "kochhill";
          case DragModel::beetstra:
            return "beetstra";
          case DragModel::gidaspow:
            return "gidaspow";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const DragCoupling coupling)
    {
      switch (coupling)
        {
          case DragCoupling::fully_implicit:
            return "implicit";
          case DragCoupling::semi_implicit:
            return "semi-implicit";
          case DragCoupling::fully_explicit:
            return "explicit";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const VANSModel model)
    {
      switch (model)
        {
          case VANSModel::modelA:
            return "modelA";
          case VANSModel::modelB:
            return "modelB";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const SubSimulationControlDEM::DEMSubIterationLogic logic)
    {
      switch (logic)
        {
          case SubSimulationControlDEM::DEMSubIterationLogic::
            fixed_number_of_iterations:
            return "number of iterations";
          case SubSimulationControlDEM::DEMSubIterationLogic::
            fixed_fraction_of_rayleigh_time_step:
            return "fraction of rayleigh time";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const FilterKernelType type)
    {
      switch (type)
        {
          case FilterKernelType::gaussian:
            return "gaussian";
          case FilterKernelType::top_hat:
            return "top-hat";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  template <int dim>
  void
  VoidFractionParameters<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("void fraction");
    prm.declare_entry(
      "mode",
      to_string(mode),
      Patterns::Selection("function|pcm|qcm|spm"),
      "Choose the method for the calculation of the void fraction");
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
      to_string(qcm_filter_type),
      Patterns::Selection("spherical|gaussian"),
      "Filter kernel used by the QCM to weigh particle contributions. With 'spherical' (default), half of 'qcm smoothing length' is the averaging-sphere radius. With 'gaussian', half of 'qcm smoothing length' is the standard deviation sigma of the Gaussian; sigma should be small compared to the QCM neighbor-cell stencil reach to avoid silent truncation bias.");
    prm.declare_entry(
      "quadrature rule",
      to_string(quadrature_rule),
      Patterns::Selection("gauss|gauss-lobatto"),
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
    qcm_smoothing_length       = prm.get_double("qcm smoothing length");
    qcm_sphere_equal_cell_volume = prm.get_bool("qcm sphere equal cell volume");

    const std::string qcm_filter_type_op = prm.get("qcm filter type");
    if (qcm_filter_type_op == "spherical")
      qcm_filter_type = Parameters::QCMFilterType::spherical;
    else if (qcm_filter_type_op == "gaussian")
      qcm_filter_type = Parameters::QCMFilterType::gaussian;
    else
      throw(std::runtime_error(
        "Invalid QCM filter type. Options are 'spherical' or 'gaussian'"));

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
                      to_string(defaults.drag_model),
                      Patterns::Selection(
                        "difelice|rong|dallavalle|kochhill|beetstra|gidaspow"),
                      "The drag model used to determine the drag coefficient");
    prm.declare_entry(
      "dem iteration control",
      to_string(defaults.dem_iteration_control),
      Patterns::Selection("number of iterations|fraction of rayleigh time"),
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
                      to_string(defaults.vans_model),
                      Patterns::Selection("modelA|modelB"),
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
      to_string(defaults.drag_coupling),
      Patterns::Selection("implicit|semi-implicit|explicit"),
      "Formulation for the drag force. Choices are implicit|semi-implicit|explicit. The default value is semi-implicit, which represents the legacy coupling method.");

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

    const std::string it_ctrl = prm.get("dem iteration control");
    if (it_ctrl == "number of iterations")
      dem_iteration_control = SubSimulationControlDEM::DEMSubIterationLogic::
        fixed_number_of_iterations;
    else if (it_ctrl == "fraction of rayleigh time")
      dem_iteration_control = SubSimulationControlDEM::DEMSubIterationLogic::
        fixed_fraction_of_rayleigh_time_step;
    else
      AssertThrow(
        false,
        ExcMessage(
          "An invalid dem iteration control strategy was parsed. Simulation will now stop."));

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

  void
  AndersonJacksonFilter::declare_parameters(ParameterHandler &prm)
  {
    const AndersonJacksonFilter defaults;
    prm.enter_subsection("anderson jackson filter");
    prm.declare_entry("kernel type",
                      to_string(defaults.kernel_type),
                      Patterns::Selection("gaussian|top-hat"),
                      "Kernel of the filter. Choices are <gaussian|top-hat>.");
    prm.declare_entry(
      "filter width",
      Patterns::Tools::Convert<double>::to_string(defaults.filter_width),
      Patterns::Double(0.),
      "Width of the filter. This is the standard deviation of the gaussian "
      "kernel or the radius of the top-hat kernel.");
    prm.declare_entry(
      "gaussian cutoff",
      Patterns::Tools::Convert<double>::to_string(defaults.gaussian_cutoff),
      Patterns::Double(0.),
      "Truncation radius of the gaussian kernel, expressed in number of "
      "standard deviations. The truncated kernel is renormalized to unit "
      "mass.");
    prm.declare_entry("output polynomial degree",
                      Patterns::Tools::Convert<unsigned int>::to_string(
                        defaults.output_degree),
                      Patterns::Integer(0),
                      "Polynomial degree of the FE_Q space on which the "
                      "filtered fields are computed. If 0, the velocity "
                      "degree of the fluid is used.");
    prm.declare_entry("quadrature points",
                      Patterns::Tools::Convert<unsigned int>::to_string(
                        defaults.n_quadrature_points),
                      Patterns::Integer(0),
                      "Number of Gauss quadrature points per direction used "
                      "to integrate the source cells. If 0, the velocity "
                      "degree of the fluid plus one is used.");
    prm.declare_entry("cut cell subdivisions",
                      Patterns::Tools::Convert<unsigned int>::to_string(
                        defaults.cut_cell_subdivisions),
                      Patterns::Integer(1),
                      "Number of subdivisions per direction of the "
                      "quadrature used on the cells cut by an immersed "
                      "solid.");
    prm.declare_entry("filter pressure",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.filter_pressure),
                      Patterns::Bool(),
                      "Filter the pressure in addition to the velocity.");
    prm.declare_entry("minimum fluid fraction",
                      Patterns::Tools::Convert<double>::to_string(
                        defaults.minimum_fluid_fraction),
                      Patterns::Double(0., 1.),
                      "Fluid volume fraction, relative to the kernel mass, "
                      "below which the phase-averaged velocity and pressure "
                      "are not defined and are set to zero.");
    prm.declare_entry("normalize at domain boundaries",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.normalize_at_domain_boundaries),
                      Patterns::Bool(),
                      "Divide the fluid and solid volume fractions by the "
                      "kernel mass. This renormalizes the kernel where it is "
                      "truncated by a domain boundary, which changes the "
                      "definition of the filter near the boundaries.");
    prm.declare_entry("extend velocity beyond walls",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.extend_velocity_beyond_walls),
                      Patterns::Bool(),
                      "Fill the part of the kernel outside of the domain, "
                      "beyond the walls where the velocity is imposed "
                      "(noslip, function and function weak boundary "
                      "conditions), with fluid moving at the velocity of the "
                      "wall when computing the phase-averaged velocity. If "
                      "false, the kernel is truncated by the walls.");
    prm.declare_entry("output folder",
                      defaults.output_folder,
                      Patterns::FileName(),
                      "Folder in which the filtered fields are written.");
    prm.declare_entry("output name",
                      defaults.output_name,
                      Patterns::FileName(),
                      "Prefix of the files in which the filtered fields are "
                      "written.");
    prm.declare_entry("verbosity",
                      to_string(defaults.verbosity),
                      Patterns::Selection("quiet|verbose|extra verbose"),
                      "Verbosity of the filter diagnostics. Choices are "
                      "<quiet|verbose|extra verbose>.");
    prm.leave_subsection();
  }

  void
  AndersonJacksonFilter::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("anderson jackson filter");
    const std::string kernel = prm.get("kernel type");
    if (kernel == "gaussian")
      kernel_type = FilterKernelType::gaussian;
    else if (kernel == "top-hat")
      kernel_type = FilterKernelType::top_hat;
    else
      AssertThrow(false,
                  ExcMessage("Invalid kernel type for the Anderson-Jackson "
                             "filter. Choices are <gaussian|top-hat>."));

    filter_width           = prm.get_double("filter width");
    gaussian_cutoff        = prm.get_double("gaussian cutoff");
    output_degree          = prm.get_integer("output polynomial degree");
    n_quadrature_points    = prm.get_integer("quadrature points");
    cut_cell_subdivisions  = prm.get_integer("cut cell subdivisions");
    filter_pressure        = prm.get_bool("filter pressure");
    minimum_fluid_fraction = prm.get_double("minimum fluid fraction");
    normalize_at_domain_boundaries =
      prm.get_bool("normalize at domain boundaries");
    extend_velocity_beyond_walls = prm.get_bool("extend velocity beyond walls");
    output_folder                = prm.get("output folder");
    output_name                  = prm.get("output name");

    const std::string op = prm.get("verbosity");
    if (op == "quiet")
      verbosity = Verbosity::quiet;
    else if (op == "verbose")
      verbosity = Verbosity::verbose;
    else if (op == "extra verbose")
      verbosity = Verbosity::extra_verbose;
    else
      AssertThrow(false,
                  ExcMessage("Invalid verbosity for the Anderson-Jackson "
                             "filter. Choices are <quiet|verbose|extra "
                             "verbose>."));
    prm.leave_subsection();

    AssertThrow(filter_width > 0.,
                ExcMessage("The filter width of the Anderson-Jackson filter "
                           "must be strictly positive."));
    AssertThrow(kernel_type != FilterKernelType::gaussian ||
                  gaussian_cutoff > 0.,
                ExcMessage("The gaussian cutoff of the Anderson-Jackson filter "
                           "must be strictly positive."));
  }
} // namespace Parameters
// Pre-compile the 2D and 3D
template class Parameters::VoidFractionParameters<2>;
template class Parameters::VoidFractionParameters<3>;
