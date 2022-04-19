#include "fem-dem/parameters_cfd_dem.h"

namespace Parameters
{
  template <int dim>
  void
  VoidFraction<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("void fraction");
    prm.declare_entry(
      "mode",
      "function",
      Patterns::Selection("function|dem"),
      "Choose the method for the calculation of the void fraction");
    prm.enter_subsection("function");
    void_fraction.declare_parameters(prm, 1);
    prm.leave_subsection();
    prm.declare_entry("read dem",
                      "false",
                      Patterns::Bool(),
                      "Define particles using a DEM simulation results file.");
    prm.declare_entry(
      "bound void fraction",
      "true",
      Patterns::Bool(),
      "Boolean for the bounding of the void fraction using an active set method.");
    prm.declare_entry("dem file name",
                      "dem",
                      Patterns::FileName(),
                      "File output dem prefix");
    prm.declare_entry("l2 smoothing factor",
                      "0.000001",
                      Patterns::Double(),
                      "The smoothing factor for void fraction L2 projection");
    prm.declare_entry("l2 lower bound",
                      "0.36",
                      Patterns::Double(),
                      "The lower bound for void fraction L2 projection");
    prm.declare_entry("l2 upper bound",
                      "1",
                      Patterns::Double(),
                      "The upper bound for void fraction L2 projection");
    prm.leave_subsection();
  }

  template <int dim>
  void
  VoidFraction<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("void fraction");
    const std::string op = prm.get("mode");
    if (op == "function")
      mode = Parameters::VoidFractionMode::function;
    else if (op == "dem")
      mode = Parameters::VoidFractionMode::dem;
    else
      throw(std::runtime_error("Invalid voidfraction model"));
    prm.enter_subsection("function");
    void_fraction.parse_parameters(prm);
    prm.leave_subsection();

    read_dem            = prm.get_bool("read dem");
    bound_void_fraction = prm.get_bool("bound void fraction");
    dem_file_name       = prm.get("dem file name");
    l2_smoothing_factor = prm.get_double("l2 smoothing factor");
    l2_lower_bound      = prm.get_double("l2 lower bound");
    l2_upper_bound      = prm.get_double("l2 upper bound");
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
    prm.declare_entry("drag force",
                      "true",
                      Patterns::Bool(),
                      "Choose whether or not to apply drag force");
    prm.declare_entry("buoyancy force",
                      "true",
                      Patterns::Bool(),
                      "Choose whether or not to apply buoyancy force");
    prm.declare_entry("shear force",
                      "false",
                      Patterns::Bool(),
                      "Choose whether or not to apply shear force");
    prm.declare_entry("pressure force",
                      "false",
                      Patterns::Bool(),
                      "Choose whether or not to apply pressure force");
    prm.declare_entry("drag model",
                      "difelice",
                      Patterns::Selection(
                        "difelice|rong|dallavalle|kochhill|beetstra"),
                      "The drag model used to determine the drag coefficient");
    prm.declare_entry("post processing",
                      "false",
                      Patterns::Bool(),
                      "Choose whether or not to apply post_processing");
    prm.declare_entry("inlet boundary id",
                      "1",
                      Patterns::Integer(),
                      "The inlet boundary of the bed");
    prm.declare_entry("outlet boundary id",
                      "2",
                      Patterns::Integer(),
                      "The outlet boundary of the bed");
    prm.declare_entry("coupling frequency",
                      "100",
                      Patterns::Integer(),
                      "dem-cfd coupling frequency");
    prm.declare_entry("vans model",
                      "modelB",
                      Patterns::Selection("modelA|modelB"),
                      "The volume averaged Navier Stokes model to be solved.");
    prm.declare_entry(
      "grad-div length scale",
      "1",
      Patterns::Double(),
      "Constant cs for the calculation of the grad-div stabilization (gamma = viscosity + cs * velocity)");
    prm.declare_entry(
      "implicit stabilization",
      "true",
      Patterns::Bool(),
      "Choose whether or not to use implicit or explicit stabilization");
    prm.leave_subsection();
  }


  void
  CFDDEM::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("cfd-dem");
    grad_div = prm.get_bool("grad div");
    void_fraction_time_derivative =
      prm.get_bool("void fraction time derivative");
    drag_force             = prm.get_bool("drag force");
    buoyancy_force         = prm.get_bool("buoyancy force");
    shear_force            = prm.get_bool("shear force");
    pressure_force         = prm.get_bool("pressure force");
    post_processing        = prm.get_bool("post processing");
    inlet_boundary_id      = prm.get_integer("inlet boundary id");
    outlet_boundary_id     = prm.get_integer("outlet boundary id");
    coupling_frequency     = prm.get_integer("coupling frequency");
    cstar                  = prm.get_double("grad-div length scale");
    implicit_stabilization = prm.get_bool("implicit stabilization");

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
    else
      throw(std::runtime_error("Invalid drag model"));

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
template class Parameters::VoidFraction<2>;
template class Parameters::VoidFraction<3>;
