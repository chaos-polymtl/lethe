#include "core/parameters_cfd_dem.h"

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
    prm.declare_entry("shock capturing",
                      "true",
                      Patterns::Bool(),
                      "Choose whether or not to apply shock_capturing");
    prm.declare_entry("grad div",
                      "false",
                      Patterns::Bool(),
                      "Choose whether or not to apply grad_div stabilization");
    prm.declare_entry(
      "full stress tensor",
      "false",
      Patterns::Bool(),
      "Choose whether or not to apply the complete stress tensor");
    prm.declare_entry("drag model",
                      "difelice",
                      Patterns::Selection("difelice|rong"),
                      "The drag model used to determine the drag coefficient");
    prm.declare_entry("reference velocity",
                      "1",
                      Patterns::Double(),
                      "The refernce velocity for the shock capturing scheme");
    prm.leave_subsection();
  }


  void
  CFDDEM::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("cfd-dem");
    shock_capturing      = prm.get_bool("shock capturing");
    grad_div             = prm.get_bool("grad div");
    full_stress_tensor   = prm.get_bool("full stress tensor");
    reference_velocity   = prm.get_double("reference velocity");
    const std::string op = prm.get("drag model");
    if (op == "difelice")
      drag_model = Parameters::DragModel::difelice;
    else if (op == "rong")
      drag_model = Parameters::DragModel::rong;
    else
      throw(std::runtime_error("Invalid drag model"));
    prm.leave_subsection();
  }


} // namespace Parameters
// Pre-compile the 2D and 3D
template class Parameters::VoidFraction<2>;
template class Parameters::VoidFraction<3>;
