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
      Patterns::Selection("function"),
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
    else
      throw(std::runtime_error("Invalid voidfraction model"));
    prm.enter_subsection("function");
    void_fraction.parse_parameters(prm);
    prm.leave_subsection();

    read_dem      = prm.get_bool("read dem");
    dem_file_name = prm.get("dem file name");

    prm.leave_subsection();
  }


} // namespace Parameters
// Pre-compile the 2D and 3D
template class Parameters::VoidFraction<2>;
template class Parameters::VoidFraction<3>;
