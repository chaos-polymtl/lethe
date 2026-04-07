#ifndef lethe_euler_euler_prm_h
#define lethe_euler_euler_prm_h

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

using namespace dealii;

template <int dim>
struct EulerEulerMeshParameters
{
  Point<dim> p1;
  Point<dim> p2;

  unsigned int nx = 4;
  unsigned int ny = 4;
  unsigned int nz = 4;

  unsigned int global_refinement = 0;

  std::string direction1 = "x";
  std::string direction2 = "y";


  static void
  declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Euler-Euler mesh");
    {
      prm.declare_entry("p1", "0, 0, 0", Patterns::Anything());
      prm.declare_entry("p2", "1, 1, 1", Patterns::Anything());
      prm.declare_entry("nx", "4", Patterns::Integer(1));
      prm.declare_entry("ny", "4", Patterns::Integer(1));
      prm.declare_entry("nz", "4", Patterns::Integer(1));
      prm.declare_entry("global refinement", "0", Patterns::Integer(0));
      prm.declare_entry("direction1", "x", Patterns::Selection("x|y|z"));
      prm.declare_entry("direction2", "y", Patterns::Selection("x|y|z"));
    }
    prm.leave_subsection();
  }

  void
  parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Euler-Euler mesh");
    {
      const auto p1_list = Utilities::split_string_list(prm.get("p1"));
      const auto p2_list = Utilities::split_string_list(prm.get("p2"));

      AssertThrow(p1_list.size() == dim,
                  ExcMessage(
                    "p1 must have the same dimension as the problem."));
      AssertThrow(p2_list.size() == dim,
                  ExcMessage(
                    "p2 must have the same dimension as the problem."));

      const auto p1_values = Utilities::string_to_double(p1_list);
      const auto p2_values = Utilities::string_to_double(p2_list);

      for (unsigned int d = 0; d < dim; ++d)
        {
          p1[d] = p1_values[d];
          p2[d] = p2_values[d];
        }

      nx                = prm.get_integer("nx");
      ny                = prm.get_integer("ny");
      nz                = prm.get_integer("nz");
      global_refinement = prm.get_integer("global refinement");
      direction1        = prm.get("direction1");
      direction2        = prm.get("direction2");
    }
    prm.leave_subsection();
  }
};

struct EulerEulerCouplingParameters
{
  bool verbose = true;

  static void
  declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Euler-Euler coupling");
    {
      prm.declare_entry("verbose", "true", Patterns::Bool());
    }
    prm.leave_subsection();
  }

  void
  parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Euler-Euler coupling");
    {
      verbose = prm.get_bool("verbose");
    }
    prm.leave_subsection();
  }
};

#endif