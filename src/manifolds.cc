#include "manifolds.h"

namespace Parameters
{
  void
  Manifolds::declareDefaultEntry(ParameterHandler &prm, unsigned int i_bc)
  {
    prm.declare_entry("type",
                      "none",
                      Patterns::Selection("none|spherical"),
                      "Type of manifold description"
                      "Choices are <none|spherical>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(i_bc, 2),
                      Patterns::Integer(),
                      "Mesh id for boundary conditions");

    prm.declare_entry("arg1",
                      "0",
                      Patterns::Double(),
                      "Argument of construction no. 1");
    prm.declare_entry("arg2",
                      "0",
                      Patterns::Double(),
                      "Argument of construction no. 2");
    prm.declare_entry("arg3",
                      "0",
                      Patterns::Double(),
                      "Argument of construction no. 3");
    prm.declare_entry("arg4",
                      "0",
                      Patterns::Double(),
                      "Argument of construction no. 4");
    prm.declare_entry("arg5",
                      "0",
                      Patterns::Double(),
                      "Argument of construction no. 5");
  }

  void
  Manifolds::parse_boundary(ParameterHandler &prm, unsigned int i_bc)
  {
    const std::string op = prm.get("type");
    if (op == "none")
      types[i_bc] = none;
    else if (op == "spherical")
      types[i_bc] = spherical;

    id[i_bc] = prm.get_integer("id");
  }

  void
  Manifolds::declare_parameters(ParameterHandler &prm)
  {
    max_size = 5;
    arg1.resize(max_size);
    arg2.resize(max_size);
    arg3.resize(max_size);
    arg4.resize(max_size);
    arg5.resize(max_size);

    prm.enter_subsection("manifolds");
    {
      prm.declare_entry("number",
                        "0",
                        Patterns::Integer(),
                        "Number of boundary conditions");
      id.resize(max_size);
      types.resize(max_size);

      prm.enter_subsection("manifold 0");
      declareDefaultEntry(prm, 0);
      prm.leave_subsection();

      prm.enter_subsection("manifold 1");
      declareDefaultEntry(prm, 1);
      prm.leave_subsection();

      prm.enter_subsection("manifold 2");
      declareDefaultEntry(prm, 2);
      prm.leave_subsection();

      prm.enter_subsection("manifold 3");
      declareDefaultEntry(prm, 3);
      prm.leave_subsection();

      prm.enter_subsection("manifold 4");
      declareDefaultEntry(prm, 4);
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  void
  Manifolds::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("manifolds");
    {
      size = prm.get_integer("number");
      types.resize(size);
      id.resize(size);
      if (size >= 1)
        {
          prm.enter_subsection("manifold 0");
          parse_boundary(prm, 0);
          prm.leave_subsection();
        }
      if (size >= 2)
        {
          prm.enter_subsection("manifold 1");
          parse_boundary(prm, 1);
          prm.leave_subsection();
        }
      if (size >= 3)
        {
          prm.enter_subsection("manifold 2");
          parse_boundary(prm, 2);
          prm.leave_subsection();
        }
      if (size >= 4)
        {
          prm.enter_subsection("manifold 3");
          parse_boundary(prm, 3);
          prm.leave_subsection();
        }
      if (size >= 5)
        {
          prm.enter_subsection("manifold 4");
          parse_boundary(prm, 4);
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }
} // namespace Parameters
