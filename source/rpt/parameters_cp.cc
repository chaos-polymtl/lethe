#include <rpt/parameters_cp.h>

void
Parameters::CPParameters::declare_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("cp parameters");
    {
        prm.declare_entry("subdivisions",
                    "1",
                    Patterns::Integer(),
                    "Cylinder subdivisions");
        prm.declare_entry("radius",
                          "1",
                          Patterns::Double(),
                          "Cylinder radius");
        prm.declare_entry("half length",
                          "1",
                          Patterns::Double(),
                          "Cylinder half length");
        prm.declare_entry("initial refinement",
                          "1",
                          Patterns::Integer(),
                          "Global refinement");
        prm.declare_entry("Tol electric field",
                          "1",
                          Patterns::Double(),
                          "Tolerance for electric field compartmentalization");
        prm.declare_entry("Tol velocity",
                          "1",
                          Patterns::Double(),
                          "Tolerance for velocity compartmentalization");  
        prm.declare_entry("velocity",
                          "0.001",
                          Patterns::Double(),
                          "Input velocity of CFD simulation");     

    }
    prm.leave_subsection();
}

void
Parameters::CPParameters::parse_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("cp parameters");
    {
        subdivisions              = prm.get_integer("subdivisions");
        cylinder_radius           = prm.get_double("radius");
        cylinder_half_length      = prm.get_double("half length");
        initial_refinement        = prm.get_integer("initial refinement");
        electric_field_tolerance  = prm.get_double("Tol electric field");
        velocity_tolerance        = prm.get_double("Tol velocity");
        CFD_input_velocity        = prm.get_double("velocity");

    }
    prm.leave_subsection();
}

void
CPCalculatingParameters::declare(ParameterHandler &prm)
{
  Parameters::CPParameters::declare_parameters(prm);
}

void
CPCalculatingParameters::parse(ParameterHandler &prm)
{
  cp_param.parse_parameters(prm);

}


