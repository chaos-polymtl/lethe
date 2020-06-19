#include "core/parameters.h"

namespace Parameters
{
  void
  SimulationControl::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("simulation control");
    {
      prm.declare_entry("method",
                        "steady",
                        Patterns::Selection(
                          "steady|bdf1|bdf2|bdf3|sdirk2|sdirk3"),
                        "The kind of solver for the linear system. "
                        "Choices are <steady|bdf1|bdf2|bdf3|sdirk2|sdirk3>.");
      prm.declare_entry("time step",
                        "1.",
                        Patterns::Double(),
                        "Time step value");
      prm.declare_entry("time end", "1", Patterns::Double(), "Time step value");
      prm.declare_entry("startup time scaling",
                        "0.4",
                        Patterns::Double(),
                        "Scaling factor used in the iterations necessary to "
                        "start-up the BDF schemes.");

      prm.declare_entry("adapt",
                        "false",
                        Patterns::Bool(),
                        "Adaptative time-stepping <true|false>");
      prm.declare_entry("number mesh adapt",
                        "0",
                        Patterns::Integer(),
                        "Number of mesh adaptation (for steady simulations)");
      prm.declare_entry("max cfl",
                        "1",
                        Patterns::Double(),
                        "Maximum CFL value");

      prm.declare_entry("adaptative time step scaling",
                        "1.1",
                        Patterns::Double(),
                        "Adaptative time step scaling");
      prm.declare_entry("output path",
                        "./",
                        Patterns::FileName(),
                        "File output prefix");

      prm.declare_entry("output name",
                        "out",
                        Patterns::FileName(),
                        "File output prefix");

      prm.declare_entry("output frequency",
                        "1",
                        Patterns::Integer(),
                        "Output frequency");

      prm.declare_entry("output time", "1", Patterns::Double(), "Output time");

      prm.declare_entry(
        "output control",
        "iteration",
        Patterns::Selection("iteration|time"),
        "The control for the output of the simulation results"
        "Results can be either outputted at constant iteration frequency or at constant time");

      prm.declare_entry("subdivision",
                        "1",
                        Patterns::Integer(),
                        "Subdivision of mesh cell in postprocessing");

      prm.declare_entry("group files",
                        "1",
                        Patterns::Integer(),
                        "Maximal number of vtu output files");
    }
    prm.leave_subsection();
  }

  void
  SimulationControl::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("simulation control");
    {
      const std::string sv = prm.get("method");
      if (sv == "steady")
        method = TimeSteppingMethod::steady;
      else if (sv == "bdf1")
        method = TimeSteppingMethod::bdf1;
      else if (sv == "bdf2")
        method = TimeSteppingMethod::bdf2;
      else if (sv == "bdf3")
        method = TimeSteppingMethod::bdf3;
      else if (sv == "sdirk2")
        method = TimeSteppingMethod::sdirk2;
      else if (sv == "sdirk3")
        method = TimeSteppingMethod::sdirk3;
      else
        {
          std::runtime_error("Invalid time stepping scheme");
        }

      const std::string osv = prm.get("output control");
      if (osv == "iteration")
        output_control = OutputControl::iteration;
      else if (osv == "time")
        output_control = OutputControl::time;
      else
        {
          std::runtime_error("Invalid output control scheme");
        }
      dt      = prm.get_double("time step");
      timeEnd = prm.get_double("time end");
      adapt   = prm.get_bool("adapt");
      maxCFL  = prm.get_double("max cfl");
      adaptative_time_step_scaling =
        prm.get_double("adaptative time step scaling");
      startup_timestep_scaling = prm.get_double("startup time scaling");
      number_mesh_adaptation   = prm.get_integer("number mesh adapt");
      output_folder            = prm.get("output path");
      output_name              = prm.get("output name");
      output_frequency         = prm.get_integer("output frequency");
      subdivision              = prm.get_integer("subdivision");
      group_files              = prm.get_integer("group files");
    }
    prm.leave_subsection();
  }

  void
  Timer::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("timer");
    {
      prm.declare_entry("type",
                        "none",
                        Patterns::Selection("none|iteration|end"),
                        "Clock monitoring methods "
                        "Choices are <none|iteration|end>.");
    }
    prm.leave_subsection();
  }

  void
  Timer::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("timer");
    {
      const std::string cl = prm.get("type");
      if (cl == "none")
        type = Type::none;
      else if (cl == "iteration")
        type = Type::iteration;
      else if (cl == "end")
        type = Type::end;
    }
    prm.leave_subsection();
  }

  void
  PhysicalProperties::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("physical properties");
    {
      prm.declare_entry("kinematic viscosity",
                        "1",
                        Patterns::Double(),
                        "Kinematic viscosity");
    }
    prm.leave_subsection();
  }

  void
  PhysicalProperties::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("physical properties");
    {
      viscosity = prm.get_double("kinematic viscosity");
    }
    prm.leave_subsection();
  }

  void
  FEM::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("FEM");
    {
      prm.declare_entry("velocity order",
                        "1",
                        Patterns::Integer(),
                        "interpolation order velocity");
      prm.declare_entry("pressure order",
                        "1",
                        Patterns::Integer(),
                        "interpolation order pressure");
      prm.declare_entry("quadrature points",
                        "0",
                        Patterns::Integer(),
                        "interpolation order pressure");
      prm.declare_entry("qmapping all",
                        "false",
                        Patterns::Bool(),
                        "Apply high order mapping everywhere");
    }
    prm.leave_subsection();
  }

  void
  FEM::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("FEM");
    {
      velocity_order           = prm.get_integer("velocity order");
      pressure_order           = prm.get_integer("pressure order");
      number_quadrature_points = prm.get_integer("quadrature points");
      qmapping_all             = prm.get_bool("qmapping all");
    }
    prm.leave_subsection();
  }

  void
  Forces::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("forces");
    {
      prm.declare_entry(
        "verbosity",
        "quiet",
        Patterns::Selection("quiet|verbose"),
        "State whether from the non-linear solver should be printed "
        "Choices are <quiet|verbose>.");
      prm.declare_entry("calculate forces",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of forces");
      prm.declare_entry("calculate torques",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of torques");
      prm.declare_entry("force name",
                        "force",
                        Patterns::FileName(),
                        "File output force prefix");
      prm.declare_entry("torque name",
                        "torque",
                        Patterns::FileName(),
                        "File output force prefix");
      prm.declare_entry("output precision",
                        "10",
                        Patterns::Integer(),
                        "Calculation frequency");
      prm.declare_entry("display precision",
                        "6",
                        Patterns::Integer(),
                        "Output frequency");
      prm.declare_entry("calculation frequency",
                        "1",
                        Patterns::Integer(),
                        "Calculation frequency");
      prm.declare_entry("output frequency",
                        "1",
                        Patterns::Integer(),
                        "Output frequency");
    }
    prm.leave_subsection();
  }

  void
  Forces::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("forces");
    {
      const std::string op = prm.get("verbosity");
      if (op == "verbose")
        verbosity = Verbosity::verbose;
      if (op == "quiet")
        verbosity = Verbosity::quiet;
      calculate_force       = prm.get_bool("calculate forces");
      calculate_torque      = prm.get_bool("calculate torques");
      force_output_name     = prm.get("force name");
      torque_output_name    = prm.get("torque name");
      display_precision     = prm.get_integer("display precision");
      output_precision      = prm.get_integer("output precision");
      calculation_frequency = prm.get_integer("calculation frequency");
      output_frequency      = prm.get_integer("output frequency");
    }
    prm.leave_subsection();
  }

  void
  PostProcessing::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("post-processing");
    {
      prm.declare_entry(
        "verbosity",
        "quiet",
        Patterns::Selection("quiet|verbose"),
        "State whether from the post-processing values should be printed "
        "Choices are <quiet|verbose>.");

      prm.declare_entry(
        "output boundaries",
        "false",
        Patterns::Bool(),
        "Output the boundaries of the domain along with their ID");

      prm.declare_entry(
        "calculate kinetic energy",
        "false",
        Patterns::Bool(),
        "Enable calculation of total kinetic energy. The total kinetic "
        "energy is calculated from the volumetric integral of the kinetic energy over the domain.");

      prm.declare_entry(
        "calculate enstrophy",
        "false",
        Patterns::Bool(),
        "Enable calculation of total enstrophy. The total enstrophy "
        "is calculated from the volumetric integral of the enstrophy over the domain.");

      prm.declare_entry("kinetic energy name",
                        "kinetic_energy",
                        Patterns::FileName(),
                        "File output kinetic energy");

      prm.declare_entry("enstrophy name",
                        "enstrophy",
                        Patterns::FileName(),
                        "File output enstrophy");

      prm.declare_entry("calculation frequency",
                        "1",
                        Patterns::Integer(),
                        "Calculation frequency");

      prm.declare_entry("output frequency",
                        "1",
                        Patterns::Integer(),
                        "Output frequency");
    }
    prm.leave_subsection();
  }

  void
  PostProcessing::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("post-processing");
    {
      const std::string op = prm.get("verbosity");
      if (op == "verbose")
        verbosity = Verbosity::verbose;
      if (op == "quiet")
        verbosity = Verbosity::quiet;

      output_boundaries          = prm.get_bool("output boundaries");
      calculate_kinetic_energy   = prm.get_bool("calculate kinetic energy");
      calculate_enstrophy        = prm.get_bool("calculate enstrophy");
      kinetic_energy_output_name = prm.get("kinetic energy name");
      enstrophy_output_name      = prm.get("enstrophy name");
      calculation_frequency      = prm.get_integer("calculation frequency");
      output_frequency           = prm.get_integer("output frequency");
    }
    prm.leave_subsection();
  }

  void
  NonLinearSolver::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("non-linear solver");
    {
      prm.declare_entry(
        "verbosity",
        "verbose",
        Patterns::Selection("quiet|verbose"),
        "State whether the outputs from the non-linear solver should be printed. "
        "Choices are <quiet|verbose>.");

      prm.declare_entry(
        "solver",
        "newton",
        Patterns::Selection("newton|skip_newton"),
        "Non-linear solver that will be used "
        "Choices are <newton|skip_newton>."
        " The newton solver is a traditional newton solver with"
        "an analytical jacobian formulation. The jacobian matrix and the preconditioner"
        "are assembled every iteration. In the skip_newton method, the jacobian matrix and"
        "the pre-conditioner are re-assembled every skip_iteration.");

      prm.declare_entry("tolerance",
                        "1e-6",
                        Patterns::Double(),
                        "Newton solver tolerance");
      prm.declare_entry("max iterations",
                        "10",
                        Patterns::Integer(),
                        "Maximum number of Newton Iterations");

      prm.declare_entry(
        "skip iterations",
        "1",
        Patterns::Integer(),
        "Non-linear iterations to skip before rebuilding the jacobian matrix "
        "and the preconditioner");

      prm.declare_entry("residual precision",
                        "4",
                        Patterns::Integer(),
                        "Number of digits displayed when showing residuals");
    }
    prm.leave_subsection();
  }

  void
  NonLinearSolver::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("non-linear solver");
    {
      const std::string op = prm.get("verbosity");
      if (op == "verbose")
        verbosity = Parameters::Verbosity::verbose;
      else if (op == "quiet")
        verbosity = Parameters::Verbosity::quiet;
      else
        throw(std::runtime_error("Invalid verbosity level"));

      const std::string str_solver = prm.get("solver");
      if (str_solver == "newton")
        solver = SolverType::newton;
      else if (str_solver == "skip_newton")
        solver = SolverType::skip_newton;
      else
        throw(std::runtime_error("Invalid non-linear solver "));

      tolerance         = prm.get_double("tolerance");
      max_iterations    = prm.get_integer("max iterations");
      skip_iterations   = prm.get_integer("skip iterations");
      display_precision = prm.get_integer("residual precision");
    }
    prm.leave_subsection();
  }

  void
  Mesh::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("mesh");
    {
      prm.declare_entry("type",
                        "dealii",
                        Patterns::Selection("gmsh|dealii"),
                        "Type of mesh "
                        "Choices are <gmsh|dealii>.");

      prm.declare_entry("file name",
                        "none",
                        Patterns::FileName(),
                        "GMSH file name");

      prm.declare_entry("initial refinement",
                        "0",
                        Patterns::Integer(),
                        "Initial refinement of the mesh");

      prm.declare_entry("grid type", "hyper_cube");
      prm.declare_entry("grid arguments", "-1 : 1 : false");
    }
    prm.leave_subsection();
  }

  void
  Mesh::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("mesh");
    {
      {
        const std::string op = prm.get("type");
        if (op == "gmsh")
          type = Type::gmsh;
        else if (op == "dealii")
          type = Type::dealii;
        else
          throw(std::runtime_error(
            "Error, invalid mesh type. Choices are gmsh and dealii"));
      }

      file_name = prm.get("file name");

      initialRefinement = prm.get_integer("initial refinement");

      grid_type      = prm.get("grid type");
      grid_arguments = prm.get("grid arguments");
    }
    prm.leave_subsection();
  }

  void
  LinearSolver::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("linear solver");
    {
      prm.declare_entry(
        "verbosity",
        "verbose",
        Patterns::Selection("quiet|verbose"),
        "State whether output from solver runs should be printed. "
        "Choices are <quiet|verbose>.");
      prm.declare_entry(
        "residual precision",
        "6",
        Patterns::Integer(),
        "Number of digits used when outputing the residual in the terminal");
      prm.declare_entry(
        "method",
        "gmres",
        Patterns::Selection("gmres|bicgstab|amg"),
        "The iterative solver for the linear system of equations. "
        "Choices are <gmres|bicgstab|amg|direct>. gmres is a GMRES iterative solver "
        "with ILU preconditioning. bicgstab is a BICGSTAB iterative solver "
        "with ILU preconditioning. "
        "amg is GMRES + AMG preconditioning with an ILU coarsener and "
        "smoother. On coarse meshes, the gmres/bicgstab solver with ILU "
        "preconditioning is more efficient. "
        "As the number of mesh elements increase, the amg solver is the most "
        "efficient. Generally, at 1M elements, the amg solver always "
        "outperforms the gmres or bicgstab");
      prm.declare_entry("relative residual",
                        "1e-3",
                        Patterns::Double(),
                        "Linear solver residual");
      prm.declare_entry("minimum residual",
                        "1e-8",
                        Patterns::Double(),
                        "Linear solver minimum residual");
      prm.declare_entry("max iters",
                        "1000",
                        Patterns::Integer(),
                        "Maximum solver iterations");

      prm.declare_entry("ilu preconditioner fill",
                        "0",
                        Patterns::Double(),
                        "Ilu preconditioner fill");

      prm.declare_entry("ilu preconditioner absolute tolerance",
                        "1e-8",
                        Patterns::Double(),
                        "Ilu preconditioner tolerance");

      prm.declare_entry("ilu preconditioner relative tolerance",
                        "1.00",
                        Patterns::Double(),
                        "Ilu relative tolerance");

      prm.declare_entry("amg preconditioner ilu fill",
                        "0",
                        Patterns::Double(),
                        "amg preconditioner ilu smoother/coarsener fill");

      prm.declare_entry(
        "amg preconditioner ilu absolute tolerance",
        "1e-12",
        Patterns::Double(),
        "amg preconditioner ilu smoother/coarsener absolute tolerance");

      prm.declare_entry(
        "amg preconditioner ilu relative tolerance",
        "1.00",
        Patterns::Double(),
        "amg preconditioner ilu smoother/coarsener relative tolerance");

      prm.declare_entry("amg aggregation threshold",
                        "1e-14",
                        Patterns::Double(),
                        "amg aggregation threshold");
      prm.declare_entry("amg n cycles",
                        "1",
                        Patterns::Integer(),
                        "amg number of cycles");
      prm.declare_entry("amg w cycles",
                        "false",
                        Patterns::Bool(),
                        "amg w cycling. If this is set to true, W cycling is "
                        "used. Otherwise, V cycling is used.");
      prm.declare_entry("amg smoother sweeps",
                        "2",
                        Patterns::Integer(),
                        "amg smoother sweeps");
      prm.declare_entry("amg smoother overlap",
                        "1",
                        Patterns::Integer(),
                        "amg smoother overlap");
    }
    prm.leave_subsection();
  }
  void
  LinearSolver::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("linear solver");
    {
      const std::string op = prm.get("verbosity");
      if (op == "verbose")
        verbosity = Parameters::Verbosity::verbose;
      else if (op == "quiet")
        verbosity = Parameters::Verbosity::quiet;
      else
        throw(
          std::runtime_error("Unknown verbosity mode for the linear solver"));

      const std::string sv = prm.get("method");
      if (sv == "amg")
        solver = SolverType::amg;
      else if (sv == "gmres")
        solver = SolverType::gmres;
      else if (sv == "bicgstab")
        solver = SolverType::bicgstab;
      else
        throw std::runtime_error(
          "Error, invalid iterative solver type. Choices are amg, gmres or bicgstab");

      residual_precision = prm.get_integer("residual precision");
      relative_residual  = prm.get_double("relative residual");
      minimum_residual   = prm.get_double("minimum residual");
      max_iterations     = prm.get_integer("max iters");
      ilu_precond_fill   = prm.get_double("ilu preconditioner fill");
      ilu_precond_atol =
        prm.get_double("ilu preconditioner absolute tolerance");
      ilu_precond_rtol =
        prm.get_double("ilu preconditioner relative tolerance");
      amg_precond_ilu_fill = prm.get_double("amg preconditioner ilu fill");
      amg_precond_ilu_atol =
        prm.get_double("amg preconditioner ilu absolute tolerance");
      amg_precond_ilu_rtol =
        prm.get_double("amg preconditioner ilu relative tolerance");
      amg_aggregation_threshold = prm.get_double("amg aggregation threshold");
      amg_n_cycles              = prm.get_integer("amg n cycles");
      amg_w_cycles              = prm.get_bool("amg w cycles");
      amg_smoother_sweeps       = prm.get_integer("amg smoother sweeps");
      amg_smoother_overlap      = prm.get_integer("amg smoother overlap");
    }
    prm.leave_subsection();
  }
  void
  MeshAdaptation::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("mesh adaptation");
    {
      prm.declare_entry("type",
                        "none",
                        Patterns::Selection("none|uniform|kelly"),
                        "Type of mesh adaptation"
                        "Choices are <none|uniform|kelly>.");

      prm.declare_entry("variable",
                        "velocity",
                        Patterns::Selection("velocity|pressure"),
                        "Variable for kelly estimation"
                        "Choices are <velocity|pressure>.");
      prm.declare_entry(
        "fraction type",
        "number",
        Patterns::Selection("number|fraction"),
        "How the fraction of refinement/coarsening are interepreted"
        "Choices are <number|fraction>.");
      prm.declare_entry("max number elements",
                        "100000000",
                        Patterns::Integer(),
                        "Maximum number of elements");
      prm.declare_entry("max refinement level",
                        "10",
                        Patterns::Integer(),
                        "Maximum refinement level");
      prm.declare_entry("min refinement level",
                        "0",
                        Patterns::Integer(),
                        "Minimum refinement level");
      prm.declare_entry("frequency",
                        "1",
                        Patterns::Integer(),
                        "Frequency of the mesh refinement");
      prm.declare_entry("fraction refinement",
                        "0.1",
                        Patterns::Double(),
                        "Fraction of refined elements");
      prm.declare_entry("fraction coarsening",
                        "0.05",
                        Patterns::Double(),
                        "Fraction of coarsened elements");
    }
    prm.leave_subsection();
  }

  void
  MeshAdaptation::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("mesh adaptation");
    {
      const std::string op = prm.get("type");
      if (op == "none")
        type = Type::none;
      if (op == "uniform")
        type = Type::uniform;
      if (op == "kelly")
        type = Type::kelly;

      const std::string vop = prm.get("variable");
      if (vop == "velocity")
        variable = Variable::velocity;
      if (vop == "pressure")
        variable = Variable::pressure;

      const std::string fop = prm.get("fraction type");
      if (fop == "number")
        fractionType = FractionType::number;
      if (fop == "fraction")
        fractionType = FractionType::fraction;
      maximum_number_elements  = prm.get_integer("max number elements");
      maximum_refinement_level = prm.get_integer("max refinement level");
      minimum_refinement_level = prm.get_integer("min refinement level");
      frequency                = prm.get_integer("frequency");
      coarsening_fraction      = prm.get_double("fraction coarsening");
      refinement_fraction      = prm.get_double("fraction refinement");
    }
    prm.leave_subsection();
  }

  void
  Testing::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("test");
    {
      prm.declare_entry("enable",
                        "false",
                        Patterns::Bool(),
                        "Enable testing mode of a solver");
    }
    prm.leave_subsection();
  }

  void
  Testing::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("test");
    {
      enabled = prm.get_bool("enable");
    }
    prm.leave_subsection();
  }

  void
  Restart::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("restart");
    {
      prm.declare_entry("filename",
                        "restart",
                        Patterns::FileName(),
                        "Prefix for the filename of checkpoints");
      prm.declare_entry("restart",
                        "false",
                        Patterns::Bool(),
                        "Frequency for checkpointing");
      prm.declare_entry("checkpoint",
                        "false",
                        Patterns::Bool(),
                        "Enable checkpointing");
      prm.declare_entry("frequency",
                        "1",
                        Patterns::Integer(),
                        "Frequency for checkpointing");
    }
    prm.leave_subsection();
  }

  void
  Restart::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("restart");
    {
      filename   = prm.get("filename");
      checkpoint = prm.get_bool("checkpoint");
      restart    = prm.get_bool("restart");
      frequency  = prm.get_integer("frequency");
    }
    prm.leave_subsection();
  }

  void
  VelocitySource::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("velocity source");
    {
      prm.declare_entry("type",
                        "none",
                        Patterns::Selection("none|srf"),
                        "Velocity-dependent source terms"
                        "Choices are <none|srf>.");
      prm.declare_entry(
        "omega_x",
        "0 ",
        Patterns::Double(),
        "X component of the angular velocity vector of the frame of reference");

      prm.declare_entry(
        "omega_y",
        "0 ",
        Patterns::Double(),
        "Y component of the angular velocity vector of the frame of reference");

      prm.declare_entry(
        "omega_z",
        "0 ",
        Patterns::Double(),
        "Z component of the angular velocity vector of the frame of reference");
    }
    prm.leave_subsection();
  }

  void
  VelocitySource::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("velocity source");
    {
      const std::string op = prm.get("type");
      if (op == "none")
        type = VelocitySourceType::none;
      else if (op == "srf")
        type = VelocitySourceType::srf;
      else
        throw("Error, invalid velocity source type");

      omega_x = prm.get_double("omega_x");
      omega_y = prm.get_double("omega_y");
      omega_z = prm.get_double("omega_z");
    }
    prm.leave_subsection();
  }
} // namespace Parameters
