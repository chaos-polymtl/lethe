#ifndef lethe_solid_phase_h
#define lethe_solid_phase_h

#include <core/boundary_conditions.h>
#include <core/pvd_handler.h>
#include <core/simulation_control.h>

#include <solvers/fluid_dynamics_matrix_based.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>



template <int dim>
struct SolidVelocityFunctionData
{
  Functions::ParsedFunction<dim> u;
  Functions::ParsedFunction<dim> v;
  Functions::ParsedFunction<dim> w;

  SolidVelocityFunctionData()
    : u(1)
    , v(1)
    , w(1)
  {}

  static void
  declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("u");
    {
      Functions::ParsedFunction<dim> tmp(1);
      tmp.declare_parameters(prm);
    }
    prm.leave_subsection();

    prm.enter_subsection("v");
    {
      Functions::ParsedFunction<dim> tmp(1);
      tmp.declare_parameters(prm);
    }
    prm.leave_subsection();

    prm.enter_subsection("w");
    {
      Functions::ParsedFunction<dim> tmp(1);
      tmp.declare_parameters(prm);
    }
    prm.leave_subsection();
  }

  void
  parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("u");
    u.parse_parameters(prm);
    prm.leave_subsection();
    prm.enter_subsection("v");
    v.parse_parameters(prm);
    prm.leave_subsection();
    prm.enter_subsection("w");
    w.parse_parameters(prm);
    prm.leave_subsection();
  }
};


template <int dim>
struct SolidInitialConditionsParameters
{
  std::shared_ptr<SolidVelocityFunctionData<dim>> velocity;
  std::shared_ptr<Functions::ParsedFunction<dim>> alpha;
  std::shared_ptr<Functions::ParsedFunction<dim>> theta;

  SolidInitialConditionsParameters()
    : velocity(std::make_shared<SolidVelocityFunctionData<dim>>())
    , alpha(std::make_shared<Functions::ParsedFunction<dim>>(1))
    , theta(std::make_shared<Functions::ParsedFunction<dim>>(1))
  {}

  static void
  declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("velocity");
    SolidVelocityFunctionData<dim>::declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("alpha");
    {
      Functions::ParsedFunction<dim> tmp(1);
      tmp.declare_parameters(prm);
    }
    prm.leave_subsection();

    prm.enter_subsection("theta");
    {
      Functions::ParsedFunction<dim> tmp(1);
      tmp.declare_parameters(prm);
    }
    prm.leave_subsection();
  }

  void
  parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("velocity");
    velocity->parse_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("alpha");
    alpha->parse_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("theta");
    theta->parse_parameters(prm);
    prm.leave_subsection();
  }
};

template <int dim>
struct SolidSourceTermParameters
{
  bool enable = false;

  // Component order:
  // 2D: S_mx ; S_my ; S_alpha
  // 3D: S_mx ; S_my ; S_mz ; S_alpha
  std::shared_ptr<Functions::ParsedFunction<dim>> function;

  SolidSourceTermParameters()
    : function(std::make_shared<Functions::ParsedFunction<dim>>(dim + 1))
  {}

  static void
  declare_parameters(ParameterHandler &prm)
  {
    prm.declare_entry("enable", "false", Patterns::Bool());

    Functions::ParsedFunction<dim> tmp(dim + 1);
    tmp.declare_parameters(prm);
  }

  void
  parse_parameters(ParameterHandler &prm)
  {
    enable = prm.get_bool("enable");

    function = std::make_shared<Functions::ParsedFunction<dim>>(dim + 1);
    function->parse_parameters(prm);
  }
};

template <int dim>
struct SolidAnalyticalSolutionParameters
{
  bool enable = false;

  std::string  filename  = "solid_mms_error.dat";
  unsigned int frequency = 1;

  // Component order:
  // 2D: u_s ; v_s ; alpha
  // 3D: u_s ; v_s ; w_s ; alpha
  std::shared_ptr<Functions::ParsedFunction<dim>> function;

  SolidAnalyticalSolutionParameters()
    : function(std::make_shared<Functions::ParsedFunction<dim>>(dim + 1))
  {}

  static void
  declare_parameters(ParameterHandler &prm)
  {
    prm.declare_entry("enable", "false", Patterns::Bool());
    prm.declare_entry("filename", "solid_mms_error.dat", Patterns::Anything());
    prm.declare_entry("frequency", "1", Patterns::Integer(1));

    Functions::ParsedFunction<dim> tmp(dim + 1);
    tmp.declare_parameters(prm);
  }

  void
  parse_parameters(ParameterHandler &prm)
  {
    enable    = prm.get_bool("enable");
    filename  = prm.get("filename");
    frequency = prm.get_integer("frequency");

    function = std::make_shared<Functions::ParsedFunction<dim>>(dim + 1);
    function->parse_parameters(prm);
  }
};


// -----------------------------------------------------------------------------
// Boundary conditions
// -----------------------------------------------------------------------------

template <int dim>
struct SolidBoundaryConditionsParameters
{
  unsigned int number_of_boundary_conditions = 0;
  bool         time_dependent                = false;

  std::map<types::boundary_id, BoundaryConditions::BoundaryType> type;

  std::map<types::boundary_id, std::shared_ptr<SolidVelocityFunctionData<dim>>>
    velocity_functions;

  std::map<types::boundary_id, std::shared_ptr<Functions::ParsedFunction<dim>>>
    alpha_functions;

  std::map<types::boundary_id, std::shared_ptr<Functions::ParsedFunction<dim>>>
    theta_functions;

  std::map<types::boundary_id, types::boundary_id> periodic_neighbor_id;
  std::map<types::boundary_id, unsigned int>       periodic_direction;

  static void
  declare_default_entry(ParameterHandler        &prm,
                        const types::boundary_id default_id)
  {
    prm.declare_entry("type",
                      "noslip",
                      Patterns::Selection(
                        "noslip|slip|function|outlet|periodic"));
    prm.declare_entry("id",
                      Utilities::int_to_string(default_id, 2),
                      Patterns::Integer());
    prm.declare_entry("periodic_id", "-1", Patterns::Integer());
    prm.declare_entry("periodic_direction", "0", Patterns::Integer());

    prm.enter_subsection("velocity");
    SolidVelocityFunctionData<dim>::declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("alpha");
    {
      Functions::ParsedFunction<dim> tmp(1);
      tmp.declare_parameters(prm);
    }
    prm.leave_subsection();

    prm.enter_subsection("theta");
    {
      Functions::ParsedFunction<dim> tmp(1);
      tmp.declare_parameters(prm);
    }
    prm.leave_subsection();
  }

  static void
  declare_parameters(ParameterHandler &prm, unsigned int max_bc = 5)
  {
    prm.enter_subsection("Boundary conditions");
    {
      prm.declare_entry("number", "0", Patterns::Integer(0));
      prm.declare_entry("time dependent", "false", Patterns::Bool());

      for (unsigned int n = 0; n < max_bc; ++n)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          declare_default_entry(prm, n);
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  void
  parse_boundary(ParameterHandler &prm)
  {
    const types::boundary_id id       = prm.get_integer("id");
    const auto               type_str = prm.get("type");

    if (type_str == "noslip")
      type[id] = BoundaryConditions::BoundaryType::noslip;
    else if (type_str == "slip")
      type[id] = BoundaryConditions::BoundaryType::slip;
    else if (type_str == "outlet")
      type[id] = BoundaryConditions::BoundaryType::outlet;
    else if (type_str == "function")
      {
        type[id] = BoundaryConditions::BoundaryType::function;

        prm.enter_subsection("velocity");
        velocity_functions[id] =
          std::make_shared<SolidVelocityFunctionData<dim>>();
        velocity_functions[id]->parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("alpha");
        alpha_functions[id] =
          std::make_shared<Functions::ParsedFunction<dim>>(1);
        alpha_functions[id]->parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("theta");
        theta_functions[id] =
          std::make_shared<Functions::ParsedFunction<dim>>(1);
        theta_functions[id]->parse_parameters(prm);
        prm.leave_subsection();
      }
    else if (type_str == "periodic")
      {
        type[id]                 = BoundaryConditions::BoundaryType::periodic;
        periodic_neighbor_id[id] = prm.get_integer("periodic_id");
        periodic_direction[id]   = prm.get_integer("periodic_direction");
      }
    else
      AssertThrow(false, ExcMessage("Unsupported solid BC type."));
  }

  void
  parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Boundary conditions");
    {
      number_of_boundary_conditions = prm.get_integer("number");
      time_dependent                = prm.get_bool("time dependent");

      for (unsigned int n = 0; n < number_of_boundary_conditions; ++n)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          parse_boundary(prm);
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }
};


// -----------------------------------------------------------------------------
// Solid phase parameters
// -----------------------------------------------------------------------------

template <int dim>
struct SolidPhaseParameters
{
  unsigned int order = 1;

  double rho_s = 1.0;
  double beta  = 50.0;

  double         particle_diameter       = 1.0e-3;
  double         restitution_coefficient = 0.9;
  double         alpha_max               = 0.63;
  Tensor<1, dim> gravity;

  SolidInitialConditionsParameters<dim>  initial_conditions;
  SolidBoundaryConditionsParameters<dim> boundary_conditions;

  SolidSourceTermParameters<dim>         source_term;
  SolidAnalyticalSolutionParameters<dim> analytical_solution;

  // Nonlinear (Picard)
  unsigned int picard_max_iterations = 10;
  double       picard_tolerance      = 1e-8;
  double       picard_relaxation     = 0.3;

  // Output
  bool         output_verbose = true;
  unsigned int output_every   = 1;
  unsigned int digits         = 4;
  std::string  output_folder  = "output/";
  std::string  output_prefix  = "solid_phase";

  // Assembly / solver
  bool verbose_assembly = true;

  bool   use_alpha_supg    = true;
  double alpha_supg_factor = 10.0;

  bool   use_momentum_supg    = false;
  double momentum_supg_factor = 10.0;

  bool   use_solid_grad_div = true;
  double grad_div_factor    = 0.1;

  unsigned int quadrature_extra = 1;

  bool         solver_verbose = true;
  double       solver_abs_tol = 1e-12;
  double       solver_rel_tol = 1e-8;
  unsigned int solver_max_it  = 10000;
  unsigned int solver_restart = 100;

  bool print_timer_each_step = true;
  bool print_timer_at_end    = true;

  // Preconditioner
  unsigned int ilu_overlap       = 1;
  bool         amg_elliptic      = false;
  unsigned int amg_sweeps        = 2;
  double       amg_agg_threshold = 0.02;
  std::string  amg_smoother_type = "ILU";

  static void
  declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Solid phase");
    {
      prm.enter_subsection("Finite element");
      prm.declare_entry("order", "1", Patterns::Integer(1));
      prm.leave_subsection();

      prm.enter_subsection("Physics");
      prm.declare_entry("rho_s", "1.0", Patterns::Double(0.0));
      prm.declare_entry("beta", "5.0", Patterns::Double(0.0));
      prm.declare_entry("particle diameter", "1.0e-3", Patterns::Double(0.0));
      prm.declare_entry("restitution coefficient",
                        "0.9",
                        Patterns::Double(0.0, 1.0));
      prm.declare_entry("alpha max", "0.63", Patterns::Double(0.0, 1.0));
      prm.declare_entry("gravity",
                        (dim == 2 ? "0.0, 0.0" : "0.0, 0.0, 0.0"),
                        Patterns::List(Patterns::Double(), dim, dim));
      prm.leave_subsection();

      prm.enter_subsection("Initial conditions");
      SolidInitialConditionsParameters<dim>::declare_parameters(prm);
      prm.leave_subsection();

      SolidBoundaryConditionsParameters<dim>::declare_parameters(prm);

      prm.enter_subsection("Source term");
      SolidSourceTermParameters<dim>::declare_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("Analytical solution");
      SolidAnalyticalSolutionParameters<dim>::declare_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("Nonlinear");
      prm.declare_entry("picard max iterations", "10", Patterns::Integer(1));
      prm.declare_entry("picard tolerance", "1e-8", Patterns::Double(0.0));
      prm.declare_entry("picard relaxation", "0.3", Patterns::Double(0.0));
      prm.leave_subsection();

      prm.enter_subsection("Preconditioner");
      prm.declare_entry("ilu overlap", "1", Patterns::Integer(0));
      prm.declare_entry("amg elliptic", "false", Patterns::Bool());
      prm.declare_entry("amg sweeps", "2", Patterns::Integer(0));
      prm.declare_entry("amg aggregation threshold",
                        "0.02",
                        Patterns::Double(0.0));
      prm.declare_entry("amg smoother type", "ILU", Patterns::Anything());
      prm.leave_subsection();

      prm.enter_subsection("Assembly");
      prm.declare_entry("verbose", "true", Patterns::Bool());

      prm.declare_entry("use alpha supg", "true", Patterns::Bool());
      prm.declare_entry("alpha supg factor", "10.0", Patterns::Double(0.0));

      prm.declare_entry("use momentum supg", "false", Patterns::Bool());
      prm.declare_entry("momentum supg factor", "10.0", Patterns::Double(0.0));

      prm.declare_entry("use solid grad-div", "true", Patterns::Bool());
      prm.declare_entry("grad-div factor", "0.1", Patterns::Double(0.0));

      prm.declare_entry("quadrature extra", "1", Patterns::Integer(0));
      prm.leave_subsection();

      prm.enter_subsection("Timers");
      prm.declare_entry("print each step", "true", Patterns::Bool());
      prm.declare_entry("print at end", "true", Patterns::Bool());
      prm.leave_subsection();

      prm.enter_subsection("Solver");
      prm.declare_entry("verbose", "true", Patterns::Bool());
      prm.declare_entry("abs tol", "1e-12", Patterns::Double(0.0));
      prm.declare_entry("rel tol", "1e-8", Patterns::Double(0.0));
      prm.declare_entry("max it", "10000", Patterns::Integer(1));
      prm.declare_entry("restart", "100", Patterns::Integer(1));
      prm.leave_subsection();

      prm.enter_subsection("Output");
      prm.declare_entry("verbose", "true", Patterns::Bool());
      prm.declare_entry("every", "10", Patterns::Integer(1));
      prm.declare_entry("digits", "4", Patterns::Integer(1));
      prm.declare_entry("folder", "output/", Patterns::Anything());
      prm.declare_entry("prefix", "solid_phase", Patterns::Anything());
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  void
  parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Solid phase");
    {
      prm.enter_subsection("Finite element");
      order = prm.get_integer("order");
      prm.leave_subsection();

      prm.enter_subsection("Physics");
      rho_s                   = prm.get_double("rho_s");
      beta                    = prm.get_double("beta");
      particle_diameter       = prm.get_double("particle diameter");
      restitution_coefficient = prm.get_double("restitution coefficient");
      alpha_max               = prm.get_double("alpha max");

      const auto gravity_components =
        Utilities::split_string_list(prm.get("gravity"));
      AssertDimension(gravity_components.size(), dim);
      for (unsigned int d = 0; d < dim; ++d)
        gravity[d] = std::stod(gravity_components[d]);
      prm.leave_subsection();

      prm.enter_subsection("Initial conditions");
      initial_conditions.parse_parameters(prm);
      prm.leave_subsection();

      boundary_conditions.parse_parameters(prm);

      prm.enter_subsection("Source term");
      source_term.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("Analytical solution");
      analytical_solution.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("Nonlinear");
      picard_max_iterations = prm.get_integer("picard max iterations");
      picard_tolerance      = prm.get_double("picard tolerance");
      picard_relaxation     = prm.get_double("picard relaxation");
      prm.leave_subsection();

      prm.enter_subsection("Preconditioner");
      ilu_overlap       = prm.get_integer("ilu overlap");
      amg_elliptic      = prm.get_bool("amg elliptic");
      amg_sweeps        = prm.get_integer("amg sweeps");
      amg_agg_threshold = prm.get_double("amg aggregation threshold");
      amg_smoother_type = prm.get("amg smoother type");
      prm.leave_subsection();

      prm.enter_subsection("Assembly");
      verbose_assembly     = prm.get_bool("verbose");
      use_alpha_supg       = prm.get_bool("use alpha supg");
      alpha_supg_factor    = prm.get_double("alpha supg factor");
      use_momentum_supg    = prm.get_bool("use momentum supg");
      momentum_supg_factor = prm.get_double("momentum supg factor");
      use_solid_grad_div   = prm.get_bool("use solid grad-div");
      grad_div_factor      = prm.get_double("grad-div factor");
      quadrature_extra     = prm.get_integer("quadrature extra");
      prm.leave_subsection();

      prm.enter_subsection("Timers");
      print_timer_each_step = prm.get_bool("print each step");
      print_timer_at_end    = prm.get_bool("print at end");
      prm.leave_subsection();

      prm.enter_subsection("Solver");
      solver_verbose = prm.get_bool("verbose");
      solver_abs_tol = prm.get_double("abs tol");
      solver_rel_tol = prm.get_double("rel tol");
      solver_max_it  = prm.get_integer("max it");
      solver_restart = prm.get_integer("restart");
      prm.leave_subsection();

      prm.enter_subsection("Output");
      output_verbose = prm.get_bool("verbose");
      output_every   = prm.get_integer("every");
      digits         = prm.get_integer("digits");
      output_folder  = prm.get("folder");
      output_prefix  = prm.get("prefix");
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
};



template <int dim>
class SolidScratchData;

template <int dim>
struct SolidCopyData;



template <int dim>
class SolidTimeAssemblerBase
{
public:
  virtual ~SolidTimeAssemblerBase() = default;

  virtual void
  update_coefficients(const double dt) = 0;

  virtual void
  assemble(const SolidScratchData<dim> &scratch,
           SolidCopyData<dim>          &copy) const = 0;
};



template <int dim>
class SolidCoreAssembler;



template <int dim>
class SolidPhaseSolver
{
public:
  SolidPhaseSolver(const SolidPhaseParameters<dim>          &parameters,
                   const std::shared_ptr<SimulationControl> &simulation_control,
                   parallel::distributed::Triangulation<dim> &tria,
                   MPI_Comm comm = MPI_COMM_WORLD);

  void
  run();
  bool
  advance_one_step();
  void
  finalize();
  bool
  finished() const;

  void
  setup();

  unsigned int
  get_step_number() const;
  double
  get_current_time() const;

  const TrilinosWrappers::MPI::Vector &
  get_solid_volume_fraction() const;

  const TrilinosWrappers::MPI::Vector &
  get_granular_temperature() const;

  void
  set_fluid_solution_field(const DoFHandler<dim>               &fluid_dh,
                           const Mapping<dim>                  &fluid_mapping,
                           const TrilinosWrappers::MPI::Vector &fluid_solution);

  const DoFHandler<dim> &
  get_dof_handler() const
  {
    return dof_handler;
  }
  const FESystem<dim> &
  get_fe() const
  {
    return fe;
  }

  const TrilinosWrappers::MPI::BlockVector &
  get_locally_relevant_solution() const
  {
    return locally_relevant_solution;
  }

  double
  get_beta() const
  {
    return beta;
  }


private:
  void
  setup_dofs();
  void
  setup_granular_temperature_dofs();
  void
  assemble_system();
  void
  assemble_granular_temperature();
  void
  solve();
  void
  solve_granular_temperature();
  void
  output_results(const double time);

  void
  make_output_dir() const;

  void
  setup_time_assembler();
  void
  setup_core_assembler();
  void
  setup_linear_preconditioners();

  void
  solid_phase_conservation_monitoring() const;

  void
  calculate_mms_error() const;

  void
  assemble_local_system(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    SolidScratchData<dim>                                &scratch,
    SolidCopyData<dim>                                   &copy);

  void
  copy_local_system_to_global(const SolidCopyData<dim> &copy);


  const SolidPhaseParameters<dim> parameters;

  MPI_Comm mpi_communicator;

  std::shared_ptr<SimulationControl> simulation_control;

  const unsigned int  order;
  const FESystem<dim> fe;
  const FE_Q<dim>     theta_fe;

  parallel::distributed::Triangulation<dim> &triangulation;
  DoFHandler<dim>                            dof_handler;
  DoFHandler<dim>                            theta_dof_handler;

  AffineConstraints<double> constraints;
  AffineConstraints<double> theta_constraints;

  const double rho_s;
  double       beta;

  TrilinosWrappers::BlockSparseMatrix system_matrix;
  TrilinosWrappers::SparseMatrix      theta_matrix;

  TrilinosWrappers::MPI::BlockVector solution;
  TrilinosWrappers::MPI::BlockVector old_solution;
  TrilinosWrappers::MPI::BlockVector older_solution;
  TrilinosWrappers::MPI::BlockVector system_rhs;

  TrilinosWrappers::MPI::BlockVector locally_relevant_solution;
  TrilinosWrappers::MPI::BlockVector locally_relevant_old_solution;
  TrilinosWrappers::MPI::BlockVector locally_relevant_older_solution;

  TrilinosWrappers::MPI::Vector theta_solution;
  TrilinosWrappers::MPI::Vector theta_old_solution;
  TrilinosWrappers::MPI::Vector theta_picard_solution;
  TrilinosWrappers::MPI::Vector theta_rhs;

  TrilinosWrappers::MPI::Vector locally_relevant_theta_solution;
  TrilinosWrappers::MPI::Vector locally_relevant_theta_old_solution;
  TrilinosWrappers::MPI::Vector locally_relevant_theta_picard_solution;

  IndexSet theta_locally_owned;
  IndexSet theta_locally_relevant;

  IndexSet              locally_owned;
  IndexSet              locally_relevant;
  std::vector<IndexSet> owned_partitioning;
  std::vector<IndexSet> relevant_partitioning;

  TrilinosWrappers::MPI::BlockVector picard_solution;
  TrilinosWrappers::MPI::BlockVector locally_relevant_picard_solution;

  const TrilinosWrappers::MPI::BlockVector *picard_solution_ptr = nullptr;

  std::shared_ptr<SolidCoreAssembler<dim>>     core_assembler;
  std::shared_ptr<SolidTimeAssemblerBase<dim>> time_assembler;

  std::shared_ptr<TrilinosWrappers::PreconditionAMG> velocity_amg;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> velocity_ilu;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> alpha_ilu;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> theta_ilu;

  const DoFHandler<dim>               *fluid_dof_handler_ptr    = nullptr;
  const Mapping<dim>                  *fluid_mapping_ptr        = nullptr;
  const TrilinosWrappers::MPI::Vector *fluid_solution_ptr       = nullptr;
  bool                                 has_fluid_velocity_field = false;

  unsigned int timestep_number = 0;

  double cfl_length_scale = 1.0;

  ConditionalOStream pcout;
  TimerOutput        computing_timer;

  PVDHandler pvd_handler;

  unsigned int output_every;
  unsigned int digits;
  std::string  output_folder;
  std::string  output_prefix;
};

#endif