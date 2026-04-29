#ifndef lethe_solid_phase_h
#define lethe_solid_phase_h

#include <core/pvd_handler.h>

#include <solvers/fluid_dynamics_matrix_based.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_system.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>



struct SolidPhaseParameters
{
  // Mesh
  unsigned int nx = 10, ny = 10, nz = 10;
  unsigned int global_refinement = 0;
  std::string  direction1        = "x";
  std::string  direction2        = "y";

  // FE
  unsigned int degree = 1;

  // Time
  double       time_step = 1e-3;
  unsigned int n_steps   = 2000;

  // Physics
  double rho_s = 1.0;
  double beta  = 50.0;

  // IC/BC
  double alpha0      = 0.1;
  double alpha_inlet = 0.9;
  double u_inlet_x   = 1.0;
  double u_inlet_y   = 1.0;
  double u_inlet_z   = 0.0;

  // Output
  unsigned int output_every  = 10;
  unsigned int digits        = 4;
  std::string  output_folder = "output/";
  std::string  output_prefix = "solid_phase";



  bool verbose_assembly = true;


  bool         solver_verbose        = true;
  double       solver_abs_tol        = 1e-12;
  double       solver_rel_tol        = 1e-8;
  unsigned int solver_max_it         = 10000;
  unsigned int solver_restart        = 100;
  bool         output_verbose        = true;
  bool         print_timer_each_step = true;
  bool         print_timer_at_end    = true;
  unsigned int ilu_overlap           = 1;

  bool         amg_elliptic      = false;
  unsigned int amg_sweeps        = 2;
  double       amg_agg_threshold = 0.02;
  std::string  amg_smoother_type = "ILU";



  static void
  declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Solid phase");
    {
      prm.enter_subsection("Mesh");
      {
        prm.declare_entry("nx", "10", Patterns::Integer(1));
        prm.declare_entry("ny", "10", Patterns::Integer(1));
        prm.declare_entry("nz", "10", Patterns::Integer(1));
        prm.declare_entry("global refinement", "0", Patterns::Integer(0));
        prm.declare_entry("direction1", "x");
        prm.declare_entry("direction2", "y");
      }
      prm.leave_subsection();

      prm.enter_subsection("Finite element");
      {
        prm.declare_entry("degree", "1", Patterns::Integer(1));
      }
      prm.leave_subsection();

      prm.enter_subsection("Time");
      {
        prm.declare_entry("time step", "1e-3", Patterns::Double(0.0));
        prm.declare_entry("number of steps", "2000", Patterns::Integer(1));
      }
      prm.leave_subsection();

      prm.enter_subsection("Physics");
      {
        prm.declare_entry("rho_s", "1.0", Patterns::Double(0.0));
        prm.declare_entry("beta", "50.0", Patterns::Double(0.0));
      }
      prm.leave_subsection();

      prm.enter_subsection("Initial conditions");
      {
        prm.declare_entry("alpha0", "0.1", Patterns::Double(0.0, 1.0));
      }
      prm.leave_subsection();

      prm.enter_subsection("Preconditioner");
      {
        prm.declare_entry("ilu overlap", "1", Patterns::Integer(0));
        prm.declare_entry("amg elliptic", "false", Patterns::Bool());
        prm.declare_entry("amg sweeps", "2", Patterns::Integer(0));
        prm.declare_entry("amg aggregation threshold",
                          "0.02",
                          Patterns::Double(0.0));
        prm.declare_entry("amg smoother type", "ILU", Patterns::Anything());
      }
      prm.leave_subsection();

      prm.enter_subsection("Boundary conditions");
      {
        prm.declare_entry("alpha inlet", "0.9", Patterns::Double(0.0, 1.0));
        prm.declare_entry("u_inlet_x", "0.0");
        prm.declare_entry("u_inlet_y", "0.0");
        prm.declare_entry("u_inlet_z", "0.0");
      }
      prm.leave_subsection();

      prm.enter_subsection("Assembly");
      {
        prm.declare_entry("verbose", "true", Patterns::Bool());
      }
      prm.leave_subsection();

      prm.enter_subsection("Timers");
      {
        prm.declare_entry("print each step", "true", Patterns::Bool());
        prm.declare_entry("print at end", "true", Patterns::Bool());
      }
      prm.leave_subsection();


      prm.enter_subsection("Solver");
      {
        prm.declare_entry("verbose", "true", Patterns::Bool());
        prm.declare_entry("abs tol", "1e-12", Patterns::Double(0.0));
        prm.declare_entry("rel tol", "1e-8", Patterns::Double(0.0));
        prm.declare_entry("max it", "10000", Patterns::Integer(1));
        prm.declare_entry("restart", "100", Patterns::Integer(1));
      }
      prm.leave_subsection();

      prm.enter_subsection("Output");
      {
        prm.declare_entry("every", "10", Patterns::Integer(1));
        prm.declare_entry("digits", "4", Patterns::Integer(1));
        prm.declare_entry("folder", "output/", Patterns::Anything());
        prm.declare_entry("prefix", "solid_phase", Patterns::Anything());
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  void
  parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Solid phase");
    {
      prm.enter_subsection("Mesh");
      {
        nx                = prm.get_integer("nx");
        ny                = prm.get_integer("ny");
        nz                = prm.get_integer("nz");
        global_refinement = prm.get_integer("global refinement");
        direction1        = prm.get("direction1");
        direction2        = prm.get("direction2");
      }
      prm.leave_subsection();

      prm.enter_subsection("Finite element");
      {
        degree = prm.get_integer("degree");
      }
      prm.leave_subsection();

      prm.enter_subsection("Preconditioner");
      {
        ilu_overlap = prm.get_integer("ilu overlap");

        amg_elliptic      = prm.get_bool("amg elliptic");
        amg_sweeps        = prm.get_integer("amg sweeps");
        amg_agg_threshold = prm.get_double("amg aggregation threshold");
        amg_smoother_type = prm.get("amg smoother type");
      }
      prm.leave_subsection();


      prm.enter_subsection("Time");
      {
        time_step = prm.get_double("time step");
        n_steps   = prm.get_integer("number of steps");
      }
      prm.leave_subsection();

      prm.enter_subsection("Physics");
      {
        rho_s = prm.get_double("rho_s");
        beta  = prm.get_double("beta");
      }
      prm.leave_subsection();

      prm.enter_subsection("Initial conditions");
      {
        alpha0 = prm.get_double("alpha0");
      }
      prm.leave_subsection();

      prm.enter_subsection("Boundary conditions");
      {
        alpha_inlet = prm.get_double("alpha inlet");
        u_inlet_x   = prm.get_double("u_inlet_x");
        u_inlet_y   = prm.get_double("u_inlet_y");
        u_inlet_z   = prm.get_double("u_inlet_z");
      }
      prm.leave_subsection();

      prm.enter_subsection("Assembly");
      {
        verbose_assembly = prm.get_bool("verbose");
      }
      prm.leave_subsection();

      prm.enter_subsection("Timers");
      {
        print_timer_each_step = prm.get_bool("print each step");
        print_timer_at_end    = prm.get_bool("print at end");
      }
      prm.leave_subsection();


      prm.enter_subsection("Solver");
      {
        solver_verbose = prm.get_bool("verbose");
        solver_abs_tol = prm.get_double("abs tol");
        solver_rel_tol = prm.get_double("rel tol");
        solver_max_it  = prm.get_integer("max it");
        solver_restart = prm.get_integer("restart");
      }
      prm.leave_subsection();


      prm.enter_subsection("Output");
      {
        output_every  = prm.get_integer("every");
        digits        = prm.get_integer("digits");
        output_folder = prm.get("folder");
        output_prefix = prm.get("prefix");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
};

// template <int dim>
// class SolidPreconditioner : public Subscriptor
// {
// public:
//   std::shared_ptr<TrilinosWrappers::PreconditionAMG> amg_u;
//   std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_u;
//   std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_a;

//   enum class Mode
//   {
//     ilu,
//     amg
//   } mode = Mode::ilu;

//   void
//   vmult(TrilinosWrappers::MPI::BlockVector       &dst,
//         const TrilinosWrappers::MPI::BlockVector &src) const
//   {
//     if (dst.block(0).size() > 0)
//       {
//         if (mode == Mode::amg && amg_u)
//           {
//             amg_u->vmult(dst.block(0), src.block(0));
//           }
//         else if (ilu_u)
//           {
//             ilu_u->vmult(dst.block(0), src.block(0));
//           }
//         else
//           {
//             dst.block(0) = src.block(0);
//           }
//       }

//     if (dst.block(1).size() > 0)
//       {
//         if (ilu_a)
//           {
//             ilu_a->vmult(dst.block(1), src.block(1));
//           }
//         else
//           {
//             dst.block(1) = src.block(1);
//           }
//       }
//   }
// };

template <int dim>
class SolidPhaseSolver
{
public:
  SolidPhaseSolver(const SolidPhaseParameters                &parameters,
                   parallel::distributed::Triangulation<dim> &tria,
                   MPI_Comm comm = MPI_COMM_WORLD);


  void
  run();


  const TrilinosWrappers::MPI::Vector &
  get_solid_volume_fraction() const;

  // const TrilinosWrappers::MPI::Vector &
  // get_solid_velocity() const;

  void
  setup();

  bool
  advance_one_step();

  void
  finalize();

  bool
  finished() const;

  unsigned int
  get_step_number() const;

  double
  get_current_time() const;


  void
  set_fluid_velocity_field(const DoFHandler<dim>               &fluid_dh,
                           const Mapping<dim>                  &fluid_mapping,
                           const TrilinosWrappers::MPI::Vector &fluid_solution);

  const DoFHandler<dim> &
  get_dof_handler() const
  {
    return dof_handler;
  }

  const TrilinosWrappers::MPI::BlockVector &
  get_locally_relevant_solution() const
  {
    return locally_relevant_solution;
  }

  const FESystem<dim> &
  get_fe() const
  {
    return fe;
  }

  double
  get_beta() const
  {
    return beta;
  }
  

private:
  // void
  // make_grid();
  void
  setup_dofs();
  void
  assemble_system();
  void
  solve();
  void
  output_results(const double time);

  void
  make_output_dir() const;

  // void
  // setup_preconditioner();
  // void
  // setup_ILU();
  // void
  // setup_AMG();

  // void
  // update_constraints();

  void
  setup_linear_preconditioners();

  std::shared_ptr<TrilinosWrappers::PreconditionAMG> velocity_amg;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> velocity_ilu;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> alpha_ilu;


  const SolidPhaseParameters parameters;


  MPI_Comm mpi_communicator;


  std::vector<unsigned int> sub;


  const unsigned int  degree;
  const FESystem<dim> fe;

  parallel::distributed::Triangulation<dim> &triangulation;
  DoFHandler<dim>                            dof_handler;


  AffineConstraints<double> constraints;


  const double rho_s;
  double       beta;


  TrilinosWrappers::BlockSparseMatrix system_matrix;
  TrilinosWrappers::MPI::BlockVector  solution, old_solution, system_rhs;

  IndexSet              locally_owned;
  IndexSet              locally_relevant;
  std::vector<IndexSet> owned_partitioning;
  std::vector<IndexSet> relevant_partitioning;

  TrilinosWrappers::MPI::BlockVector locally_relevant_solution;
  TrilinosWrappers::MPI::BlockVector locally_relevant_old_solution;

  // std::shared_ptr<SolidPreconditioner<dim>> preconditioner;

  // time
  double time_step;
  // double       old_time_step;
  unsigned int timestep_number;

  // BC data
  Tensor<1, dim> inlet_velocity;

  // io/timers
  ConditionalOStream pcout;
  TimerOutput        computing_timer;

  // double max_inlet_velocity = 0.0;
  double cfl_length_scale = 1.0;


  // output
  PVDHandler   pvd_handler;
  unsigned int output_every;
  unsigned int group_files;
  unsigned int digits;
  std::string  output_folder;
  std::string  output_prefix;


  const DoFHandler<dim>               *fluid_dof_handler_ptr    = nullptr;
  const Mapping<dim>                  *fluid_mapping_ptr        = nullptr;
  const TrilinosWrappers::MPI::Vector *fluid_solution_ptr       = nullptr;
  bool                                 has_fluid_velocity_field = false;
};

#endif