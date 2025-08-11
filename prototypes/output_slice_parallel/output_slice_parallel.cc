// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// This prototype tests the DataOutResample feature of deal.II
// It allows one to output the solution on a different triangulation than the
// one used for the simulation. This is useful for example when one is
// interested on looking at a slice of the domain.

#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_resample.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

using namespace dealii;

// Define all the parameters that can be specified in the .prm file
struct Settings
{
  bool
  try_parse(const std::string &prm_filename);

  enum GeometryType
  {
    hyperball,
    hypercube,
    hyperrectangle
  };

  GeometryType geometry;
  GeometryType slice_plane = hyperrectangle;

  unsigned int dimension;
  unsigned int element_order;
  unsigned int initial_refinement;
  unsigned int repetitions;
  bool         output;
  bool         output_slice;
  std::string  output_name;
  std::string  output_path;
};

bool
Settings::try_parse(const std::string &prm_filename)
{
  ParameterHandler prm;
  prm.declare_entry(
    "dim",
    "3",
    Patterns::Integer(),
    "The problem dimension <3>. This program only works in 3d due to the slice geometry used");
  prm.declare_entry("element order",
                    "1",
                    Patterns::Integer(),
                    "Order of FE element <1|2|3>");
  prm.declare_entry("geometry",
                    "hyperball",
                    Patterns::Selection("hyperball|hypercube|hyperrectangle"),
                    "Geometry <hyperball|hypercube|hyperrectangle>");
  prm.declare_entry("initial refinement",
                    "1",
                    Patterns::Integer(),
                    "Global refinement 1st cycle");
  prm.declare_entry("repetitions",
                    "1",
                    Patterns::Integer(),
                    "Repetitions in z direction for the hyper rectangle");
  prm.declare_entry("output",
                    "true",
                    Patterns::Bool(),
                    "Output vtu files <true|false>");
  prm.declare_entry("output slice",
                    "true",
                    Patterns::Bool(),
                    "Output slice vtu files <true|false>");
  prm.declare_entry("output name",
                    "solution",
                    Patterns::FileName(),
                    "Name for vtu files");
  prm.declare_entry("output path",
                    "./",
                    Patterns::FileName(),
                    "Path for vtu output files");

  if (prm_filename.size() == 0)
    {
      std::cout
        << "****  Error: No input file provided!\n"
        << "****  Error: Call this program as './output_slice_parallel input.prm\n"
        << '\n'
        << "****  You may want to use one of the input files in this\n"
        << "****  directory, or use the following default values\n"
        << "****  to create an input file:\n";
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
#if DEAL_II_VERSION_GTE(9, 7, 0)
        prm.print_parameters(std::cout, ParameterHandler::DefaultStyle);
#else
        prm.print_parameters(std::cout, ParameterHandler::Text);
#endif
      return false;
    }

  try
    {
      prm.parse_input(prm_filename);
    }
  catch (std::exception &e)
    {
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        std::cerr << e.what() << std::endl;
      return false;
    }

  if (prm.get("geometry") == "hyperball")
    this->geometry = hyperball;
  else if (prm.get("geometry") == "hypercube")
    this->geometry = hypercube;
  else if (prm.get("geometry") == "hyperrectangle")
    this->geometry = hyperrectangle;
  else
    AssertThrow(false, ExcNotImplemented());

  this->dimension = prm.get_integer("dim");

  AssertThrow(this->dimension == 3,
              ExcMessage("This program only works in 3d."));

  this->element_order      = prm.get_integer("element order");
  this->initial_refinement = prm.get_integer("initial refinement");
  this->repetitions        = prm.get_integer("repetitions");
  this->output             = prm.get_bool("output");
  this->output_slice       = prm.get_bool("output slice");
  this->output_name        = prm.get("output name");
  this->output_path        = prm.get("output path");

  return true;
}

template <int dim>
class SolutionFunction : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p,
        const unsigned int /* component */ = 0) const override
  {
    double val = 0.;
    for (unsigned int d = 0; d < dim; ++d)
      {
        val += std::pow(p[d], 2.);
      }
    return val;
  }
};

template <int dim, int fe_degree>
class OutputSliceParallel
{
public:
  OutputSliceParallel(const Settings &parameters);

  void
  run();

private:
  void
  make_grid();

  void
  make_slice_grid();

  void
  setup_system();

  void
  output_results();

  void
  output_slice_results();

  using VectorType = TrilinosWrappers::MPI::Vector;
  using MatrixType = TrilinosWrappers::SparseMatrix;

  parallel::distributed::Triangulation<dim>          triangulation;
  parallel::distributed::Triangulation<dim - 1, dim> patch_triangulation;
  const MappingQ<dim>                                mapping;
  const MappingQ<dim - 1, dim>                       patch_mapping;

  FE_Q<dim>       fe;
  DoFHandler<dim> dof_handler;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  VectorType solution;

  ConditionalOStream pcout;
  TimerOutput        computing_timer;
  MPI_Comm           mpi_communicator;

  Settings parameters;
};

template <int dim, int fe_degree>
OutputSliceParallel<dim, fe_degree>::OutputSliceParallel(
  const Settings &parameters)
  : triangulation(MPI_COMM_WORLD,
                  Triangulation<dim>::limit_level_difference_at_vertices,
                  parallel::distributed::Triangulation<dim>::default_setting)
  , patch_triangulation(MPI_COMM_WORLD)
  , mapping(fe_degree)
  , patch_mapping(fe_degree)
  , fe(fe_degree)
  , dof_handler(triangulation)
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , computing_timer(MPI_COMM_WORLD,
                    pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
  , mpi_communicator(MPI_COMM_WORLD)
  , parameters(parameters)
{}

template <int dim, int fe_degree>
void
OutputSliceParallel<dim, fe_degree>::make_grid()
{
  TimerOutput::Scope t(computing_timer, "make grid");

  switch (parameters.geometry)
    {
      case Settings::hyperball:
        {
          SphericalManifold<dim>                boundary_manifold;
          TransfiniteInterpolationManifold<dim> inner_manifold;

          GridGenerator::hyper_ball(triangulation);

          triangulation.set_all_manifold_ids(1);
          triangulation.set_all_manifold_ids_on_boundary(0);

          triangulation.set_manifold(0, boundary_manifold);

          inner_manifold.initialize(triangulation);
          triangulation.set_manifold(1, inner_manifold);

          break;
        }
      case Settings::hypercube:
        {
          GridGenerator::hyper_cube(triangulation);
          break;
        }
      case Settings::hyperrectangle:
        {
          std::vector<unsigned int> repetitions(dim);
          for (unsigned int i = 0; i < dim - 1; i++)
            {
              repetitions[i] = 1;
            }
          repetitions[dim - 1] = parameters.repetitions;

          GridGenerator::subdivided_hyper_rectangle(
            triangulation,
            repetitions,
            Point<dim>(-1, -1, -1.),
            Point<dim>(1., 1., parameters.repetitions - 1));
          break;
        }
    }
  triangulation.refine_global(parameters.initial_refinement);
}

template <int dim, int fe_degree>
void
OutputSliceParallel<dim, fe_degree>::make_slice_grid()
{
  TimerOutput::Scope t(computing_timer, "make slice grid");
  GridGenerator::subdivided_hyper_rectangle(patch_triangulation,
                                            std::vector<unsigned int>(dim, 1),
                                            Point<dim - 1>(-1., -1.),
                                            Point<dim - 1>(1., 1.));
  patch_triangulation.refine_global(parameters.initial_refinement);
}

template <int dim, int fe_degree>
void
OutputSliceParallel<dim, fe_degree>::setup_system()
{
  TimerOutput::Scope t(computing_timer, "setup system");

  dof_handler.distribute_dofs(fe);

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
}

template <int dim, int fe_degree>
void
OutputSliceParallel<dim, fe_degree>::output_results()
{
  TimerOutput::Scope t(computing_timer, "output complete");

  solution.update_ghost_values();

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation({DataComponentInterpretation::component_is_scalar});

  data_out.add_data_vector(solution,
                           "solution",
                           DataOut<dim>::type_dof_data,
                           interpretation);

  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    {
      subdomain(i) = triangulation.locally_owned_subdomain();
    }
  data_out.add_data_vector(subdomain,
                           "subdomain",
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.build_patches(mapping, fe.degree, DataOut<dim>::curved_inner_cells);

  DataOutBase::VtkFlags flags;
  data_out.set_flags(flags);
  data_out.write_vtu_in_parallel(parameters.output_path +
                                   parameters.output_name + ".vtu",
                                 MPI_COMM_WORLD);
}

template <int dim, int fe_degree>
void
OutputSliceParallel<dim, fe_degree>::output_slice_results()
{
  TimerOutput::Scope t(computing_timer, "output slice");

  solution.update_ghost_values();

  DataOutResample<dim, dim - 1, dim> data_out(patch_triangulation,
                                              patch_mapping);
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(dof_handler, solution, "solution");

  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    {
      subdomain(i) = triangulation.locally_owned_subdomain();
    }
  data_out.add_data_vector(subdomain, "subdomain");

  data_out.build_patches(mapping,
                         fe.degree,
                         DataOut<dim - 1, dim>::curved_inner_cells);
  DataOutBase::VtkFlags flags;
  data_out.set_flags(flags);
  data_out.write_vtu_in_parallel(parameters.output_path +
                                   parameters.output_name + "_slice.vtu",
                                 MPI_COMM_WORLD);
}
template <int dim, int fe_degree>
void
OutputSliceParallel<dim, fe_degree>::run()
{
  const unsigned int n_ranks = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  std::string DAT_header = "START DATE: " + Utilities::System::get_date() +
                           ", TIME: " + Utilities::System::get_time();
  std::string MPI_header = "Running with " + std::to_string(n_ranks) +
                           " MPI process" + (n_ranks > 1 ? "es" : "");
  std::string REFINE_header =
    "Initial refinement: " + std::to_string(parameters.initial_refinement);
  std::string REPETITIONS_header =
    "Repetitions in z: " + std::to_string(parameters.repetitions);
  pcout << std::string(80, '=') << std::endl;

  pcout << std::string(80, '=') << std::endl;
  pcout << DAT_header << std::endl;
  pcout << std::string(80, '-') << std::endl;
  pcout << MPI_header << std::endl;
  pcout << REFINE_header << std::endl;
  pcout << REPETITIONS_header << std::endl;

  Timer timer;

  make_grid();
  setup_system();

  pcout << "   Triangulation: " << triangulation.n_global_active_cells()
        << " cells" << std::endl;
  pcout << "   DoFHandler:    " << dof_handler.n_dofs() << " DoFs" << std::endl;

  SolutionFunction<dim> solution_function;
  VectorTools::interpolate(dof_handler, solution_function, solution);

  if (parameters.output)
    output_results();

  if (parameters.output_slice)
    {
      make_slice_grid();
      pcout << "   Triangulation for slice: "
            << patch_triangulation.n_global_active_cells() << " cells"
            << std::endl;
      output_slice_results();
    }

  timer.stop();
  computing_timer.print_summary();
  computing_timer.reset();
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize       mpi_init(argc, argv, 1);
  dealii::Utilities::System::MemoryStats stats;
  dealii::Utilities::System::get_memory_stats(stats);

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  Settings parameters;
  if (!parameters.try_parse((argc > 1) ? (argv[1]) : ""))
    return 0;
  try
    {
      switch (parameters.element_order)
        {
          case 1:
            {
              OutputSliceParallel<3, 1> output_slice_parallel(parameters);
              output_slice_parallel.run();

              break;
            }
          case 2:
            {
              OutputSliceParallel<3, 2> output_slice_parallel(parameters);
              output_slice_parallel.run();

              break;
            }
          case 3:
            {
              OutputSliceParallel<3, 3> output_slice_parallel(parameters);
              output_slice_parallel.run();

              break;
            }
        }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
