#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/numerics/data_out.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <rpt/particle_detector_interactions.h>
#include <rpt/rpt_fem_reconstruction.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace dealii;


template <int dim>
RPTFEMReconstruction<dim>::AssemblyScratchData::AssemblyScratchData(
  const FiniteElement<dim> &fe,
  const unsigned int        no_detector)
  : fe_values(fe,
              QGaussSimplex<dim>(fe.degree + 1),
              update_values | update_quadrature_points | update_JxW_values)
  , detector_id(no_detector)
{}

template <int dim>
RPTFEMReconstruction<dim>::AssemblyScratchData::AssemblyScratchData(
  const AssemblyScratchData &scratch)
  : fe_values(scratch.fe_values.get_fe(),
              scratch.fe_values.get_quadrature(),
              update_values | update_quadrature_points | update_JxW_values)
  , detector_id(scratch.detector_id)
{}


template <int dim>
void
RPTFEMReconstruction<dim>::setup_triangulation()
{
  TimerOutput::Scope t(computing_timer, "setting_up_grid");

  if (rpt_parameters.fem_reconstruction_param.mesh_type ==
      Parameters::RPTFEMReconstructionParameters::FEMMeshType::gmsh)
    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream input_file(
        rpt_parameters.fem_reconstruction_param.mesh_file);
      grid_in.read_msh(input_file);

      const CylindricalManifold<dim> manifold(2);
      triangulation.set_all_manifold_ids(0);
      triangulation.set_manifold(0, manifold);
    }
  else
    {
      Triangulation<dim> temp_triangulation;
      Triangulation<dim> flat_temp_triangulation;

      GridGenerator::subdivided_cylinder(
        temp_triangulation,
        rpt_parameters.fem_reconstruction_param.z_subdivisions,
        rpt_parameters.rpt_param.reactor_radius,
        rpt_parameters.rpt_param.reactor_height * 0.5);
      temp_triangulation.refine_global(
        rpt_parameters.fem_reconstruction_param.mesh_refinement);

      // Flatten the triangulation
      GridGenerator::flatten_triangulation(temp_triangulation,
                                           flat_temp_triangulation);
      // Convert to simplex elements
      GridGenerator::convert_hypercube_to_simplex_mesh(flat_temp_triangulation,
                                                       triangulation);
      triangulation.set_all_manifold_ids(0);

      // Grid transformation
      Tensor<1, dim, double> axis({0, 1, 0});
      GridTools::rotate(axis, M_PI_2, triangulation);
      Tensor<1, dim> shift_vector(
        {0, 0, rpt_parameters.rpt_param.reactor_height * 0.5});
      GridTools::shift(shift_vector, triangulation);
    }
}

template <int dim>
void
RPTFEMReconstruction<dim>::setup_system()
{
  TimerOutput::Scope t(computing_timer, "setup_system");

  dof_handler.distribute_dofs(fe);

  system_rhs.reinit(dof_handler.n_dofs());
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);

  nodal_counts.resize(n_detector);
  for (auto &nodal_counts_for_one_detector : nodal_counts)
    {
      nodal_counts_for_one_detector.reinit(dof_handler.n_dofs());
    }

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();
}

template <int dim>
void
RPTFEMReconstruction<dim>::solve_linear_system(unsigned detector_no)
{
  TimerOutput::Scope t(computing_timer,
                       "solve_linear_system_" +
                         Utilities::to_string(detector_no, 2));

  SolverControl                          solver_control(1000, 1e-12);
  SolverCG<Vector<double>>               solver(solver_control);
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);
  solver.solve(system_matrix,
               nodal_counts[detector_no],
               system_rhs,
               preconditioner);
  constraints.distribute(nodal_counts[detector_no]);
}


template <int dim>
void
RPTFEMReconstruction<dim>::assemble_local_system(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  AssemblyScratchData &                                 sd,
  AssemblyCopyData &                                    copy_data)
{
  unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
  copy_data.cell_rhs.reinit(dofs_per_cell);

  copy_data.local_dof_indices.resize(dofs_per_cell);

  auto &fe_values = sd.fe_values;

  fe_values.reinit(cell);
  for (const unsigned int q_index : fe_values.quadrature_point_indices())
    {
      Point<dim>         q_point_position = fe_values.quadrature_point(q_index);
      RadioParticle<dim> particle(q_point_position, 0);

      ParticleDetectorInteractions<dim> p_q_interaction(
        particle, detectors[sd.detector_id], rpt_parameters.rpt_param);

      double count = p_q_interaction.calculate_count();

      for (const unsigned int i : fe_values.dof_indices())
        {
          for (const unsigned int j : fe_values.dof_indices())
            {
              copy_data.cell_matrix(i, j) +=
                (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                 fe_values.shape_value(j, q_index) * // phi_j(x_q)
                 fe_values.JxW(q_index));            // dx
            }
          copy_data.cell_rhs(i) +=
            (count *                             // f(x)
             fe_values.shape_value(i, q_index) * // phi_i(x_q)
             fe_values.JxW(q_index));            // dx
        }
    }
  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
RPTFEMReconstruction<dim>::copy_local_to_global(
  const AssemblyCopyData &copy_data)
{
  constraints.distribute_local_to_global(copy_data.cell_matrix,
                                         copy_data.cell_rhs,
                                         copy_data.local_dof_indices,
                                         system_matrix,
                                         system_rhs);
}

template <int dim>
void
RPTFEMReconstruction<dim>::assemble_system(unsigned no_detector)
{
  TimerOutput::Scope t(computing_timer,
                       "assemble_system_" +
                         Utilities::to_string(no_detector, 2));
  system_rhs    = 0;
  system_matrix = 0;

  WorkStream::run(dof_handler.begin_active(),
                  dof_handler.end(),
                  *this,
                  &RPTFEMReconstruction::assemble_local_system,
                  &RPTFEMReconstruction::copy_local_to_global,
                  AssemblyScratchData(fe, no_detector),
                  AssemblyCopyData());
}


template <int dim>
void
RPTFEMReconstruction<dim>::assign_detector_positions()
{
  TimerOutput::Scope t(computing_timer, "assigning_detector_positions");

  // Read text file with detector positions and store it in vector
  std::ifstream detector_file(
    rpt_parameters.detector_param.detector_positions_file);

  std::vector<double> values;
  std::copy(std::istream_iterator<double>(detector_file),
            std::istream_iterator<double>(),
            std::back_inserter(values));

  // Get the number of detector (2 positions for 1 detector, face and middle)
  int number_of_detector = values.size() / (2 * dim);

  // Extract positions, create point objects and detectors
  for (int i = 0; i < number_of_detector; i++)
    {
      Point<dim> face_point(values[2 * dim * i],
                            values[2 * dim * i + 1],
                            values[2 * dim * i + 2]);
      Point<dim> middle_point(values[2 * dim * i + dim],
                              values[2 * dim * i + dim + 1],
                              values[2 * dim * i + dim + 2]);

      Detector<dim> detector(rpt_parameters.detector_param,
                             i,
                             face_point,
                             middle_point);

      detectors.push_back(detector);
    }
}

template <int dim>
void
RPTFEMReconstruction<dim>::output_results()
{
  TimerOutput::Scope t(computing_timer, "output_results_vtu");

  // Export ".vtu" file with the nodal counts
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  for (unsigned d = 0; d < n_detector; ++d)
    {
      data_out.add_data_vector(nodal_counts[d],
                               "detector_" + Utilities::to_string(d, 2));
    }
  data_out.build_patches();
  std::ofstream output("solution.vtu");
  data_out.write_vtu(output);
}

template <int dim>
void
RPTFEMReconstruction<dim>::output_raw_results()
{
  TimerOutput::Scope t(computing_timer, "output_results_raw");

  std::map<types::global_dof_index, Point<dim>> dof_index_and_location;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      for (unsigned int v = 0; v < cell->n_vertices(); ++v)
        {
          auto dof_index       = cell->vertex_dof_index(v, 0);
          auto vertex_location = cell->vertex(v);
          std::pair<types::global_dof_index, Point<dim>> dof_and_vertex(
            dof_index, vertex_location);
          dof_index_and_location.insert(dof_and_vertex);
        }
    }

  std::string filename = "raw_counts.dat";

  // Open a file
  std::ofstream myfile;
  std::string   sep = " ";
  myfile.open(filename);
  myfile << "vertex_positions_x vertex_position_y vertex_position_z ";
  for (unsigned int i = 0; i < n_detector; ++i)
    myfile << "detector_" + Utilities::to_string(i, 2) + sep;
  myfile << std::endl;

  // Output in file and on terminal
  if (rpt_parameters.rpt_param.verbosity == Parameters::Verbosity::verbose)
    {
      for (auto it = dof_index_and_location.begin();
           it != dof_index_and_location.end();
           ++it)
        {
          for (unsigned p = 0; p < dim; ++p)
            myfile << it->second[p] << sep;

          for (unsigned int d = 0; d < n_detector; ++d)
            {
              myfile << nodal_counts[d][it->first] << sep;
              std::cout << nodal_counts[d][it->first] << sep;
            }

          myfile << "\n";
          std::cout << std::endl;
        }
    }
  else // Output in file only
    {
      for (auto it = dof_index_and_location.begin();
           it != dof_index_and_location.end();
           ++it)
        {
          for (unsigned p = 0; p < dim; ++p)
            myfile << it->second[p] << sep;

          for (unsigned int d = 0; d < n_detector; ++d)
            {
              myfile << nodal_counts[d][it->first] << sep;
            }

          myfile << "\n";
        }
    }
  myfile.close();
}


template <int dim>
void
RPTFEMReconstruction<dim>::L2_project()
{
  MultithreadInfo::set_thread_limit(1);
  std::cout << "***********************************************" << std::endl;
  std::cout << "Assigning detector positions" << std::endl;
  assign_detector_positions();

  n_detector = detectors.size();
  std::cout << "Number of detectors identified: " << n_detector << std::endl;
  std::cout << "***********************************************" << std::endl;
  std::cout << "Setting up the grid" << std::endl;
  setup_triangulation();
  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
  std::cout << "***********************************************" << std::endl;

  setup_system();

  for (unsigned d = 0; d < n_detector; ++d)
    {
      std::cout << "Detector_id: " << Utilities::to_string(d, 2) << std::endl;
      std::cout << "Assembling system" << std::endl;
      assemble_system(d);
      std::cout << "Solving system" << std::endl;
      solve_linear_system(d);
      std::cout << "System solved" << std::endl;
      std::cout << "-----------------------------------------------"
                << std::endl;
    }
  std::cout << "Outputting results" << std::endl;
  output_results();
  output_raw_results();
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "Saving dof handler and nodal counts" << std::endl;
  checkpoint();
  std::cout << "***********************************************" << std::endl;
  std::cout << "Done!" << std::endl;
  std::cout << "***********************************************" << std::endl;

  // Disable the output of time clock
  if (!rpt_parameters.fem_reconstruction_param.verbose_clock_fem_reconstruction)
    computing_timer.disable_output();
}


template <int dim>
Vector<double>
assemble_matrix_and_rhs(
  std::vector<std::vector<double>> &vertex_count,
  std::vector<double> &             experimental_count,
  Parameters::RPTFEMReconstructionParameters::FEMCostFunction
    &cost_function_type)
{
  Vector<double>           reference_location;
  unsigned int             detector_size = vertex_count.size();
  double                   sigma         = 0;
  LAPACKFullMatrix<double> sys_matrix;
  Vector<double>           sys_rhs;

  sys_matrix.reinit(dim);
  sys_rhs.reinit(dim);

  if (cost_function_type ==
      Parameters::RPTFEMReconstructionParameters::FEMCostFunction::absolute)
    {
      // Assembling sys_matrix
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int j = 0; j < dim; ++j)
            {
              for (unsigned int d = 0; d < detector_size; ++d)
                {
                  sigma += (-vertex_count[d][0] + vertex_count[d][j + 1]) *
                           (-vertex_count[d][0] + vertex_count[d][i + 1]);
                }
              sys_matrix.set(i, j, sigma);
              sigma = 0;
            }
        }

      // Assembling sys_rhs
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int d = 0; d < detector_size; ++d)
            {
              sigma += (vertex_count[d][0] - experimental_count[d]) *
                       (-vertex_count[d][0] + vertex_count[d][i + 1]);
            }
          sys_rhs[i] = -sigma;
          sigma      = 0;
        }
    }
  else if (cost_function_type == Parameters::RPTFEMReconstructionParameters::
                                   FEMCostFunction::relative)
    {
      std::vector<double> denom(detector_size);

      for (unsigned int d = 0; d < experimental_count.size(); ++d)
        denom[d] = 1 / (experimental_count[d] * experimental_count[d]);

      // Assembling sys_matrix
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int j = 0; j < dim; ++j)
            {
              for (unsigned int d = 0; d < detector_size; ++d)
                {
                  sigma += (-vertex_count[d][0] + vertex_count[d][j + 1]) *
                           (-vertex_count[d][0] + vertex_count[d][i + 1]) *
                           denom[d];
                }
              sys_matrix.set(i, j, sigma);
              sigma = 0;
            }
        }

      // Assembling sys_rhs
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int d = 0; d < detector_size; ++d)
            {
              sigma += (vertex_count[d][0] - experimental_count[d]) *
                       (-vertex_count[d][0] + vertex_count[d][i + 1]) *
                       denom[d];
            }
          sys_rhs[i] = -sigma;
          sigma      = 0;
        }
    }

  // Setup and solve linear system
  sys_matrix.set_property(LAPACKSupport::general);
  sys_matrix.compute_lu_factorization();
  sys_matrix.solve(sys_rhs);

  reference_location = sys_rhs;
  return reference_location;
}


template <int dim>
double
RPTFEMReconstruction<dim>::calculate_reference_location_error(
  Vector<double> &reference_location,
  const double &  last_constraint)
{
  std::vector<double> err_coordinates(dim + 1, 0);
  double              norm_error_coordinates = 0;

  // Check if the location is a valid one
  for (unsigned int i = 0; i < dim; ++i)
    {
      if (reference_location[i] > 1)
        err_coordinates[i] = (reference_location[i] - 1);

      if (reference_location[i] < 0)
        err_coordinates[i] = (0 - reference_location[i]);
    }

  // Check the last constraint
  if (last_constraint < 0)
    err_coordinates[3] = std::abs(last_constraint);
  else if (last_constraint > 1)
    err_coordinates[3] = 1 - last_constraint;

  // Calculate the norm of the error coordinate vector
  for (const double &error : err_coordinates)
    {
      norm_error_coordinates += error * error;
    }

  return norm_error_coordinates;
}


template <int dim>
double
RPTFEMReconstruction<dim>::calculate_cost(
  const TriaActiveIterator<DoFCellAccessor<dim, dim, false>> &cell,
  Vector<double> &     reference_location,
  const double &       last_constraint,
  std::vector<double> &experimental_count)
{
  double cost  = 0;
  double count = 0;

  if (rpt_parameters.fem_reconstruction_param.fem_cost_function ==
      Parameters::RPTFEMReconstructionParameters::FEMCostFunction::absolute)
    {
      for (unsigned int d = 0; d < n_detector; ++d)
        {
          count +=
            (nodal_counts[d][cell->vertex_dof_index(0, 0)] * last_constraint);

          for (unsigned int i = 0; i < dim; ++i)
            {
              count += (nodal_counts[d][cell->vertex_dof_index(i + 1, 0)] *
                        reference_location[i]);
            }
          count -= experimental_count[d];
          cost += count * count;
          count = 0;
        }
    }
  else if (rpt_parameters.fem_reconstruction_param.fem_cost_function ==
           Parameters::RPTFEMReconstructionParameters::FEMCostFunction::
             relative)
    {
      std::vector<double> denom(n_detector);

      for (unsigned int d = 0; d < n_detector; ++d)
        {
          count +=
            (nodal_counts[d][cell->vertex_dof_index(0, 0)] * last_constraint);

          denom[d] = 1 / (experimental_count[d] * experimental_count[d]);

          for (unsigned int i = 0; i < dim; ++i)
            {
              count += (nodal_counts[d][cell->vertex_dof_index(i + 1, 0)] *
                        reference_location[i]);
            }

          count -= experimental_count[d];
          cost += count * count * denom[d];
          count = 0;
        }
    }

  return cost;
}


template <int dim>
void
RPTFEMReconstruction<dim>::find_cell(std::vector<double> &experimental_count,
                                     const double tol_reference_location)
{
  double                           max_cost_function = DBL_MAX;
  double                           last_constraint_reference_location;
  double                           norm_error_coordinates;
  double                           calculated_cost;
  Point<dim>                       real_location;
  Vector<double>                   reference_location;
  std::vector<std::vector<double>> count_from_all_detectors(
    n_detector, std::vector<double>(4));

  // Loop over cell, loop over detectors and get nodal count values for each
  // vertex of the cell
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      for (unsigned int d = 0; d < n_detector; ++d)
        {
          for (unsigned int v = 0; v < cell->n_vertices(); ++v)
            {
              auto dof_index                 = cell->vertex_dof_index(v, 0);
              count_from_all_detectors[d][v] = nodal_counts[d][dof_index];
            }
        }

      // Solve linear system to find the location in reference coordinates
      reference_location = assemble_matrix_and_rhs<dim>(
        count_from_all_detectors,
        experimental_count,
        rpt_parameters.fem_reconstruction_param.fem_cost_function);

      // 4th constraint on the location of the particle in reference coordinates
      last_constraint_reference_location = 1 - reference_location[0] -
                                           reference_location[1] -
                                           reference_location[2];

      // Evaluate the error of the reference position (Is it outside the
      // reference tetrahedron ?)
      norm_error_coordinates =
        calculate_reference_location_error(reference_location,
                                           last_constraint_reference_location);

      // Extrapolation limit
      if (norm_error_coordinates < tol_reference_location)
        {
          // Calculate cost with the selected cost function
          calculated_cost = calculate_cost(cell,
                                           reference_location,
                                           last_constraint_reference_location,
                                           experimental_count);
          // Check if the new cost is lower than the previously stored
          // real_location's
          if (calculated_cost < max_cost_function)
            {
              max_cost_function      = calculated_cost;
              previous_position_cell = cell;

              // Evaluate the real location of the particle
              for (unsigned int v = 0; v < cell->n_vertices(); ++v)
                {
                  if (v == 0)
                    real_location = cell->vertex(v);
                  else
                    real_location += reference_location[v - 1] *
                                     (cell->vertex(v) - cell->vertex(0));
                }
            }
        }
    }
  found_positions.push_back(real_location);
}

template <int dim>
bool
RPTFEMReconstruction<dim>::find_in_adjacent_cells(
  std::vector<double> &                                 experimental_count,
  const double                                          tol_reference_location,
  const typename DoFHandler<dim>::active_cell_iterator &cell)
{
  double                           max_cost_function = DBL_MAX;
  double                           last_constraint_reference_location;
  double                           norm_error_coordinates;
  double                           calculated_cost;
  bool                             position_found = false;
  Point<dim>                       real_location;
  Vector<double>                   reference_location;
  std::vector<std::vector<double>> count_from_all_detectors(
    n_detector, std::vector<double>(4));

  // Find cells adjacent to the cell were found the previous position
  std::vector<typename DoFHandler<dim>::active_cell_iterator> adjacent_cells;
  std::set<typename DoFHandler<dim>::active_cell_iterator> all_adjacent_cells;

  // Add previous solution's cell
  all_adjacent_cells.insert(cell);

  // Search for adjacent cells and stock them in the set container
  for (unsigned int v = 0; v < cell->n_vertices(); ++v)
    {
      auto v_index = cell->vertex_index(v);
      adjacent_cells =
        GridTools::find_cells_adjacent_to_vertex(dof_handler, v_index);

      for (const auto &adjacent_cell : adjacent_cells)
        all_adjacent_cells.insert(adjacent_cell);
    }

  // Add adjacent cells according to the level of proximity
  if (rpt_parameters.fem_reconstruction_param.search_proximity_level > 1)
    {
      std::set<typename DoFHandler<dim>::active_cell_iterator>
        previous_adjacent_cells;
      std::set<typename DoFHandler<dim>::active_cell_iterator>
        previous_main_cells;

      previous_main_cells.insert(cell);

      for (unsigned int proximity_level = 2;
           proximity_level <
           rpt_parameters.fem_reconstruction_param.search_proximity_level + 1;
           ++proximity_level)
        {
          previous_adjacent_cells = all_adjacent_cells;

          for (const auto &current_main_cell : previous_adjacent_cells)
            {
              if (*previous_main_cells.find(current_main_cell) !=
                  current_main_cell)
                {
                  for (unsigned int v = 0; v < current_main_cell->n_vertices();
                       ++v)
                    {
                      auto v_index = current_main_cell->vertex_index(v);
                      adjacent_cells =
                        GridTools::find_cells_adjacent_to_vertex(dof_handler,
                                                                 v_index);

                      for (const auto &adjacent_cell : adjacent_cells)
                        all_adjacent_cells.insert(adjacent_cell);
                    }
                }
              previous_main_cells = previous_adjacent_cells;
            }
        }
    }

  // Loop over adjacent cells, loop over detectors and get nodal count values
  // for each vertex of the cell
  for (const auto &adjacent_cell : all_adjacent_cells)
    {
      for (unsigned int d = 0; d < n_detector; ++d)
        {
          for (unsigned int v = 0; v < adjacent_cell->n_vertices(); ++v)
            {
              auto dof_index = adjacent_cell->vertex_dof_index(v, 0);
              count_from_all_detectors[d][v] = nodal_counts[d][dof_index];
            }
        }

      // Solve linear system to find the location in reference coordinates
      reference_location = assemble_matrix_and_rhs<dim>(
        count_from_all_detectors,
        experimental_count,
        rpt_parameters.fem_reconstruction_param.fem_cost_function);

      // 4th constraint on the location of the particle in reference coordinates
      last_constraint_reference_location = 1 - reference_location[0] -
                                           reference_location[1] -
                                           reference_location[2];

      // Evaluate the error of the reference position (Is it outside the
      // reference tetrahedron ?)
      norm_error_coordinates =
        calculate_reference_location_error(reference_location,
                                           last_constraint_reference_location);

      // Extrapolation limit
      if (norm_error_coordinates < tol_reference_location)
        {
          // Calculate cost with the selected cost function
          calculated_cost = calculate_cost(adjacent_cell,
                                           reference_location,
                                           last_constraint_reference_location,
                                           experimental_count);
          // Check if the new cost is lower than the previously stored
          // real_location's
          if (calculated_cost < max_cost_function)
            {
              position_found         = true;
              max_cost_function      = calculated_cost;
              previous_position_cell = adjacent_cell;

              // Evaluate the real location of the particle
              for (unsigned int v = 0; v < adjacent_cell->n_vertices(); ++v)
                {
                  if (v == 0)
                    real_location = adjacent_cell->vertex(v);
                  else
                    real_location +=
                      reference_location[v - 1] *
                      (adjacent_cell->vertex(v) - adjacent_cell->vertex(0));
                }
            }
        }
    }

  if (position_found)
    {
      found_positions.push_back(real_location);
      return true;
    }
  else
    return false;
}

template <int dim>
void
RPTFEMReconstruction<dim>::read_experimental_counts(
  std::vector<std::vector<double>> &all_experimental_counts)
{
  TimerOutput::Scope t(computing_timer, "read_experimental_counts");

  std::ifstream in(
    rpt_parameters.fem_reconstruction_param.experimental_counts_file);
  double value;

  if (in)
    {
      std::string line;

      while (std::getline(in, line))
        {
          // Ignore empty lines
          if (line == "")
            continue;

          all_experimental_counts.emplace_back(std::vector<double>());

          // Break down the row into column values
          std::stringstream split(line);

          while (split >> value)
            all_experimental_counts.back().push_back(value);
        }
    }
}


template <int dim>
void
RPTFEMReconstruction<dim>::trajectory()
{
  // Tolerance or the extrapolation limit in the reference space for a found
  // position
  double tol_reference_location = 0.005;

  if (rpt_parameters.fem_reconstruction_param.mesh_type ==
      Parameters::RPTFEMReconstructionParameters::FEMMeshType::dealii)
    {
      const unsigned int power =
        pow(2, rpt_parameters.fem_reconstruction_param.mesh_refinement);
      const unsigned int n_cell_z =
        2 * rpt_parameters.fem_reconstruction_param.z_subdivisions * power;
      tol_reference_location =
        rpt_parameters.rpt_param.reactor_height / n_cell_z * 1.15;
    }

  std::cout << "tol: " << tol_reference_location << std::endl;

  // Read and store all experimental counts
  std::vector<std::vector<double>> all_experimental_counts;
  read_experimental_counts(all_experimental_counts);

  {
    TimerOutput::Scope t(computing_timer, "find_particle_positions");

    if (rpt_parameters.fem_reconstruction_param.search_type ==
        Parameters::RPTFEMReconstructionParameters::FEMSearchType::local)
      {
        bool adjacent_cell_search = false;

        // Find the position of the particle with the experimental counts by
        // prioritizing local search
        for (std::vector<double> &experimental_counts : all_experimental_counts)
          {
            if (adjacent_cell_search)
              adjacent_cell_search =
                find_in_adjacent_cells(experimental_counts,
                                       tol_reference_location,
                                       previous_position_cell);

            if (!adjacent_cell_search)
              {
                find_cell(experimental_counts, tol_reference_location);
                adjacent_cell_search = true;
              }
          }
      }
    else
      {
        // Find the position of the particle with the experimental counts by
        // global search only
        for (std::vector<double> &experimental_counts : all_experimental_counts)
          {
            find_cell(experimental_counts, tol_reference_location);
          }
      }
  }
}

template <int dim>
void
RPTFEMReconstruction<dim>::checkpoint()
{
  TimerOutput::Scope t(computing_timer, "checkpoint");

  // Save dof_handler object
  {
    std::ofstream                 ofs("temp_dof_handler.dof");
    boost::archive::text_oarchive oa(ofs);
    dof_handler.save(oa, 0);
  }

  // Save nodal_counts_per_detector
  {
    for (unsigned int i = 0; i < n_detector; ++i)
      {
        std::string                   filename("temp_nodal_counts_detector" +
                             Utilities::to_string(i, 2) + ".counts");
        std::ofstream                 ofs(filename);
        boost::archive::text_oarchive oa(ofs);
        nodal_counts[i].save(oa, 0);
      }
  }
}


template <int dim>
void
RPTFEMReconstruction<dim>::load_from_checkpoint()
{
  TimerOutput::Scope t(computing_timer, "load_from_checkpoint");

  n_detector = rpt_parameters.fem_reconstruction_param.nodal_counts_file.size();

  // Import dof handler
  {
    dof_handler.distribute_dofs(fe);
    std::ifstream ifs(rpt_parameters.fem_reconstruction_param.dof_handler_file);
    boost::archive::text_iarchive ia(ifs);
    dof_handler.load(ia, 0);
  }
  // Import nodal counts
  {
    Vector<double> counts_per_detector;

    nodal_counts.resize(n_detector);

    for (unsigned int i = 0; i < n_detector; ++i)
      {
        std::ifstream ifs(
          rpt_parameters.fem_reconstruction_param.nodal_counts_file[i]);
        boost::archive::text_iarchive ia(ifs);
        counts_per_detector.load(ia, 0);
        nodal_counts[i] = counts_per_detector;
      }
  }
}

template <int dim>
void
RPTFEMReconstruction<dim>::export_found_positions()
{
  TimerOutput::Scope t(computing_timer, "export_found_positions");

  std::string filename =
    rpt_parameters.fem_reconstruction_param.export_positions_file;

  // Look if extension of the exporting file is specified. If it is not
  // specified, the ".csv" extension is added as a default format.
  std::size_t csv_file = filename.find(".csv");
  std::size_t dat_file = filename.find(".dat");

  if ((csv_file == std::string::npos) && (dat_file == std::string::npos))
    filename += ".csv";

  // Open a file
  std::ofstream myfile;
  myfile.open(filename);

  // Output in file and in terminal
  if (rpt_parameters.rpt_param.verbosity == Parameters::Verbosity::verbose)
    {
      if (filename.substr(filename.find_last_of(".") + 1) == ".dat")
        {
          myfile << "position_x position_y position_z " << std::endl;
          for (const Point<dim> &position : found_positions)
            {
              myfile << position << std::endl;
              std::cout << position << std::endl;
            }
        }
      else
        {
          myfile << "position_x,position_y,position_z " << std::endl;
          std::string sep = ",";

          for (const Point<dim> &position : found_positions)
            {
              for (unsigned int i = 0; i < dim; ++i)
                myfile << position[i] << sep;
              myfile << std::endl;
              std::cout << position << std::endl;
            }
        }
    }
  else // Output only in file
    {
      if (filename.substr(filename.find_last_of(".") + 1) == ".dat")
        {
          myfile << "position_x position_y position_z " << std::endl;
          for (const Point<dim> &position : found_positions)
            {
              myfile << position << std::endl;
            }
        }
      else
        {
          myfile << "position_x,position_y,position_z " << std::endl;
          std::string sep = ",";

          for (const Point<dim> &position : found_positions)
            {
              for (unsigned int i = 0; i < dim; ++i)
                myfile << position[i] << sep;
              myfile << std::endl;
            }
        }
    }

  myfile.close();
}

template <int dim>
void
RPTFEMReconstruction<dim>::rpt_fem_reconstruct()
{
  std::cout << "***********************************************" << std::endl;
  std::cout << "Setting up the grid" << std::endl;
  setup_triangulation();
  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "Loading dof handler and nodal counts from " << std::endl;
  std::cout << "saved files " << std::endl;
  load_from_checkpoint();
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "Finding particle positions " << std::endl;
  trajectory();
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "Exporting particle positions " << std::endl;
  export_found_positions();
  std::cout << "***********************************************" << std::endl;
  std::cout << "Done!" << std::endl;
  std::cout << "***********************************************" << std::endl;

  // Disable the output of time clock
  if (!rpt_parameters.fem_reconstruction_param.verbose_clock_fem_reconstruction)
    computing_timer.disable_output();
}


template class RPTFEMReconstruction<3>;
