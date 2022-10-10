#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>
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

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace dealii;

template <int dim>
void
RPTL2Projection<dim>::setup_triangulation()
{
  if (fem_reconstruction_parameters.mesh_type ==
      Parameters::RPTFEMReconstructionParameters::FEMMeshType::gmsh)
    {
      // Warn user
      pcout << "-----------------------------------------------" << std::endl;
      pcout << "At the moment, the application calculates" << std::endl;
      pcout << "the gamma-ray attenuation in cylindrical" << std::endl;
      pcout << "vessels only. Using another type of geometry" << std::endl;
      pcout << "would yield in committing an error when" << std::endl;
      pcout << "calculating photon counts. Furthermore, the" << std::endl;
      pcout << "axis of symmetry of the cylinder should" << std::endl;
      pcout << "be oriented in z-direction." << std::endl;
      pcout << "-----------------------------------------------" << std::endl;

      // Import and attach mesh to triangulation
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream input_file(fem_reconstruction_parameters.mesh_file);
      grid_in.read_msh(input_file);

      // Attach a cylindrical manifold object to the triangulation
      const CylindricalManifold<dim> manifold(2); // Cylinder along the z-axis
      triangulation.set_all_manifold_ids(0);
      triangulation.set_manifold(0, manifold);

      // Check if all cells are tetrahedral elements
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          AssertThrow(
            cell->n_vertices() == 4,
            ExcMessage(
              "The imported mesh has to have only tetrahedral elements"));
        }
    }
  else
    {
      Triangulation<dim> temp_triangulation;
      Triangulation<dim> flat_temp_triangulation;

      GridGenerator::subdivided_cylinder(
        temp_triangulation,
        fem_reconstruction_parameters.z_subdivisions,
        parameters.reactor_radius,
        parameters.reactor_height * 0.5);
      temp_triangulation.refine_global(
        fem_reconstruction_parameters.mesh_refinement);

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
      Tensor<1, dim> shift_vector({0, 0, parameters.reactor_height * 0.5});
      GridTools::shift(shift_vector, triangulation);
    }
}

template <int dim>
void
RPTFEMReconstruction<dim>::setup_triangulation()
{
  TimerOutput::Scope t(computing_timer, "setting_up_grid");

  triangulation.clear();

  if (fem_reconstruction_parameters.mesh_type ==
      Parameters::RPTFEMReconstructionParameters::FEMMeshType::gmsh)
    {
      // Warn user
      std::cout << "-----------------------------------------------"
                << std::endl;
      std::cout << "At the moment, the application calculates" << std::endl;
      std::cout << "the gamma-ray attenuation in cylindrical" << std::endl;
      std::cout << "vessels only. Using another type of geometry" << std::endl;
      std::cout << "would yield in committing an error when" << std::endl;
      std::cout << "calculating photon counts. Furthermore, the" << std::endl;
      std::cout << "axis of symmetry of the cylinder should" << std::endl;
      std::cout << "be oriented in z-direction." << std::endl;
      std::cout << "-----------------------------------------------"
                << std::endl;

      // Import and attach mesh to triangulation
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream input_file(fem_reconstruction_parameters.mesh_file);
      grid_in.read_msh(input_file);

      // Attach a cylindrical manifold object to the triangulation
      const CylindricalManifold<dim> manifold(2); // Cylinder along the z-axis
      triangulation.set_all_manifold_ids(0);
      triangulation.set_manifold(0, manifold);

      // Check if all cells are tetrahedral elements
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          AssertThrow(
            cell->n_vertices() == 4,
            ExcMessage(
              "The imported mesh has to have only tetrahedral elements"));
        }
    }
  else
    {
      Triangulation<dim> temp_triangulation;
      Triangulation<dim> flat_temp_triangulation;

      GridGenerator::subdivided_cylinder(
        temp_triangulation,
        fem_reconstruction_parameters.z_subdivisions,
        parameters.reactor_radius,
        parameters.reactor_height * 0.5);
      temp_triangulation.refine_global(
        fem_reconstruction_parameters.mesh_refinement);

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
      Tensor<1, dim> shift_vector({0, 0, parameters.reactor_height * 0.5});
      GridTools::shift(shift_vector, triangulation);
    }
}

template <int dim>
void
RPTL2Projection<dim>::setup_system()
{
  TimerOutput::Scope t(computing_timer, "setup_system");

  // Partition triangulation between the processes
  GridTools::partition_triangulation(
    Utilities::MPI::n_mpi_processes(mpi_communicator), triangulation);

  dof_handler.distribute_dofs(fe);
  DoFRenumbering::subdomain_wise(dof_handler);

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

  const std::vector<IndexSet> locally_owned_dofs_per_proc =
    DoFTools::locally_owned_dofs_per_subdomain(dof_handler);

  // Get the DoFs owned by each process
  const IndexSet locally_owned_dofs =
    locally_owned_dofs_per_proc[Utilities::MPI::this_mpi_process(
      mpi_communicator)];

  nodal_counts.resize(n_detector);
  for (Vector<double> &nodal_counts_for_one_detector : nodal_counts)
    nodal_counts_for_one_detector.reinit(locally_owned_dofs.size());

  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);
  solution.reinit(locally_owned_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
}

template <int dim>
void
RPTL2Projection<dim>::solve_linear_system(unsigned detector_no)
{
  TimerOutput::Scope t(computing_timer, "solve_linear_system");

  this->pcout << "Norm of RHS is : " << system_rhs.l2_norm() << std::endl;

  SolverControl   solver_control(10000000, std::max(1e-6,system_rhs.l2_norm()*1e-6);
  LA::SolverGMRES solver(solver_control);

  LA::MPI::PreconditionILU                 preconditioner;
  LA::MPI::PreconditionILU::AdditionalData data(0,1e-10,1,0);
  preconditioner.initialize(system_matrix, data);

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  Vector<double> localized_solution(
    solution); // constructor Trilinos MPI Vector -> deal.II Vector

  constraints.distribute(localized_solution);
  nodal_counts[detector_no] = localized_solution;
}

template <int dim>
void
RPTL2Projection<dim>::assemble_system(unsigned no_detector)
{
  TimerOutput::Scope t(computing_timer, "assemble_system");

  system_rhs    = 0;
  system_matrix = 0;

  const QGaussSimplex<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->subdomain_id() ==
          Utilities::MPI::this_mpi_process(mpi_communicator))
        {
          cell_matrix = 0;
          cell_rhs    = 0;
          fe_values.reinit(cell);

          for (const unsigned int q_index :
               fe_values.quadrature_point_indices())
            {
              Point<dim> q_point_position = fe_values.quadrature_point(q_index);
              RadioParticle<dim> particle(q_point_position, 0);

              ParticleDetectorInteractions<dim> p_q_interaction(
                particle, detectors[no_detector], parameters);

              double count = p_q_interaction.calculate_count();

              for (const unsigned int i : fe_values.dof_indices())
                {
                  for (const unsigned int j : fe_values.dof_indices())
                    {
                      cell_matrix(i, j) +=
                        (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                         fe_values.shape_value(j, q_index) * // phi_j(x_q)
                         fe_values.JxW(q_index));            // dx
                    }
                  cell_rhs(i) +=
                    (count *                             // f(x)
                     fe_values.shape_value(i, q_index) * // phi_i(x_q)
                     fe_values.JxW(q_index));            // dx
                }
            }
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 cell_rhs,
                                                 local_dof_indices,
                                                 system_matrix,
                                                 system_rhs);
        }
    }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
RPTL2Projection<dim>::output_results()
{
  TimerOutput::Scope t(computing_timer, "output_results_vtu");

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::ofstream output("solution.vtk");

      DataOut<dim> data_out;
      data_out.attach_dof_handler(dof_handler);

      for (unsigned d = 0; d < n_detector; ++d)
        {
          data_out.add_data_vector(nodal_counts[d],
                                   "detector_" + Utilities::to_string(d, 2));
        }
      data_out.build_patches();
      data_out.write_vtk(output);
    }
}

template <int dim>
void
RPTL2Projection<dim>::output_raw_results()
{
  TimerOutput::Scope t(computing_timer, "output_results_raw");

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
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
      if (parameters.verbosity == Parameters::Verbosity::verbose)
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
}

template <int dim>
void
RPTL2Projection<dim>::checkpoint()
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
RPTL2Projection<dim>::L2_project()
{
  pcout << "***********************************************" << std::endl;
  pcout << "Assigning detector positions" << std::endl;
  {
    TimerOutput::Scope t(computing_timer, "assigning_detector_positions");
    detectors = assign_detector_positions<dim>(detector_parameters);
  }
  n_detector = detectors.size();
  pcout << "Number of detectors identified: " << n_detector << std::endl;
  pcout << "***********************************************" << std::endl;
  pcout << "Setting up the grid" << std::endl;
  setup_triangulation();
  pcout << "Number of active cells: " << triangulation.n_active_cells()
        << std::endl;
  pcout << "***********************************************" << std::endl;

  setup_system();

  for (unsigned d = 0; d < n_detector; ++d)
    {
      pcout << "Detector_id: " << Utilities::to_string(d, 2) << std::endl;
      pcout << "Assembling system" << std::endl;
      assemble_system(d);
      pcout << "Solving system" << std::endl;
      solve_linear_system(d);
      pcout << "System solved" << std::endl;
      pcout << "-----------------------------------------------" << std::endl;
    }
  pcout << "Outputting results" << std::endl;
  output_results();
  output_raw_results();
  pcout << "-----------------------------------------------" << std::endl;
  pcout << "Saving dof handler and nodal counts" << std::endl;
  checkpoint();
  pcout << "***********************************************" << std::endl;
  pcout << "Done!" << std::endl;
  pcout << "***********************************************" << std::endl;

  // Disable the output of time clock
  if (!fem_reconstruction_parameters.verbose_clock_fem_reconstruction)
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
      // Assembling sys_matrix (Jacobian)
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
      // e_1 = \xi - 1; e_2 = \eta - 1; e_3 = \zeta - 1;
      if (reference_location[i] > 1)
        err_coordinates[i] = (reference_location[i] - 1);

      // e_1 = - \xi; e_2 = - \eta; e_3 = - \zeta;
      if (reference_location[i] < 0)
        err_coordinates[i] = std::abs(reference_location[i]);
    }

  // Check the last constraint
  if (last_constraint < 0)
    {
      // e_4 = |(1 - \xi - \eta - \zeta)|
      err_coordinates[3] = std::abs(last_constraint);
    }
  else if (last_constraint > 1)
    {
      // e_4 = (1 - \xi - \eta - \zeta) - 1
      err_coordinates[3] = last_constraint - 1;
    }

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

  if (fem_reconstruction_parameters.fem_cost_function ==
      Parameters::RPTFEMReconstructionParameters::FEMCostFunction::absolute)
    {
      // cost = \sum_{d=1}^{n_detector} [C_{0,d} * (1 - \xi - \eta - \zeta) +
      // C_{1,d} * \xi + C_{2,d} * \eta + C_{3,d} * \zeta - C_{exp,d}]^2
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
  else if (fem_reconstruction_parameters.fem_cost_function ==
           Parameters::RPTFEMReconstructionParameters::FEMCostFunction::
             relative)
    {
      std::vector<double> denom(n_detector);

      // cost = \sum_{d=1}^{n_detector} [C_{0,d} * (1 - \xi - \eta - \zeta) +
      // C_{1,d} * \xi + C_{2,d} * \eta + C_{3,d} * \zeta - C_{exp,d}]^2 / (
      // C_{exp,d}^2)
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
RPTFEMReconstruction<dim>::search_position_in_reference_space(
  std::vector<std::vector<double>> &count_from_all_detectors,
  std::vector<double> &             experimental_count,
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const double &                                        tol_reference_location,
  double &                                              max_cost,
  Vector<double> &                                      reference_location,
  bool &                                                position_found,
  Point<dim> &                                          real_location)
{
  double last_constraint_reference_location;
  double norm_error_coordinates;
  double calculated_cost;


  // Solve linear system to find the location in reference coordinates
  reference_location = assemble_matrix_and_rhs<dim>(
    count_from_all_detectors,
    experimental_count,
    fem_reconstruction_parameters.fem_cost_function);

  // 4th constraint on the location of the particle in reference coordinates
  last_constraint_reference_location =
    1 - reference_location[0] - reference_location[1] - reference_location[2];

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
      if (calculated_cost < max_cost)
        {
          max_cost               = calculated_cost;
          position_found         = true;
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

template <int dim>
bool
RPTFEMReconstruction<dim>::find_position_global_search(
  std::vector<double> &experimental_counts,
  const double         tol_reference_location)
{
  double                           max_cost       = DBL_MAX;
  bool                             position_found = false;
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
      search_position_in_reference_space(count_from_all_detectors,
                                         experimental_counts,
                                         cell,
                                         tol_reference_location,
                                         max_cost,
                                         reference_location,
                                         position_found,
                                         real_location);
    }
  found_positions.push_back(real_location);
  return position_found;
}

template <int dim>
bool
RPTFEMReconstruction<dim>::find_position_local_search(
  std::vector<double> &experimental_counts,
  const double         tol_reference_location)
{
  double                           max_cost       = DBL_MAX;
  bool                             position_found = false;
  Point<dim>                       real_location;
  Vector<double>                   reference_location;
  std::vector<std::vector<double>> count_from_all_detectors(
    n_detector, std::vector<double>(4));
  const typename DoFHandler<dim>::active_cell_iterator cell =
    previous_position_cell;

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
  if (fem_reconstruction_parameters.search_proximity_level > 1)
    {
      std::set<typename DoFHandler<dim>::active_cell_iterator>
        previous_adjacent_cells;
      std::set<typename DoFHandler<dim>::active_cell_iterator>
        previous_main_cells;

      previous_main_cells.insert(cell);

      for (unsigned int proximity_level = 2;
           proximity_level <
           fem_reconstruction_parameters.search_proximity_level + 1;
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

      search_position_in_reference_space(count_from_all_detectors,
                                         experimental_counts,
                                         adjacent_cell,
                                         tol_reference_location,
                                         max_cost,
                                         reference_location,
                                         position_found,
                                         real_location);
    }

  if (position_found)
    {
      found_positions.push_back(real_location);
      return position_found;
    }
  else
    {
      position_found = find_position_global_search(experimental_counts,
                                                   tol_reference_location);
      return position_found;
    }
}

template <int dim>
void
RPTFEMReconstruction<dim>::trajectory()
{
  // Tolerance or the extrapolation limit in the reference space for a found
  // position
  double tol_reference_location =
    fem_reconstruction_parameters.extrapolation_tolerance;

  if (fem_reconstruction_parameters.mesh_type ==
      Parameters::RPTFEMReconstructionParameters::FEMMeshType::dealii)
    {
      if (fem_reconstruction_parameters.extrapolation_tolerance == 0.005)
        {
          const unsigned int power =
            pow(2, fem_reconstruction_parameters.mesh_refinement);
          const unsigned int n_cell_z =
            2 * fem_reconstruction_parameters.z_subdivisions * power;
          tol_reference_location = parameters.reactor_height / n_cell_z * 1.15;
        }
    }

  // Read and store all experimental counts
  std::vector<std::vector<double>> all_experimental_counts;
  {
    TimerOutput::Scope t(computing_timer, "input_experimental_counts_file");
    all_experimental_counts = read_detectors_counts<dim>(
      fem_reconstruction_parameters.experimental_counts_file, n_detector);
  }

  {
    TimerOutput::Scope t(computing_timer, "find_particle_positions");

    if (fem_reconstruction_parameters.search_type ==
        Parameters::RPTFEMReconstructionParameters::FEMSearchType::local)
      {
        // It is set to "false" to force the global search on the first search
        bool adjacent_cell_search = false;

        // Find the position of the particle with the experimental counts by
        // prioritizing local search
        for (std::vector<double> &experimental_counts : all_experimental_counts)
          {
            if (adjacent_cell_search)
              adjacent_cell_search =
                find_position_local_search(experimental_counts,
                                           tol_reference_location);
            else
              adjacent_cell_search =
                find_position_global_search(experimental_counts,
                                            tol_reference_location);
          }
      }
    else
      {
        // Find the position of the particle with the experimental counts by
        // global search only
        for (std::vector<double> &experimental_counts : all_experimental_counts)
          {
            find_position_global_search(experimental_counts,
                                        tol_reference_location);
          }
      }
  }
}



template <int dim>
void
RPTFEMReconstruction<dim>::load_from_checkpoint()
{
  TimerOutput::Scope t(computing_timer, "load_from_checkpoint");

  n_detector = fem_reconstruction_parameters.nodal_counts_file.size();

  // Import dof handler
  {
    dof_handler.distribute_dofs(fe);
    std::ifstream ifs(fem_reconstruction_parameters.dof_handler_file);
    boost::archive::text_iarchive ia(ifs);
    dof_handler.load(ia, 0);
  }
  // Import nodal counts
  {
    Vector<double> counts_per_detector;

    nodal_counts.resize(n_detector);

    for (unsigned int i = 0; i < n_detector; ++i)
      {
        std::ifstream ifs(fem_reconstruction_parameters.nodal_counts_file[i]);
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

  std::string filename = fem_reconstruction_parameters.export_positions_file;

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
  if (parameters.verbosity == Parameters::Verbosity::verbose)
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
                myfile << position[0] << sep << position[1] << sep << position[2] << std::endl;
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
  if (!fem_reconstruction_parameters.verbose_clock_fem_reconstruction)
    computing_timer.disable_output();
}

template class RPTL2Projection<3>;
template class RPTFEMReconstruction<3>;
