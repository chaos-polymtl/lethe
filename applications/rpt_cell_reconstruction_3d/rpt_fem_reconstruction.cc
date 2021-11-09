#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/numerics/data_out.h>

#include <rpt/particle_detector_interactions.h>
#include <rpt/rpt_fem_reconstruction.h>



using namespace dealii;



template <int dim>
void
RPTFEMReconstruction<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  system_rhs.reinit(dof_handler.n_dofs());
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);

  nodal_counts.resize(detectors.size());
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
RPTFEMReconstruction<dim>::assemble_system(unsigned detector_no)
{
  system_rhs    = 0;
  system_matrix = 0;

  const QGauss<dim>  quadrature_formula(fe.degree + 1);
  FEValues<dim>      fe_values(fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  unsigned int                         cell_calculated = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // std::cout << " Cell id : " << ++cell_calculated << std::endl;
      cell_matrix = 0;
      cell_rhs    = 0;
      fe_values.reinit(cell);
      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          Point<dim> q_point_position = fe_values.quadrature_point(q_index);
          RadioParticle<dim> particle(q_point_position, 0);

          ParticleDetectorInteractions<dim> p_q_interaction(
            particle, detectors[detector_no], rpt_parameters);

          double count = p_q_interaction.calculate_count();

          for (const unsigned int i : fe_values.dof_indices())
            {
              for (const unsigned int j : fe_values.dof_indices())
                cell_matrix(i, j) +=
                  (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                   fe_values.shape_value(j, q_index) * // phi_j(x_q)
                   fe_values.JxW(q_index));            // dx

              cell_rhs(i) += (count *                             // f(x)
                              fe_values.shape_value(i, q_index) * // phi_i(x_q)
                              fe_values.JxW(q_index));            // dx
            }
        }
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
}


template <int dim>
void
RPTFEMReconstruction<dim>::assign_detector_positions()
{
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

      {
        detectors.push_back(detector);
      }
    }
}

template <int dim>
void
RPTFEMReconstruction<dim>::output_results()
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  for (unsigned d = 0; d < detectors.size(); ++d)
    {
      data_out.add_data_vector(nodal_counts[d],
                               "detector_" + std::to_string(d));
    }
  data_out.build_patches();
  std::ofstream output("solution.vtu");
  data_out.write_vtu(output);
}

template <int dim>
void
RPTFEMReconstruction<dim>::L2_project()
{
  std::cout << "Assigning detector positions" << std::endl;
  assign_detector_positions();
  std::cout << "Number of detectors identified : " << detectors.size()
            << std::endl;

  GridGenerator::cylinder(triangulation, 0.1, 0.2);
  triangulation.refine_global(3);


  GridOut grid_out;
  {
    std::ofstream output_file("original_triangulation.vtk");
    grid_out.write_vtk(triangulation, output_file);
  }

  GridTools::rotate(1.57078, 1, triangulation);
  Tensor<1, dim> shift_vector({0, 0, 0.2});
  GridTools::shift(shift_vector, triangulation);
  {
    std::ofstream output_file("rotated_triangulation.vtk");
    grid_out.write_vtk(triangulation, output_file);
  }


  setup_system();
  for (unsigned d = 0; d < detectors.size(); ++d)
    {
      std::cout << "Assembling system" << std::endl;
      assemble_system(d);
      std::cout << "Solving system" << std::endl;
      solve_linear_system(d);
    }
  // should I add a function that instead of out putting the result, save it in
  // a form of matrix?
  output_results();
  // values();
  // assemble_jacobian();
  // assemble_rhs(initial_guess,cell);
  // Calculate_Jacobian_1(initial_guess,cell);
  solve();
  // assemble_jacobian_for_Newton_method();
}

template <int dim>
void
RPTFEMReconstruction<dim>::values()
{
  Vector<double> c1(
    {212.332, 236.828, 211.65, 236.008, 211.546, 236.099, 210.901, 235.401});
  Vector<double> c2(
    {71.1969, 76.7869, 66.4289, 71.641, 72.9663, 78.6936, 67.9966, 73.4223});
  Vector<double> c3(
    {22.852, 22.273, 20.858, 20.3253, 23.6624, 23.0588, 21.5689, 21.05});
  for (int i = 0; i < 8; i++)
    {
      c[0][i] = c1[i];
    }
  for (int i = 0; i < 8; i++)
    {
      c[1][i] = c2[i];
    }

  for (int i = 0; i < 8; i++)
    {
      c[2][i] = c3[i];
    }
}


template <int dim>
void
RPTFEMReconstruction<dim>::assemble_jacobian_for_Newton_method()
{
  jacobian_matrix[0][0] = Calculate_Jacobian_1();
  jacobian_matrix[0][1] = Calculate_Jacobian_2();
  jacobian_matrix[0][2] = Calculate_Jacobian_3();
  jacobian_matrix[1][0] = Calculate_Jacobian_4();
  jacobian_matrix[1][1] = Calculate_Jacobian_5();
  jacobian_matrix[1][2] = Calculate_Jacobian_6();
  jacobian_matrix[2][0] = Calculate_Jacobian_7();
  jacobian_matrix[2][1] = Calculate_Jacobian_8();
  jacobian_matrix[2][2] = Calculate_Jacobian_9();
}


template <int dim>
double
RPTFEMReconstruction<dim>::Calculate_Jacobian_1()
{
  double sigma = 0;



  for (int i = 0; i < 3; i++)
    {
      sigma += (-0.125 * c[i][0] * (1 - unknown[2]) * (1 - unknown[1]) +
                0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) -
                0.125 * c[i][2] * (1 - unknown[2]) * (unknown[1] + 1) +
                0.125 * c[i][3] * (1 - unknown[2]) * (unknown[1] + 1) -
                0.125 * c[i][4] * (1 - unknown[1]) * (unknown[2] + 1) +
                0.125 * c[i][5] * (1 - unknown[1]) * (unknown[2] + 1) -
                0.125 * c[i][6] * (unknown[2] + 1) * (unknown[1] + 1) +
                0.125 * c[i][7] * (unknown[2] + 1) * (unknown[1] + 1)) *
               (c[i][0] * (1 - unknown[1]) * (0.125 * unknown[2] - 0.125) +
                c[i][1] * (0.125 - 0.125 * unknown[2]) * (1 - unknown[1]) +
                c[i][2] * (1 - unknown[2]) * (-0.125 * unknown[1] - 0.125) +
                c[i][3] * (1 - unknown[2]) * (0.125 * unknown[1] + 0.125) +
                c[i][4] * (1 - unknown[1]) * (-0.125 * unknown[2] - 0.125) +
                c[i][5] * (1 - unknown[1]) * (0.125 * unknown[2] + 0.125) -
                c[i][6] * (0.125 * unknown[2] + 0.125) * (unknown[1] + 1) +
                c[i][7] * (0.125 * unknown[2] + 0.125) * (unknown[1] + 1));
    }
  return sigma;
}

template <int dim>
double
RPTFEMReconstruction<dim>::Calculate_Jacobian_2()
{
  double sigma = 0;



  for (unsigned int i = 0; i < 3; i++)
    {
      sigma += (-c[i][0] * (0.125 * unknown[2] - 0.125) -
                c[i][1] * (0.125 - 0.125 * unknown[2]) -
                0.125 * c[i][2] * (1 - unknown[2]) +
                0.125 * c[i][3] * (1 - unknown[2]) -
                c[i][4] * (-0.125 * unknown[2] - 0.125) -
                c[i][5] * (0.125 * unknown[2] + 0.125) -
                c[i][6] * (0.125 * unknown[2] + 0.125) +
                c[i][7] * (0.125 * unknown[2] + 0.125)) *
                 (0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (1 - unknown[1]) +
                  0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) *
                    (unknown[0] + 1) +
                  0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (unknown[1] + 1) +
                  0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) *
                    (unknown[2] + 1) +
                  0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) *
                    (unknown[2] + 1) +
                  0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1) *
                    (unknown[1] + 1) -
                  experimental_count[i]) +
               (-0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) -
                0.125 * c[i][1] * (1 - unknown[2]) * (unknown[0] + 1) +
                0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) +
                0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) -
                0.125 * c[i][4] * (1 - unknown[0]) * (unknown[2] + 1) -
                0.125 * c[i][5] * (unknown[0] + 1) * (unknown[2] + 1) +
                0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) +
                0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1)) *
                 (c[i][0] * (1 - unknown[1]) * (0.125 * unknown[2] - 0.125) +
                  c[i][1] * (0.125 - 0.125 * unknown[2]) * (1 - unknown[1]) +
                  c[i][2] * (1 - unknown[2]) * (-0.125 * unknown[1] - 0.125) +
                  c[i][3] * (1 - unknown[2]) * (0.125 * unknown[1] + 0.125) +
                  c[i][4] * (1 - unknown[1]) * (-0.125 * unknown[2] - 0.125) +
                  c[i][5] * (1 - unknown[1]) * (0.125 * unknown[2] + 0.125) -
                  c[i][6] * (0.125 * unknown[2] + 0.125) * (unknown[1] + 1) +
                  c[i][7] * (0.125 * unknown[2] + 0.125) * (unknown[1] + 1));
    }
  std::cout << sigma << std::endl;
  return sigma;
}

template <int dim>
double
RPTFEMReconstruction<dim>::Calculate_Jacobian_3()
{
  double sigma = 0;



  for (unsigned int i = 0; i < 3; i++)
    {
      sigma += (0.125 * c[i][0] * (1 - unknown[1]) -
                0.125 * c[i][1] * (1 - unknown[1]) -
                c[i][2] * (-0.125 * unknown[1] - 0.125) -
                c[i][3] * (0.125 * unknown[1] + 0.125) -
                0.125 * c[i][4] * (1 - unknown[1]) +
                0.125 * c[i][5] * (1 - unknown[1]) -
                0.125 * c[i][6] * (unknown[1] + 1) +
                0.125 * c[i][7] * (unknown[1] + 1)) *
                 (0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (1 - unknown[1]) +
                  0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) *
                    (unknown[0] + 1) +
                  0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (unknown[1] + 1) +
                  0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) *
                    (unknown[2] + 1) +
                  0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) *
                    (unknown[2] + 1) +
                  0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1) *
                    (unknown[1] + 1) -
                  experimental_count[i]) +
               (-0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[1]) -
                0.125 * c[i][1] * (1 - unknown[1]) * (unknown[0] + 1) -
                0.125 * c[i][2] * (1 - unknown[0]) * (unknown[1] + 1) -
                0.125 * c[i][3] * (unknown[0] + 1) * (unknown[1] + 1) +
                0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) +
                0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) +
                0.125 * c[i][6] * (1 - unknown[0]) * (unknown[1] + 1) +
                0.125 * c[i][7] * (unknown[0] + 1) * (unknown[1] + 1)) *
                 (c[i][0] * (1 - unknown[1]) * (0.125 * unknown[2] - 0.125) +
                  c[i][1] * (0.125 - 0.125 * unknown[2]) * (1 - unknown[1]) +
                  c[i][2] * (1 - unknown[2]) * (-0.125 * unknown[1] - 0.125) +
                  c[i][3] * (1 - unknown[2]) * (0.125 * unknown[1] + 0.125) +
                  c[i][4] * (1 - unknown[1]) * (-0.125 * unknown[2] - 0.125) +
                  c[i][5] * (1 - unknown[1]) * (0.125 * unknown[2] + 0.125) -
                  c[i][6] * (0.125 * unknown[2] + 0.125) * (unknown[1] + 1) +
                  c[i][7] * (0.125 * unknown[2] + 0.125) * (unknown[1] + 1));
    }
  std::cout << sigma << std::endl;
  return sigma;
}

template <int dim>
double
RPTFEMReconstruction<dim>::Calculate_Jacobian_4()
{
  double sigma = 0;



  for (unsigned int i = 0; i < 3; i++)
    {
      sigma += (0.125 * c[i][0] * (1 - unknown[2]) -
                0.125 * c[i][1] * (1 - unknown[2]) -
                0.125 * c[i][2] * (1 - unknown[2]) +
                0.125 * c[i][3] * (1 - unknown[2]) +
                c[i][4] * (0.125 * unknown[2] + 0.125) -
                0.125 * c[i][5] * (unknown[2] + 1) -
                c[i][6] * (0.125 * unknown[2] + 0.125) +
                0.125 * c[i][7] * (unknown[2] + 1)) *
                 (0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (1 - unknown[1]) +
                  0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) *
                    (unknown[0] + 1) +
                  0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (unknown[1] + 1) +
                  0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) *
                    (unknown[2] + 1) +
                  0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) *
                    (unknown[2] + 1) +
                  0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1) *
                    (unknown[1] + 1) -
                  experimental_count[i]) +
               (-c[i][0] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[2]) -
                c[i][1] * (1 - unknown[2]) * (0.125 * unknown[0] + 0.125) +
                c[i][2] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[2]) +
                c[i][3] * (1 - unknown[2]) * (0.125 * unknown[0] + 0.125) -
                c[i][4] * (1 - unknown[0]) * (0.125 * unknown[2] + 0.125) -
                c[i][5] * (0.125 * unknown[0] + 0.125) * (unknown[2] + 1) +
                c[i][6] * (1 - unknown[0]) * (0.125 * unknown[2] + 0.125) +
                c[i][7] * (0.125 * unknown[0] + 0.125) * (unknown[2] + 1)) *
                 (-0.125 * c[i][0] * (1 - unknown[2]) * (1 - unknown[1]) +
                  0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) -
                  0.125 * c[i][2] * (1 - unknown[2]) * (unknown[1] + 1) +
                  0.125 * c[i][3] * (1 - unknown[2]) * (unknown[1] + 1) -
                  0.125 * c[i][4] * (1 - unknown[1]) * (unknown[2] + 1) +
                  0.125 * c[i][5] * (1 - unknown[1]) * (unknown[2] + 1) -
                  0.125 * c[i][6] * (unknown[2] + 1) * (unknown[1] + 1) +
                  0.125 * c[i][7] * (unknown[2] + 1) * (unknown[1] + 1));
    }
}

template <int dim>
double
RPTFEMReconstruction<dim>::Calculate_Jacobian_5()
{
  double sigma = 0;



  for (unsigned int i = 0; i < 3; i++)
    {
      sigma += (-c[i][0] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[2]) -
                c[i][1] * (1 - unknown[2]) * (0.125 * unknown[0] + 0.125) +
                c[i][2] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[2]) +
                c[i][3] * (1 - unknown[2]) * (0.125 * unknown[0] + 0.125) -
                c[i][4] * (1 - unknown[0]) * (0.125 * unknown[2] + 0.125) -
                c[i][5] * (0.125 * unknown[0] + 0.125) * (unknown[2] + 1) +
                c[i][6] * (1 - unknown[0]) * (0.125 * unknown[2] + 0.125) +
                c[i][7] * (0.125 * unknown[0] + 0.125) * (unknown[2] + 1)) *
               (-0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) -
                0.125 * c[i][1] * (1 - unknown[2]) * (unknown[0] + 1) +
                0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) +
                0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) -
                0.125 * c[i][4] * (1 - unknown[0]) * (unknown[2] + 1) -
                0.125 * c[i][5] * (unknown[0] + 1) * (unknown[2] + 1) +
                0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) +
                0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1));
    }
  std::cout << sigma << std::endl;
  return sigma;
}


template <int dim>
double
RPTFEMReconstruction<dim>::Calculate_Jacobian_6()
{
  double sigma = 0;



  for (unsigned int i = 0; i < 3; i++)
    {
      sigma += (c[i][0] * (0.125 - 0.125 * unknown[0]) +
                c[i][1] * (0.125 * unknown[0] + 0.125) -
                c[i][2] * (0.125 - 0.125 * unknown[0]) -
                c[i][3] * (0.125 * unknown[0] + 0.125) -
                0.125 * c[i][4] * (1 - unknown[0]) -
                c[i][5] * (0.125 * unknown[0] + 0.125) +
                0.125 * c[i][6] * (1 - unknown[0]) +
                c[i][7] * (0.125 * unknown[0] + 0.125)) *
                 (0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (1 - unknown[1]) +
                  0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) *
                    (unknown[0] + 1) +
                  0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (unknown[1] + 1) +
                  0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) *
                    (unknown[2] + 1) +
                  0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) *
                    (unknown[2] + 1) +
                  0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1) *
                    (unknown[1] + 1) -
                  experimental_count[i]) +
               (-c[i][0] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[2]) -
                c[i][1] * (1 - unknown[2]) * (0.125 * unknown[0] + 0.125) +
                c[i][2] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[2]) +
                c[i][3] * (1 - unknown[2]) * (0.125 * unknown[0] + 0.125) -
                c[i][4] * (1 - unknown[0]) * (0.125 * unknown[2] + 0.125) -
                c[i][5] * (0.125 * unknown[0] + 0.125) * (unknown[2] + 1) +
                c[i][6] * (1 - unknown[0]) * (0.125 * unknown[2] + 0.125) +
                c[i][7] * (0.125 * unknown[0] + 0.125) * (unknown[2] + 1)) *
                 (-0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[1]) -
                  0.125 * c[i][1] * (1 - unknown[1]) * (unknown[0] + 1) -
                  0.125 * c[i][2] * (1 - unknown[0]) * (unknown[1] + 1) -
                  0.125 * c[i][3] * (unknown[0] + 1) * (unknown[1] + 1) +
                  0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) +
                  0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) +
                  0.125 * c[i][6] * (1 - unknown[0]) * (unknown[1] + 1) +
                  0.125 * c[i][7] * (unknown[0] + 1) * (unknown[1] + 1));
    }
  std::cout << sigma << std::endl;
  return sigma;
}

template <int dim>
double
RPTFEMReconstruction<dim>::Calculate_Jacobian_7()
{
  double sigma = 0;



  for (unsigned int i = 0; i < 3; i++)
    {
      sigma += (0.125 * c[i][0] * (1 - unknown[1]) -
                0.125 * c[i][1] * (1 - unknown[1]) +
                c[i][2] * (0.125 * unknown[1] + 0.125) -
                0.125 * c[i][3] * (unknown[1] + 1) -
                0.125 * c[i][4] * (1 - unknown[1]) +
                0.125 * c[i][5] * (1 - unknown[1]) -
                c[i][6] * (0.125 * unknown[1] + 0.125) +
                0.125 * c[i][7] * (unknown[1] + 1)) *
                 (0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (1 - unknown[1]) +
                  0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) *
                    (unknown[0] + 1) +
                  0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (unknown[1] + 1) +
                  0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) *
                    (unknown[2] + 1) +
                  0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) *
                    (unknown[2] + 1) +
                  0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1) *
                    (unknown[1] + 1) -
                  experimental_count[i]) +
               (-0.125 * c[i][0] * (1 - unknown[2]) * (1 - unknown[1]) +
                0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) -
                0.125 * c[i][2] * (1 - unknown[2]) * (unknown[1] + 1) +
                0.125 * c[i][3] * (1 - unknown[2]) * (unknown[1] + 1) -
                0.125 * c[i][4] * (1 - unknown[1]) * (unknown[2] + 1) +
                0.125 * c[i][5] * (1 - unknown[1]) * (unknown[2] + 1) -
                0.125 * c[i][6] * (unknown[2] + 1) * (unknown[1] + 1) +
                0.125 * c[i][7] * (unknown[2] + 1) * (unknown[1] + 1)) *
                 (c[i][0] * (1 - unknown[1]) * (0.125 * unknown[0] - 0.125) +
                  c[i][1] * (1 - unknown[1]) * (-0.125 * unknown[0] - 0.125) -
                  c[i][2] * (1 - unknown[0]) * (0.125 * unknown[1] + 0.125) -
                  c[i][3] * (0.125 * unknown[0] + 0.125) * (unknown[1] + 1) +
                  c[i][4] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[1]) +
                  c[i][5] * (1 - unknown[1]) * (0.125 * unknown[0] + 0.125) +
                  c[i][6] * (1 - unknown[0]) * (0.125 * unknown[1] + 0.125) +
                  c[i][7] * (0.125 * unknown[0] + 0.125) * (unknown[1] + 1));
    }
  std::cout << sigma << std::endl;
  return sigma;
}

template <int dim>
double
RPTFEMReconstruction<dim>::Calculate_Jacobian_8()
{
  double sigma = 0;



  for (unsigned int i = 0; i < 3; i++)
    {
      sigma += (-c[i][0] * (0.125 * unknown[0] - 0.125) -
                c[i][1] * (-0.125 * unknown[0] - 0.125) -
                0.125 * c[i][2] * (1 - unknown[0]) -
                c[i][3] * (0.125 * unknown[0] + 0.125) -
                c[i][4] * (0.125 - 0.125 * unknown[0]) -
                c[i][5] * (0.125 * unknown[0] + 0.125) +
                0.125 * c[i][6] * (1 - unknown[0]) +
                c[i][7] * (0.125 * unknown[0] + 0.125)) *
                 (0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (1 - unknown[1]) +
                  0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) *
                    (unknown[0] + 1) +
                  0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) *
                    (unknown[1] + 1) +
                  0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) *
                    (unknown[2] + 1) +
                  0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) *
                    (unknown[2] + 1) +
                  0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) *
                    (unknown[1] + 1) +
                  0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1) *
                    (unknown[1] + 1) -
                  experimental_count[i]) +
               (-0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) -
                0.125 * c[i][1] * (1 - unknown[2]) * (unknown[0] + 1) +
                0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) +
                0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) -
                0.125 * c[i][4] * (1 - unknown[0]) * (unknown[2] + 1) -
                0.125 * c[i][5] * (unknown[0] + 1) * (unknown[2] + 1) +
                0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) +
                0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1)) *
                 (c[i][0] * (1 - unknown[1]) * (0.125 * unknown[0] - 0.125) +
                  c[i][1] * (1 - unknown[1]) * (-0.125 * unknown[0] - 0.125) -
                  c[i][2] * (1 - unknown[0]) * (0.125 * unknown[1] + 0.125) -
                  c[i][3] * (0.125 * unknown[0] + 0.125) * (unknown[1] + 1) +
                  c[i][4] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[1]) +
                  c[i][5] * (1 - unknown[1]) * (0.125 * unknown[0] + 0.125) +
                  c[i][6] * (1 - unknown[0]) * (0.125 * unknown[1] + 0.125) +
                  c[i][7] * (0.125 * unknown[0] + 0.125) * (unknown[1] + 1));
    }
  std::cout << sigma << std::endl;
  return sigma;
}

template <int dim>
double
RPTFEMReconstruction<dim>::Calculate_Jacobian_9()
{
  double sigma = 0;



  for (unsigned int i = 0; i < 3; i++)
    {
      sigma += (-0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[1]) -
                0.125 * c[i][1] * (1 - unknown[1]) * (unknown[0] + 1) -
                0.125 * c[i][2] * (1 - unknown[0]) * (unknown[1] + 1) -
                0.125 * c[i][3] * (unknown[0] + 1) * (unknown[1] + 1) +
                0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) +
                0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) +
                0.125 * c[i][6] * (1 - unknown[0]) * (unknown[1] + 1) +
                0.125 * c[i][7] * (unknown[0] + 1) * (unknown[1] + 1)) *
               (c[i][0] * (1 - unknown[1]) * (0.125 * unknown[0] - 0.125) +
                c[i][1] * (1 - unknown[1]) * (-0.125 * unknown[0] - 0.125) -
                c[i][2] * (1 - unknown[0]) * (0.125 * unknown[1] + 0.125) -
                c[i][3] * (0.125 * unknown[0] + 0.125) * (unknown[1] + 1) +
                c[i][4] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[1]) +
                c[i][5] * (1 - unknown[1]) * (0.125 * unknown[0] + 0.125) +
                c[i][6] * (1 - unknown[0]) * (0.125 * unknown[1] + 0.125) +
                c[i][7] * (0.125 * unknown[0] + 0.125) * (unknown[1] + 1));
    }
  std::cout << sigma << std::endl;
  return sigma;
}

template <int dim>
double
RPTFEMReconstruction<dim>::f1()
{
  double sigma = 0;



  for (unsigned int i = 0; i < 3; i++)
    {
      sigma += (c[i][0] * (1 - unknown[1]) * (0.125 * unknown[2] - 0.125) +
                c[i][1] * (0.125 - 0.125 * unknown[2]) * (1 - unknown[1]) +
                c[i][2] * (1 - unknown[2]) * (-0.125 * unknown[1] - 0.125) +
                c[i][3] * (1 - unknown[2]) * (0.125 * unknown[1] + 0.125) +
                c[i][4] * (1 - unknown[1]) * (-0.125 * unknown[2] - 0.125) +
                c[i][5] * (1 - unknown[1]) * (0.125 * unknown[2] + 0.125) -
                c[i][6] * (0.125 * unknown[2] + 0.125) * (unknown[1] + 1) +
                c[i][7] * (0.125 * unknown[2] + 0.125) * (unknown[1] + 1)) *
               (0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) *
                  (1 - unknown[1]) +
                0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) *
                  (unknown[0] + 1) +
                0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) *
                  (unknown[1] + 1) +
                0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) *
                  (unknown[1] + 1) +
                0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) *
                  (unknown[2] + 1) +
                0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) *
                  (unknown[2] + 1) +
                0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) *
                  (unknown[1] + 1) +
                0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1) *
                  (unknown[1] + 1) -
                experimental_count[i]);
    }
  std::cout << sigma << std::endl;
  return sigma;
}

template <int dim>
double
RPTFEMReconstruction<dim>::f2()
{
  double sigma = 0;



  for (unsigned int i = 0; i < 3; i++)
    {
      sigma += (-c[i][0] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[2]) -
                c[i][1] * (1 - unknown[2]) * (0.125 * unknown[0] + 0.125) +
                c[i][2] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[2]) +
                c[i][3] * (1 - unknown[2]) * (0.125 * unknown[0] + 0.125) -
                c[i][4] * (1 - unknown[0]) * (0.125 * unknown[2] + 0.125) -
                c[i][5] * (0.125 * unknown[0] + 0.125) * (unknown[2] + 1) +
                c[i][6] * (1 - unknown[0]) * (0.125 * unknown[2] + 0.125) +
                c[i][7] * (0.125 * unknown[0] + 0.125) * (unknown[2] + 1)) *
               (0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) *
                  (1 - unknown[1]) +
                0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) *
                  (unknown[0] + 1) +
                0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) *
                  (unknown[1] + 1) +
                0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) *
                  (unknown[1] + 1) +
                0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) *
                  (unknown[2] + 1) +
                0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) *
                  (unknown[2] + 1) +
                0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) *
                  (unknown[1] + 1) +
                0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1) *
                  (unknown[1] + 1) -
                experimental_count[i]);
      ;
    }
  std::cout << sigma << std::endl;
  return sigma;
}

template <int dim>
double
RPTFEMReconstruction<dim>::f3()
{
  double sigma = 0;



  for (unsigned int i = 0; i < 3; i++)
    {
      sigma += (c[i][0] * (1 - unknown[1]) * (0.125 * unknown[0] - 0.125) +
                c[i][1] * (1 - unknown[1]) * (-0.125 * unknown[0] - 0.125) -
                c[i][2] * (1 - unknown[0]) * (0.125 * unknown[1] + 0.125) -
                c[i][3] * (0.125 * unknown[0] + 0.125) * (unknown[1] + 1) +
                c[i][4] * (0.125 - 0.125 * unknown[0]) * (1 - unknown[1]) +
                c[i][5] * (1 - unknown[1]) * (0.125 * unknown[0] + 0.125) +
                c[i][6] * (1 - unknown[0]) * (0.125 * unknown[1] + 0.125) +
                c[i][7] * (0.125 * unknown[0] + 0.125) * (unknown[1] + 1)) *
               (0.125 * c[i][0] * (1 - unknown[0]) * (1 - unknown[2]) *
                  (1 - unknown[1]) +
                0.125 * c[i][1] * (1 - unknown[2]) * (1 - unknown[1]) *
                  (unknown[0] + 1) +
                0.125 * c[i][2] * (1 - unknown[0]) * (1 - unknown[2]) *
                  (unknown[1] + 1) +
                0.125 * c[i][3] * (1 - unknown[2]) * (unknown[0] + 1) *
                  (unknown[1] + 1) +
                0.125 * c[i][4] * (1 - unknown[0]) * (1 - unknown[1]) *
                  (unknown[2] + 1) +
                0.125 * c[i][5] * (1 - unknown[1]) * (unknown[0] + 1) *
                  (unknown[2] + 1) +
                0.125 * c[i][6] * (1 - unknown[0]) * (unknown[2] + 1) *
                  (unknown[1] + 1) +
                0.125 * c[i][7] * (unknown[0] + 1) * (unknown[2] + 1) *
                  (unknown[1] + 1) -
                experimental_count[i]);
    }
  std::cout << sigma << std::endl;
  return sigma;
}

template <int dim>
void
RPTFEMReconstruction<dim>::assemble_rhs()

{
  rhs_matrix[0] = f1();
  rhs_matrix[1] = f2();
  rhs_matrix[2] = f3();
}


template <int dim>
void
RPTFEMReconstruction<dim>::solve()
{
  experimental_count[0] = 247.096;
  experimental_count[1] = 79.2302;
  experimental_count[2] = 22.0657;
  unknown[0]            = 0.1;
  unknown[1]            = 0.1;
  unknown[2]            = 0.1;
  values();
  assemble_jacobian_for_Newton_method();
  assemble_rhs();

  double tol     = 1e-4;
  double norm_dx = 1;
  while (norm_dx > tol)
    {
      assemble_jacobian_for_Newton_method();
      assemble_rhs();
      Tensor<2, 3, double> inverse_jacobian = invert(jacobian_matrix);
      Tensor<1, 3, double> dx               = -inverse_jacobian * rhs_matrix;
      norm_dx                               = dx.norm();
      unknown += dx;
    }
  std::cout << unknown << std::endl;
}



// From here is the main program
/*
template <int dim>
void
RPTFEMReconstruction<dim>::Loop_over_cells()
{



    unsigned int level=0;
    for (const auto &cell:
         dof_handler.cell_iterators_on_level(level))
    {//active_cell_iterator or cell_iterator?



        //assemble_jacobian_for_Newton_method();
        //assemble_rhs(initial_guess,cell);
        //Calculate_Jacobian_1(initial_guess,cell);
        //assemble_rhs(initial_guess,cell);

    }


}
template<int dim>
void
RPTFEMReconstruction<dim>::assemble_jacobian_for_Newton_method(Tensor<1,dim>
unknown)
{
    jacobian_matrix[0][0]=2*unknown[0];
    jacobian_matrix[0][1]=2*unknown[1];
    jacobian_matrix[0][2]=2*unknown[2];
    jacobian_matrix[1][0]=2*unknown[0];
    jacobian_matrix[1][1]=-2*unknown[1];
    jacobian_matrix[1][2]=4*unknown[2];
    jacobian_matrix[2][0]=4*unknown[0];
    jacobian_matrix[2][1]=2*unknown[1];
    jacobian_matrix[2][2]=-2*unknown[2];


}


template<int dim>
void
RPTFEMReconstruction<dim>::assemble_rhs(Tensor<1,dim> unknown)
//Vector<double>unknown, const typename DoFHandler<dim>::active_cell_iterator
&cell
{


    rhs_matrix[0]=pow(unknown[0],2)+pow(unknown[1],2)+pow(unknown[2],2)-6;
    rhs_matrix[1]=pow(unknown[0],2)-pow(unknown[1],2)+2*pow(unknown[2],2)-2;
    rhs_matrix[2]=2*pow(unknown[0],2)+pow(unknown[1],2)-pow(unknown[2],2)-3;









}

template<int dim>
void
RPTFEMReconstruction<dim>::Calculate_Jacobian_1(Tensor<1,dim> unknown,
                                                const typename
DoFHandler<dim>::active_cell_iterator &cell){


    double sigma=0;
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);

    for (unsigned i=0;i<detectors.size();i++){


        sigma+=

}


}

RPTFEMReconstruction<dim>::Calculate_Jacobian_2(Tensor<1,dim> unknown,
                                                const typename
DoFHandler<dim>::active_cell_iterator &cell){


    double sigma=0;
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);

    for (unsigned i=0;i<detectors.size();i++){


        sigma+=


}


}

template<int dim>
void
RPTFEMReconstruction<dim>::solve()
{
    initial_guess[0]=1;
    initial_guess[1]=1;
    initial_guess[2]=1;
    double tol=1e-4;
    double norm_dx=1;
    std::cout<<"ok1"<<std::endl;
    while(norm_dx>tol){
        assemble_jacobian_for_Newton_method(initial_guess);
        std::cout<<jacobian_matrix<<std::endl;
        assemble_rhs(initial_guess);
        std::cout<<rhs_matrix<<std::endl;
        Tensor<2,3,double> inverse_jacobian= invert(jacobian_matrix);
        std::cout<<inverse_jacobian<<std::endl;
        Tensor<1,3,double> dx= -inverse_jacobian*rhs_matrix;
        std::cout<<dx<<std::endl;
        norm_dx=dx.norm();
        //std::cout<<norm_dx<<std::endl;
        initial_guess+=dx;
         std::cout<<dx<<std::endl;
        std::cout<<norm_dx<<std::endl;


    }
    std::cout<<initial_guess<<std::endl;



}
*/

template class RPTFEMReconstruction<3>;
