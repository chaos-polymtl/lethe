#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/numerics/data_out.h>

#include <rpt/particle_detector_interactions.h>
#include <rpt/rpt_fem_reconstruction.h>

#include <fstream>

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

  const QGaussSimplex<dim> quadrature_formula(fe.degree + 1);
  FEValues<dim>            fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  const unsigned int       dofs_per_cell = fe.n_dofs_per_cell();
  FullMatrix<double>       cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>           cell_rhs(dofs_per_cell);
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
            particle, detectors[detector_no], rpt_parameters.rpt_param);

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

      detectors.push_back(detector);
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
RPTFEMReconstruction<dim>::output_raw_results_per_level()
{
  for (unsigned int level = 0; level < triangulation.n_levels(); ++level)
    {
      std::map<types::global_dof_index, Point<dim>> dof_index_and_location;

      for (const auto &cell : dof_handler.cell_iterators_on_level(level))
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
      output_counts_on_level(level, dof_index_and_location);
    }
}

template <int dim>
void
RPTFEMReconstruction<dim>::output_counts_on_level(
  unsigned int                                   level,
  std::map<types::global_dof_index, Point<dim>> &dof_index_and_location)
{
  std::cout << "Ouputing on level : " << level << std::endl;
  std::string filename =
    "raw_counts_" + Utilities::int_to_string(level) + ".dat";

  // Open a file
  std::ofstream myfile;
  std::string   sep;
  myfile.open(filename);
  myfile << "vertex_positions_x vertex_position_y vertex_position_z ";
  for (unsigned int i = 0; i < detectors.size(); ++i)
    myfile << "detector_" + Utilities::int_to_string(i) + " ";
  myfile << std::endl;
  sep = " ";

  // showing contents:
  for (auto it = dof_index_and_location.begin();
       it != dof_index_and_location.end();
       ++it)
    {
      for (unsigned d = 0; d < dim; ++d)
        myfile << it->second[d] << sep;
      for (unsigned int i = 0; i < detectors.size(); ++i)
        myfile << nodal_counts[i][it->first] << sep;

      myfile << "\n";
    }
  myfile.close();
}

template <int dim>
void
RPTFEMReconstruction<dim>::L2_project()
{
  std::cout << "Assigning detector positions" << std::endl;
  assign_detector_positions();
  std::cout << "Number of detectors identified : " << detectors.size()
            << std::endl;

  // flatten the triangulation
  Triangulation<dim> temp_triangulation;
  Triangulation<dim> flat_temp_triangulation;
  GridGenerator::cylinder(temp_triangulation, 0.1, 0.2);
  temp_triangulation.refine_global(3);

  GridGenerator::flatten_triangulation(temp_triangulation,
                                       flat_temp_triangulation);

  GridGenerator::convert_hypercube_to_simplex_mesh(flat_temp_triangulation,
                                                   triangulation);

  triangulation.set_all_manifold_ids(0);


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
      std::cout << "System solved" << std::endl;
    }

  std::cout << "Outputting results" << std::endl;
  output_results();
  output_raw_results_per_level();
}


template class RPTFEMReconstruction<3>;


//This new function takes two arguments, a matrix of nodal counts from all the detectors
//and a vector of tracer counts from all the detectors
//This function assemble system_matrix and system_rhs
template <int dim>
void
assemble_matrix_and_rhs(std::vector<std::vector<double>> vertex_count, std::vector<double> experimental_count,int detector_size)
{

    //Assembling system_matrix


    std::vector<std::vector<double>>system_matrix;

    double sigma=0;
    for ( int i=0; i<detector_size;i++){
        sigma+=(-vertex_count[i][0]+vertex_count[i][1])*(-vertex_count[i][0]+vertex_count[i][1]);
    }
    system_matrix[0][0]=sigma;


    sigma=0;
    for (int i=0; i<detector_size;i++){
        sigma+=(-vertex_count[i][0]+vertex_count[i][2])*(-vertex_count[i][0]+vertex_count[i][1]);
    }
    system_matrix[0][1]=sigma;



    sigma=0;
    for (int i=0; i<detector_size;i++){
        sigma+=(-vertex_count[i][0]+vertex_count[i][3])*(-vertex_count[i][0]+vertex_count[i][1]);
    }
    system_matrix[0][2]=sigma;


    sigma=0;
    for (int i=0; i<detector_size;i++){
        sigma+=(-vertex_count[i][0]+vertex_count[i][2])*(-vertex_count[i][0]+vertex_count[i][1]);
    }
    system_matrix[1][0]=sigma;

    sigma=0;
    for (int i=0; i<detector_size;i++){
        sigma+=(-vertex_count[i][0]+vertex_count[i][2])*(-vertex_count[i][0]+vertex_count[i][2]);
    }

    system_matrix[1][1]=sigma;


    sigma=0;
    for(int i=0; i<detector_size;i++){
        sigma+=(-vertex_count[i][0]+vertex_count[i][3])*(-vertex_count[i][0]+vertex_count[i][2]);
    }
    system_matrix[1][2]=sigma;

    sigma=0;
    for(int i=0;i<detector_size;i++){
        sigma+=(-vertex_count[i][0]+vertex_count[i][3])*(-vertex_count[i][0]+vertex_count[i][1]);
    }
    system_matrix[2][0]=sigma;

    sigma=0;
    for(int i=0; i<detector_size;i++){
        sigma+=(-vertex_count[i][0]+vertex_count[i][3])*(-vertex_count[i][0]+vertex_count[i][2]);
    }
    system_matrix[2][1]=sigma;

    sigma=0;
    for (int i=0;i<detector_size;i++){
        sigma+=(-vertex_count[i][0]+vertex_count[i][3])*(-vertex_count[i][0]+vertex_count[i][3]);
    }
    system_matrix[2][2]=sigma;


    //Assembling system_rhs
    Vector<double> system_rhs;

    sigma=0;
    for (int i=0; i<detector_size;i++){
        sigma+=(vertex_count[i][0]-experimental_count[i])*(-vertex_count[i][0]+vertex_count[i][1]);
    }
    system_rhs[0]=sigma;

    sigma=0;
    for (int i=0; i<detector_size;i++){
        sigma+=(vertex_count[i][0]-experimental_count[i])*(-vertex_count[i][0]+vertex_count[i][2]);
    }
    system_rhs[1]=sigma;

    sigma=0;
    for (int i=0; i<detector_size;i++){
        sigma+=(vertex_count[i][0]-experimental_count[i])*(-vertex_count[i][0]+vertex_count[i][3]);
    }
    system_rhs[2]=sigma;
    //SOLVE

    Vector<double> solution;
    SolverControl solver_control(1000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

}


