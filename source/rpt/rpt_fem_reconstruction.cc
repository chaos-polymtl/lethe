#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/numerics/data_out.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <rpt/particle_detector_interactions.h>
#include <rpt/rpt_fem_reconstruction.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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
    //for (unsigned d = 0; d < detectors.size(); ++d)
    //ajouté
    for (unsigned d = 0; d < n_detector; ++d)
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
    std::cout << "Outputting on level : " << level << std::endl;
    std::string filename =
            "raw_counts_" + Utilities::int_to_string(level) + ".dat";

    // Open a file
    std::ofstream myfile;
    std::string   sep;
    myfile.open(filename);
    myfile << "vertex_positions_x vertex_position_y vertex_position_z ";
    for (unsigned int i = 0; i < n_detector; ++i)
        myfile << "detector_" + Utilities::int_to_string(i) + " ";
    myfile << std::endl;
    sep = " ";

    // showing contents:
    for (auto it = dof_index_and_location.begin();
         it != dof_index_and_location.end();
         ++it)
    {
        for (unsigned d = 0; d < dim; ++d){
            myfile << it->second[d] << sep;
        }
        myfile << "\n";


        //for (unsigned int i = 0; i < detectors.size(); ++i)
        //myfile << nodal_counts[i][it->first] << sep;

        //myfile << "\n";
    }
    myfile.close();
}

template <int dim>
void
RPTFEMReconstruction<dim>::L2_project()
{
    std::cout << "Assigning detector positions" << std::endl;
    assign_detector_positions();
    //ajouté
    n_detector = detectors.size();
    std::cout << "Number of detectors identified : " << detectors.size()
              << std::endl;

    // flatten the triangulation
    Triangulation<dim> temp_triangulation;
    Triangulation<dim> flat_temp_triangulation;
    GridGenerator::cylinder(temp_triangulation, 0.1, 0.07);
    //To make the mesh finer in z diraction
    //GridGenerator::subdivided_cylinder(temp_triangulation, 3,0.1, 0.13);
    temp_triangulation.refine_global(2);

    GridGenerator::flatten_triangulation(temp_triangulation,
                                         flat_temp_triangulation);

    GridGenerator::convert_hypercube_to_simplex_mesh(flat_temp_triangulation,
                                                     triangulation);

    triangulation.set_all_manifold_ids(0);

    // ajouté
    std::cout << "n_cell:" << this->triangulation.n_cells() << std::endl;

    GridOut grid_out;
    {
        std::ofstream output_file("original_triangulation.vtk");
        grid_out.write_vtk(triangulation, output_file);
    }

    GridTools::rotate(1.57078, 1, triangulation);
    Tensor<1, dim> shift_vector({0, 0, 0.065});
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
    checkpoint();
    //output_raw_results_per_level();
    //trajectory();
}

Vector<double>
assemble_matrix_and_rhs(
        std::vector<std::vector<double>> &vertex_count,
        std::vector<double> &             experimental_count)
{
    // Assembling system_matrix
    Vector<double> reference_location;



    LAPACKFullMatrix<double> system_matrix;
    system_matrix.reinit(3);

    unsigned int detector_size = vertex_count.size();
    // std::cout << "Detector size is " << detector_size << std::endl;


    // TODO, this could be refactored to just three loops...
    double sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (-vertex_count[i][0] + vertex_count[i][1]) *
                 (-vertex_count[i][0] + vertex_count[i][1]);
    }
    system_matrix.set(0, 0, sigma);


    sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (-vertex_count[i][0] + vertex_count[i][2]) *
                 (-vertex_count[i][0] + vertex_count[i][1]);
    }
    system_matrix.set(0, 1, sigma);



    sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (-vertex_count[i][0] + vertex_count[i][3]) *
                 (-vertex_count[i][0] + vertex_count[i][1]);
    }
    system_matrix.set(0, 2, sigma);


    sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (-vertex_count[i][0] + vertex_count[i][2]) *
                 (-vertex_count[i][0] + vertex_count[i][1]);
    }
    system_matrix.set(1, 0, sigma);

    sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (-vertex_count[i][0] + vertex_count[i][2]) *
                 (-vertex_count[i][0] + vertex_count[i][2]);
    }

    system_matrix.set(1, 1, sigma);


    sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (-vertex_count[i][0] + vertex_count[i][3]) *
                 (-vertex_count[i][0] + vertex_count[i][2]);
    }
    system_matrix.set(1, 2, sigma);

    sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (-vertex_count[i][0] + vertex_count[i][3]) *
                 (-vertex_count[i][0] + vertex_count[i][1]);
    }
    system_matrix.set(2, 0, sigma);

    sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (-vertex_count[i][0] + vertex_count[i][3]) *
                 (-vertex_count[i][0] + vertex_count[i][2]);
    }
    system_matrix.set(2, 1, sigma);

    sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (-vertex_count[i][0] + vertex_count[i][3]) *
                 (-vertex_count[i][0] + vertex_count[i][3]);
    }
    system_matrix.set(2, 2, sigma);

    // std::cout << "Printing matrix " << std::endl;
    // system_matrix.print_formatted(std::cout);


    // Assembling system_rhs
    Vector<double> system_rhs;
    system_rhs.reinit(3);

    sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (vertex_count[i][0] - experimental_count[i]) *
                 (-vertex_count[i][0] + vertex_count[i][1]);
    }
    system_rhs[0] = -sigma;

    sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (vertex_count[i][0] - experimental_count[i]) *
                 (-vertex_count[i][0] + vertex_count[i][2]);
    }
    system_rhs[1] = -sigma;

    sigma = 0;
    for (unsigned int i = 0; i < detector_size; i++)
    {
        sigma += (vertex_count[i][0] - experimental_count[i]) *
                 (-vertex_count[i][0] + vertex_count[i][3]);
    }
    system_rhs[2] = -sigma;



    system_matrix.set_property(LAPACKSupport::general);
    system_matrix.compute_lu_factorization();
    system_matrix.solve(system_rhs);




    reference_location = system_rhs;
    return reference_location;
}



template <int dim>
void
RPTFEMReconstruction<dim>::find_cell(std::vector<double> experimental_count)
{
    // Loop over cell in the finest level, loop over detetors get nodal values,
    // solve first loop over cell

    // for now I define experimental counts manually


    double max_reconstruction_error = DBL_MAX;
    double max_cost_function=DBL_MAX;
    double comparing_cost_function;
    Point<dim> result;
    Point<dim> final_result;
    std::vector<double> check;




    const auto &cell_iterator = this->dof_handler.active_cell_iterators();
    for (const auto &cell : cell_iterator)
    {

        std::vector<std::vector<double>> count_from_all_detectors;

        for (unsigned int i = 0; i < n_detector; ++i)
        {
            std::vector<double> detectorCount;


            for (unsigned int v = 0; v < cell->n_vertices(); ++v)
            {
                auto dof_index = cell->vertex_dof_index(v, 0);
                detectorCount.push_back(nodal_counts[i][dof_index]);
            }
            count_from_all_detectors.push_back(detectorCount);
            detectorCount.clear();
        }


        Vector<double> reference_location =
                assemble_matrix_and_rhs(count_from_all_detectors, experimental_count);

        // Check if the location is a valid one
        std::vector<double> err_coordinates(4);
        //what about the ones that they meet all four conditions? We are considering they are not valid but they are so close!

        for (unsigned int i = 0; i < 3; ++i)
        {
            err_coordinates[i] = 0;
            if (reference_location[i] > 1)
                err_coordinates[i] = (reference_location[i] - 1);

            if (reference_location[i] < 0)
                err_coordinates[i] = (0 - reference_location[i]);
        }

        // fourth error block
        {
            err_coordinates[3] = 0;

            if ((1 - reference_location[0] - reference_location[1] -
                 reference_location[2]) < 0)
                err_coordinates[3] =
                        std::abs(1 - reference_location[0] - reference_location[1] -
                                 reference_location[2]);
            else if ((1 - reference_location[0] - reference_location[1] -
                      reference_location[2]) > 1)
                err_coordinates[3] =
                        1 - std::abs(1 - reference_location[0] - reference_location[1] -
                                     reference_location[2]);
        }


        double norm_error_coordinates = 0;
        for (unsigned int i = 0; i < 4; ++i)
        {
            norm_error_coordinates += err_coordinates[i] * err_coordinates[i];
        }



        if (norm_error_coordinates < 0.05)
        {
            //std::cout << "Printing reference solution " << std::endl;
            //reference_location.print(std::cout);
            Point<dim> real_location;
            //std::cout << "Printing vertex " << std::endl;
            for (unsigned int v = 0; v < cell->n_vertices(); ++v)
            {
                auto vertex_location = cell->vertex(v);
                //std::cout << vertex_location << std::endl;
                if (v == 0)
                    real_location = cell->vertex(v);
                else
                    real_location =
                            real_location + reference_location[v - 1] *
                                            (cell->vertex(v) - cell->vertex(0));
            }
            max_reconstruction_error = norm_error_coordinates;
            /*std::cout << "The reconstruction error is : "
                      << norm_error_coordinates << std::endl;*/
            //std::cout << "The real location is : " << real_location << std::endl;

            //Cost function must be calculated here

            double cost_function1=0;

            for (unsigned int i = 0; i < n_detector; ++i){

                double count1=0;

                count1+=(nodal_counts[i][cell->vertex_dof_index(0, 0)]*(1-reference_location[0]-reference_location[1]
                                                                        -reference_location[2]))+
                        (nodal_counts[i][cell->vertex_dof_index(1, 0)]*(reference_location[0]))+
                        (nodal_counts[i][cell->vertex_dof_index(2, 0)]*(reference_location[1]))+
                        (nodal_counts[i][cell->vertex_dof_index(3, 0)]*(reference_location[2]));

                count1=count1-experimental_count[i];
                cost_function1+=count1*count1;

            }
            //std::cout<<cost_function1<<std::endl;
            comparing_cost_function=cost_function1;
            result=real_location;
            if (comparing_cost_function<max_cost_function){
                max_cost_function=comparing_cost_function;
                final_result=result;
            }
        }
        count_from_all_detectors.clear();
    }

    //std::cout<<"this: "<<max_cost_function<<std::endl;
    std::cout<<final_result<<std::endl;


}



template <int dim>
void
RPTFEMReconstruction<dim>::trajectory()
{

    std::ifstream in("Test.txt");
    std::vector<std::vector<double> > v;

    if (in) {
        std::string line;

        while (std::getline(in, line)) {
            v.push_back(std::vector<double>());

            // Break down the row into column values
            std::stringstream split(line);
            double value;

            while (split >> value)
                v.back().push_back(value);
        }
    }
    for(int i=0; i<v.size();i++){

        std::vector<double> experimental_count=v[i];
        find_cell(experimental_count);


    }

}

template <int dim>
void
RPTFEMReconstruction<dim>::checkpoint()
{
    // save triangulation
    {
        std::ofstream ofs("temp_tria.tria");
        boost::archive::text_oarchive oa(ofs);
        triangulation.save(oa, 0);
    }

    // save dof_handler
    {
        std::ofstream ofs("temp_dof_handler.dof");
        boost::archive::text_oarchive oa(ofs);
        dof_handler.save(oa, 0);
    }

    // save nodal_counts_per_detector
    {
        for (unsigned int i = 0; i < detectors.size(); ++i)
        {
            std::string filename("temp_nodal_counts_detector" + Utilities::to_string(i,2) + ".counts");
            std::ofstream ofs(filename);
            boost::archive::text_oarchive oa(ofs);
            nodal_counts[i].save(oa, 0);
        }
    }


}


template <int dim>
void
RPTFEMReconstruction<dim>::load_from_checkpoint()
{
  n_detector = rpt_parameters.fem_reconstruction_param.nodal_counts_file.size();
  // flatten the triangulation
  Triangulation<dim> temp_triangulation;
  Triangulation<dim> flat_temp_triangulation;
  GridGenerator::cylinder(temp_triangulation, 0.1, 0.07);
  temp_triangulation.refine_global(2);

  GridGenerator::flatten_triangulation(temp_triangulation,
                                       flat_temp_triangulation);

  GridGenerator::convert_hypercube_to_simplex_mesh(flat_temp_triangulation,
                                                   triangulation);

  triangulation.set_all_manifold_ids(0);

  // import triangulation
  {
    //print number of cells
    std::cout << "n_cell:" << this->triangulation.n_cells() << std::endl;

    std::ifstream ifs(rpt_parameters.fem_reconstruction_param.triangulation_file);
    boost::archive::text_iarchive ia(ifs);
    triangulation.load(ia, 0);
  }

  // import dof handler
  {
    dof_handler.distribute_dofs(fe);
    std::ifstream ifs(rpt_parameters.fem_reconstruction_param.dof_handler_file);
    boost::archive::text_iarchive ia(ifs);
    dof_handler.load(ia, 0);
  }

  // import nodal counts
  {
    Vector<double>     counts_per_detector;

    nodal_counts.resize(n_detector);

    for (unsigned int i = 0; i < n_detector; ++i)
      {
        std::ifstream ifs(rpt_parameters.fem_reconstruction_param
                            .nodal_counts_file[i]);
        boost::archive::text_iarchive ia(ifs);
        counts_per_detector.load(ia, 0);
        nodal_counts[i] = counts_per_detector;
      }
  }
  // verification
  /*
  GridOut grid_out;
  std::ofstream output_file("loaded_rotated_triangulation.vtk");
  grid_out.write_vtk(triangulation, output_file);

  output_results();
  */
}

template <int dim>
void
RPTFEMReconstruction<dim>::rpt_fem_reconstruct()
{
  load_from_checkpoint();
  output_raw_results_per_level();
  trajectory();
}

template class RPTFEMReconstruction<3>;
