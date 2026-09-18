// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Deal.II includes
#include <deal.II/base/function_signed_distance.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/cell_data.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/interface_tools.h>
#include <core/parameters.h>
#include <core/simulation_control.h>

#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  /* This test checks the InterfaceTools::SignedDistanceSolver class from within
  an arbitrary structure that replicates the usual physical solver structure,
  denoted by the background prefix. This test performs verification of the
  signed distance computations for a sphere.
  */
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
  ConditionalOStream pcout(
        std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

  // Triangulation (as a shared_ptr to reproduce an arbitrary background solver
  // architecture)
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
    background_triangulation =
      std::make_shared<parallel::distributed::Triangulation<dim>>(
        mpi_communicator,
        typename Triangulation<dim>::MeshSmoothing(
          Triangulation<dim>::smoothing_on_refinement |
          Triangulation<dim>::smoothing_on_coarsening));

  Point<dim> p_0 = Point<dim>();
  Point<dim> p_1 = Point<dim>();
  for (int n = 0; n < dim; n++)
    {
      p_0[n] = 0.0;
      p_1[n] = 1.0;
    }
  GridGenerator::hyper_rectangle(*background_triangulation, p_0, p_1);
  background_triangulation->refine_global(3);


  // Discretize the domain as in an arbitrary background solver
  std::shared_ptr<FiniteElement<dim>> background_fe =
    std::make_shared<FE_Q<dim>>(1);
  std::shared_ptr<MappingQ<dim>> background_mapping =
    std::make_shared<MappingQ<dim>>(background_fe->degree);

  Vector<double> error(4);
  Vector<double> error_volume(4);

  
  for (unsigned int m = 0; m < 4; m++)
  {
    DoFHandler<dim> background_dof_handler =
      DoFHandler(*background_triangulation);
    background_dof_handler.distribute_dofs(*background_fe);

    IndexSet background_locally_owned_dofs =
      background_dof_handler.locally_owned_dofs();
    IndexSet background_locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(background_dof_handler);

    // Set the background level-set field of the sphere
    GlobalVectorType background_level_set(background_locally_owned_dofs,
                                          background_locally_relevant_dofs,
                                          mpi_communicator);


    const double sphere_radius = 0.23;
    Point<dim>   sphere_center = Point<dim>();
    for (int n = 0; n < dim; n++)
      {
        sphere_center[n] = 0.51;
      }

    auto sphere_function = Functions::SignedDistance::Sphere<dim>(sphere_center, sphere_radius);
    VectorTools::interpolate(
      *background_mapping,
      background_dof_handler,
      sphere_function,
      background_level_set);
      

    double h = 1.0/std::pow(2,(3+m));
    double band = 4.0*h;

    // Instanciate the SignedDistanceSolver as it would as a member of the
    // background solver class
    double max_reinitialization_distance = band*2;
    std::shared_ptr<InterfaceTools::SignedDistanceSolver<dim, GlobalVectorType>>
      signed_distance_solver = std::make_shared<
        InterfaceTools::SignedDistanceSolver<dim, GlobalVectorType>>(
        background_triangulation,
        background_fe,
        max_reinitialization_distance,
        0.0,
        1.0,
        Parameters::Verbosity::quiet);

    signed_distance_solver->setup_dofs();

    auto background_level_set_0 = background_level_set;

    GlobalVectorType background_signed_distance(background_locally_owned_dofs,
                                                    background_locally_relevant_dofs,
                                                    mpi_communicator);

    double volume_0, surface_0;

      std::tie(volume_0, surface_0) =
        InterfaceTools::integrate_volume_and_surface(background_dof_handler,
                                                     *background_fe,
                                                     background_level_set);
    

    std::vector<std::pair<double, std::string>> times_and_names;
    std::string filename = (dim == 2 ?
            "solution_2d_ref0" + std::to_string(m): "solution_3d_ref0" + std::to_string(m));

    for (int k = 0; k < 100; k++)     
      {
        signed_distance_solver->set_level_set_from_background_mesh(
          background_dof_handler, background_level_set);

        // Solve the signed_distance field.
        signed_distance_solver->solve();

        // Get the signed_distance field from the SignedDistanceSolver
        auto &signed_distance = signed_distance_solver->get_signed_distance();

        // Interpolate the signed_distance solution from the SignedDistanceSolver
        // DoFHandler to the main solver DoFHandler.
        

        TrilinosWrappers::MPI::Vector tmp_background_signed_distance(
          background_locally_owned_dofs, mpi_communicator);

        FETools::interpolate(signed_distance_solver->dof_handler,
                            signed_distance,
                            background_dof_handler,
                            tmp_background_signed_distance);

        background_signed_distance = tmp_background_signed_distance;

        DataOut<dim> data_out;
        data_out.attach_dof_handler(background_dof_handler);
        data_out.add_data_vector(background_signed_distance, "solution");
      
        data_out.build_patches();

        
        //   std::ofstream output(filename.c_str());

        const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record("./",filename, k,mpi_communicator,4);

       
        times_and_names.emplace_back(k, pvtu_filename);
        
        background_level_set = background_signed_distance;
      }

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          std::ofstream pvd_output(filename + ".pvd");
          DataOutBase::write_pvd_record(pvd_output, times_and_names);
        }

    double volume_1, surface_1;


      std::tie(volume_1, surface_1) =
        InterfaceTools::integrate_volume_and_surface(background_dof_handler,
                                                     *background_fe,
                                                     background_signed_distance);
    
    pcout << "volume_1 = " << volume_1 << ", volume_0 = " << volume_0 << std::endl;
    if (dim == 3)                                                 
      error_volume[m] = std::abs(volume_1 - 4.0 * M_PI * std::pow(sphere_radius, 3) / 3.0);
    
    if (dim == 2)
      error_volume[m] = std::abs(volume_1 - M_PI * std::pow(sphere_radius, 2));

    error_volume[m] = std::abs(volume_1 - volume_0);
    QGauss<dim>   quadrature(background_fe->degree + 1);

    FEValues<dim> fe_values(*background_mapping, *background_fe, quadrature, update_values | update_JxW_values | update_quadrature_points);
    Vector<double> cell_dof_values(background_fe->dofs_per_cell);
    Vector<double> cell_dof_values_analytical(background_fe->dofs_per_cell);

    std::vector<double> cell_quadrature_values(fe_values.n_quadrature_points);
    std::vector<double> cell_quadrature_values_analytical(fe_values.n_quadrature_points);

    double error_L2 = 0.0;
    double volume = 0.0;

    for (const auto &cell : background_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
            fe_values.reinit(cell);

            cell->get_dof_values(background_signed_distance,cell_dof_values.begin(),
                                      cell_dof_values.end());
            cell->get_dof_values(background_level_set_0,cell_dof_values_analytical.begin(),
                                      cell_dof_values_analytical.end());

            bool in_band = true;

            for (unsigned int j = 0; j < background_fe->dofs_per_cell; ++j)
              {
                if (std::abs(cell_dof_values[j]) > band)
                {
                  in_band = false;
                }
              }   
            
            if (!in_band)
              continue;

            fe_values.get_function_values(background_signed_distance, cell_quadrature_values);
            fe_values.get_function_values(background_level_set_0, cell_quadrature_values_analytical);

            for (const unsigned int q :
                 fe_values.quadrature_point_indices())
              {

                double extact_distance = sphere_function.value(fe_values.quadrature_point(q));
                // double extact_distance = cell_quadrature_values_analytical[q];

                error_L2 += (cell_quadrature_values[q]-extact_distance)*(cell_quadrature_values[q]-extact_distance)*fe_values.JxW(q);

                volume += fe_values.JxW(q);
              }
        }

    }

    error_L2 = Utilities::MPI::sum(error_L2, mpi_communicator);
    volume = Utilities::MPI::sum(volume, mpi_communicator);

    error_L2 = sqrt(error_L2/volume);

    

    // Compute the L2 norm of the error.
    Vector<float> error_per_cell(background_triangulation->n_active_cells());

    VectorTools::integrate_difference(
      *background_mapping,
      background_dof_handler,
      background_signed_distance,
      Functions::SignedDistance::Sphere<dim>(sphere_center, sphere_radius),
      error_per_cell,
      QGauss<dim>(background_fe->degree + 1),
      VectorTools::L2_norm);

    double error_L2_ff =
      VectorTools::compute_global_error(*background_triangulation,
                                        error_per_cell,
                                        VectorTools::L2_norm);

    // background_signed_distance.sadd(-1.0,background_level_set);

    // background_signed_distance.print(std::cout);

    // double error_L2 = sqrt(background_signed_distance.norm_sqr()/background_level_set.norm_sqr());


    deallog << "The L2 norm of the signed distance error in " << dim
            << "D is: " << error_L2 << std::endl;

    pcout << "The L2 norm of the signed distance error in " << dim
            << "D is: " << error_L2 << std::endl;
    pcout << "The volume error in " << dim
            << "D is: " << error_volume[m] << std::endl;

    error[m] = error_L2;

    

    background_triangulation->refine_global(1);
  }

  double convergence_order =
    log(error[3] / error[2]) / log(0.5);

  deallog << "The convergence is: " << convergence_order << std::endl;
  pcout << "The convergence is: " << convergence_order << std::endl;

  convergence_order =
    log(error_volume[3] / error_volume[2]) / log(0.5);
  pcout << "The convergence of the volume is: " << convergence_order << std::endl;

}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test<2>();
      test<3>();
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
