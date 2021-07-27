#include <core/solutions_output.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/data_out.h>

#include <rpt/rpt_map.h>

#include <iostream>

template <int dim>
RPTMap<dim>::RPTMap(std::vector<Detector<dim>> &detectors,
                    RPTCalculatingParameters &  rpt_parameters)
  : fe(FE_Q<dim>(1), detectors.size())
  , mapping(1, 0)
  , cell_quadrature(1)
  , face_quadrature(1 + 1)
  , parameters(rpt_parameters.rpt_param)
  , reconstruction_parameters(rpt_parameters.reconstruction_param)
  , detectors_positions(detectors)
{}

template <int dim>
void
RPTMap<dim>::create_grid()
{
  // Generate cylinder (needs rotation and shift to get origin at the bottom
  // with z towards top)
  GridGenerator::cylinder(triangulation,
                          parameters.reactor_radius,
                          parameters.reactor_height / 2.);
  GridTools::rotate(M_PI_2, 1, triangulation);
  const Tensor<1, dim> shift_vector({0, 0, parameters.reactor_height / 2.});
  GridTools::shift(shift_vector, triangulation);

  // Add cylindrical manifold
  const Tensor<1, dim>                direction({0, 0, 1});
  const CylindricalManifold<dim, dim> manifold(direction, {0, 0, 1});
  triangulation.set_manifold(0, manifold);

  // Refine the grid !!!!!!!!!!! Fix for MPI
  // triangulation->set_all_refine_flags();
  triangulation.prepare_coarsening_and_refinement();
  triangulation.refine_global(reconstruction_parameters.reactor_refinement);
}

template <int dim>
std::vector<Point<dim>>
RPTMap<dim>::get_positions()
{
  // Associate triangulation to dof_handler and distribute fe
  dof_handler.clear();
  dof_handler.reinit(triangulation);
  dof_handler.distribute_dofs(fe);

  // Find the positions associate to the sorted dofs and keep unique positions
  // in order to generate the vector of counts organized for data out.
  std::vector<std::pair<types::global_dof_index, Point<dim>>>
    dof_index_location;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      for (unsigned int i = 0; i < cell->n_vertices(); i++)
        {
          types::global_dof_index dof_index    = cell->vertex_dof_index(i, 0);
          Point<dim>              dof_location = cell->vertex(i);

          dof_index_location.push_back(std::make_pair(dof_index, dof_location));
        }
    }

  std::sort(dof_index_location.begin(),
            dof_index_location.end(),
            [](std::pair<types::global_dof_index, Point<dim>> &pair1,
               std::pair<types::global_dof_index, Point<dim>> &pair2) {
              return (pair1.first < pair2.first);
            });

  auto last = std::unique(dof_index_location.begin(), dof_index_location.end());
  dof_index_location.erase(last, dof_index_location.end());

  for (unsigned int i = 0; i < dof_index_location.size(); i++)
    positions.push_back(dof_index_location[i].second);

  return positions;
}

template <int dim>
void
RPTMap<dim>::output_results()
{
  FEValues<dim> fe_values(mapping,
                          fe,
                          cell_quadrature,
                          update_values | update_quadrature_points |
                            update_gradients);

  const unsigned int                   dofs_per_cell = fe.n_dofs_per_cell();
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


  std::vector<std::string> detector_ids = {"detector_0"};
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      detectors_positions.size(),
      DataComponentInterpretation::component_is_scalar);
  if (detectors_positions.size() > 1)
    {
      for (unsigned int i = 1; i < detectors_positions.size(); i++)
        detector_ids.push_back("detector_" + std::to_string(i));
    }

  DataOut<dim> data_out;
  // Attach the solution data to data_out object
  data_out.attach_dof_handler(dof_handler);

  data_out.add_data_vector(counts,
                           detector_ids,
                           DataOut<dim>::type_dof_data,
                           data_component_interpretation);

  data_out.build_patches(mapping);

  write_vtu_and_pvd(
    this->pvd_handler, data_out, "./", "test_output", 0, 0, 1, MPI_COMM_WORLD);
}

template <int dim>
void
RPTMap<dim>::add_calculated_counts(std::vector<double> calculated_counts)
{
  Vector<double> solution(calculated_counts.begin(), calculated_counts.end());
  counts = solution;

  output_results();
}


template class RPTMap<3>;