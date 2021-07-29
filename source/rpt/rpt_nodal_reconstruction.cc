
#include <core/solutions_output.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/data_out.h>

#include <rpt/particle_detector_interactions.h>
#include <rpt/rpt.h>
#include <rpt/rpt_nodal_reconstruction.h>

#include <fstream>

template <int dim>
RPTNodalReconstruction<dim>::RPTNodalReconstruction(
  std::vector<Detector<dim>> &detectors,
  RPTCalculatingParameters   &rpt_parameters)
  : fe(FE_Q<dim>(1), detectors.size())
  , mapping(1, 0)
  , cell_quadrature(1)
  , face_quadrature(1 + 1)
  , parameters(rpt_parameters)
  , reconstruction_parameters(rpt_parameters.reconstruction_param)
  , detectors(detectors)
{
  // Read counts for reconstruction
  const std::string filename =
    rpt_parameters.reconstruction_param.reconstruction_counts_file;
  std::ifstream counts_file(filename);

  std::copy(std::istream_iterator<double>(counts_file),
            std::istream_iterator<double>(),
            std::back_inserter(reconstruction_counts));
}

template <int dim>
void
RPTNodalReconstruction<dim>::execute_nodal_reconstruction()
{
  create_grid();
  find_unknown_position();
}

template <int dim>
void
RPTNodalReconstruction<dim>::find_unknown_position()
{
  unsigned int     level             = 0;
  std::vector<int> closer_cell_index = {-1, 0};
  double           volume            = 1;
  Point<dim>       pos;
  double           min_volume = pow(0.001, 3);

  while (volume > min_volume &&
         level < parameters.reconstruction_param.reactor_refinement)
    {
      IteratorRange<TriaIterator<CellAccessor<dim, dim>>> cell_iterators =
        triangulation.cell_iterators_on_level(level);
      get_positions(cell_iterators, closer_cell_index);
      calculate_counts(map_vertices_index);
      closer_cell_index = find_closer_cell(cell_iterators, closer_cell_index);
      for (const auto &cell : cell_iterators)
        {
          if (cell->index() == closer_cell_index[0])
            {
              volume = cell->measure();
              pos    = cell->center(true);
              std::cout << "Volume : " << volume << " position : " << pos
                        << std::endl;
            }
        }
      level++;
    }
}

template <int dim>
void
RPTNodalReconstruction<dim>::create_grid()
{
  // Generate cylinder (needs rotation and shift to get origin at the bottom
  // with z towards top)
  GridGenerator::cylinder(triangulation,
                          parameters.rpt_param.reactor_radius,
                          parameters.rpt_param.reactor_height / 2.);
  GridTools::rotate(M_PI_2, 1, triangulation);
  const Tensor<1, dim> shift_vector(
    {0, 0, parameters.rpt_param.reactor_height / 2.});
  GridTools::shift(shift_vector, triangulation);

  // Add cylindrical manifold
  const Tensor<1, dim>                direction({0, 0, 1});
  const CylindricalManifold<dim, dim> manifold(direction, {0, 0, 1});

  // Set manifold and refine global
  triangulation.set_manifold(0, manifold);
  triangulation.prepare_coarsening_and_refinement();
  triangulation.refine_global(reconstruction_parameters.reactor_refinement);
}

template <int dim>
std::vector<Point<dim>>
RPTNodalReconstruction<dim>::get_positions(
  IteratorRange<TriaIterator<CellAccessor<dim, dim>>> &cell_iterators,
  std::vector<int> parent_cell_indexes /* {-1, 0, 0} */)
{
  // Find the positions associated to the vertices and keep unique positions
  // in order to generate the map vertex_index-position.
  std::vector<std::pair<unsigned int, Point<dim>>> vertices_index_location;

  for (const auto &cell : cell_iterators)
    {
      if (parent_cell_indexes[0] == -1 ||
          cell->parent()->index() == parent_cell_indexes[0] ||
          cell->parent()->index() == parent_cell_indexes[1])
        {
          for (unsigned int i = 0; i < cell->n_vertices(); i++)
            {
              unsigned int vertex_index    = cell->vertex_index(i);
              Point<dim>   vertex_location = cell->vertex(i);

              vertices_index_location.push_back(
                std::make_pair(vertex_index, vertex_location));
            }
        }
    }

  std::sort(vertices_index_location.begin(),
            vertices_index_location.end(),
            [](std::pair<unsigned int, Point<dim>> &pair1,
               std::pair<unsigned int, Point<dim>> &pair2) {
              return (pair1.first < pair2.first);
            });

  // Keep unique position to avoid redundant counts calculation
  auto last =
    std::unique(vertices_index_location.begin(), vertices_index_location.end());
  vertices_index_location.erase(last, vertices_index_location.end());

  // Put vertices indexes to a container and position in a vector
  map_vertices_index.clear();
  for (unsigned int i = 0; i < vertices_index_location.size(); i++)
    {
      map_vertices_index[vertices_index_location[i].first];
      positions.push_back(vertices_index_location[i].second);
      std::cout << vertices_index_location[i].first << " "
                << vertices_index_location[i].second << std::endl;
    }

  return positions;
}

template <int dim>
std::vector<int>
RPTNodalReconstruction<dim>::find_closer_cell(
  IteratorRange<TriaIterator<CellAccessor<dim, dim>>> &cell_iterators,
  std::vector<int> parent_cell_indexes /* {-1, 0, 0} */)
{
  std::vector<std::pair<int, double>> cellid_error;

  // Calculate least squared error cell by cell for all counts per detector and
  // all vertices of the cell
  for (const auto &cell : cell_iterators)
    {
      if (parent_cell_indexes[0] == -1 ||
          cell->parent()->index() == parent_cell_indexes[0] ||
          cell->parent()->index() == parent_cell_indexes[1])
        {
          double error = 0;
          for (unsigned int i = 0; i < cell->n_vertices(); i++)
            {
              for (unsigned int j = 0; j < detectors.size(); j++)
                {
                  unsigned int vertex_index = cell->vertex_index(i);
                  error +=
                    std::pow(std::fabs(reconstruction_counts[j] -
                                       map_vertices_index[vertex_index][j]),
                             2);
                }
            }
          cellid_error.push_back(std::make_pair(cell->index(), error));
        }
    }


  // Sorted index of cell
  std::sort(cellid_error.begin(),
            cellid_error.end(),
            [](std::pair<int, double> &pair1, std::pair<int, double> &pair2) {
              return (pair1.second < pair2.second);
            });

  for (unsigned int i = 0; i < cellid_error.size(); i++)
    std::cout << cellid_error[i].first << " " << cellid_error[i].second
              << std::endl;

  std::vector<int> cellids{cellid_error[0].first, cellid_error[1].first};

  return cellids;
}

/*
template <int dim>
void
RPTNodalReconstruction<dim>::output_results()
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
    detectors.size(),
    DataComponentInterpretation::component_is_scalar);
  if (detectors.size() > 1)
  {
    for (unsigned int i = 1; i < detectors.size(); i++)
      detector_ids.push_back("detector_" + std::to_string(i));
  }

  DataOut<dim> data_out;
  // Attach the solution data to data_out object
  data_out.attach_dof_handler(dof_handler);

  data_out.add_data_vector(map_counts,
                           detector_ids,
                           DataOut<dim>::type_dof_data,
                           data_component_interpretation);

  data_out.build_patches(mapping);

  write_vtu_and_pvd(
    this->pvd_handler, data_out, "./", "test_output", 0, 0, 1, MPI_COMM_WORLD);
} */


template <int dim>
void
RPTNodalReconstruction<dim>::calculate_counts(
  std::map<unsigned int, std::vector<double>> &index_count)
{
  std::vector<double> calculated_counts;
  // Create radioactive particle with positions
  particle_positions.clear();
  for (unsigned int i = 0; i < positions.size(); i++)
    {
      RadioParticle<dim> particle(positions[i], i);
      particle_positions.push_back(particle);
    }


  // Calculate count for every particle-detector pair
  for (unsigned int i_particle = 0; i_particle < positions.size(); i_particle++)
    {
      for (unsigned int i_detector = 0; i_detector < detectors.size();
           i_detector++)
        {
          // Create the particle-detector interaction object
          ParticleDetectorInteractions<dim> particle_detector_interactions(
            particle_positions[i_particle], detectors[i_detector], parameters);


          double count = particle_detector_interactions.calculate_count();
          calculated_counts.push_back(count);
        }
    }

  // Transfer counts to the map of vertices id
  unsigned int i = 0;
  for (auto &it : index_count)
    {
      std::vector<double>::const_iterator first = calculated_counts.begin() + i;
      std::vector<double>::const_iterator last =
        calculated_counts.begin() + i + detectors.size();
      std::vector<double> counts(first, last);
      it.second = counts;
      i++;
    }
}

template class RPTNodalReconstruction<3>;