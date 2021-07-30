
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
  set_coarse_mesh_map();
  find_unknown_position();
}

template <int dim>
void
RPTNodalReconstruction<dim>::set_coarse_mesh_map()
{
  unsigned int coarse_level = 1;

  get_positions(triangulation.cell_iterators_on_level(coarse_level));
  calculate_counts(map_vertices_index);
  map_vertices_index_coarse_mesh = map_vertices_index;
  cells_indexes_coarse_mesh =
    find_cells(triangulation.cell_iterators_on_level(coarse_level));
}

template <int dim>
void
RPTNodalReconstruction<dim>::find_unknown_position()
{
  unsigned int     level = 1;
  std::vector<int> cells_indexes, new_cells_indexes;
  double           volume = 1;
  Point<dim>       pos;
  double           min_volume = parameters.reconstruction_param.minimum_volume;
  bool             still_new_candidates = true;

  cells_indexes = cells_indexes_coarse_mesh;

  while (volume > min_volume &&
         level <= parameters.reconstruction_param.reactor_refinement &&
         still_new_candidates)
    {
      level++;
      get_positions(triangulation.cell_iterators_on_level(level),
                    cells_indexes);
      calculate_counts(map_vertices_index);
      new_cells_indexes =
        find_cells(triangulation.cell_iterators_on_level(level), cells_indexes);

      // If no refined cells are candidates research has to stop and need to
      // have the cell iterators for the level of parents
      if (new_cells_indexes[0] == cells_indexes[0])
        {
          still_new_candidates = false;
          level--;
        }

      for (const auto &cell : triangulation.cell_iterators_on_level(level))
        {
          if (cell->index() == new_cells_indexes[0])
            {
              volume = cell->measure();
              pos    = cell->center(true);
            }
        }
      cells_indexes = new_cells_indexes;
    }

  // If after reaching criteria there's many candidate, the least squared method
  // is used
  int best_cell;
  if (cells_indexes.size() > 1)
    {
      best_cell = find_best_cell(triangulation.cell_iterators_on_level(level),
                                 cells_indexes);
      for (const auto &cell : triangulation.cell_iterators_on_level(level))
        {
          if (cell->index() == best_cell)
            {
              volume = cell->measure();
              pos    = cell->center(true);
            }
        }
    }

  std::cout << "Volume : " << volume << " position : " << pos << std::endl;
}

template <int dim>
void
RPTNodalReconstruction<dim>::create_grid()
{
  // Generate cylinder (needs rotation and shift to get origin at the bottom
  // with z towards top)
  GridGenerator::subdivided_cylinder(triangulation,
                                     5,
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
void
RPTNodalReconstruction<dim>::get_positions(
  IteratorRange<TriaIterator<CellAccessor<dim, dim>>> cell_iterators,
  std::vector<int> parent_cell_indexes /* {-1} */)
{
  // Find the positions associated to the vertices and keep unique positions
  // in order to generate the map vertex_index-position.
  std::vector<std::pair<unsigned int, Point<dim>>> vertices_index_location;

  for (const auto &cell : cell_iterators)
    {
      bool parent_is_candidate = false;
      if (cell->level() > 1)
        {
          for (auto &id : parent_cell_indexes)
            {
              if (id == cell->parent()->index())
                parent_is_candidate = true;
            }
        }
      else
        parent_is_candidate = true; // Level 0 : no parent

      if (parent_is_candidate)
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
  for (auto &it_vertices : vertices_index_location)
    {
      std::pair<Point<dim>, std::vector<double>> point_empty_vector(
        it_vertices.second, {});
      map_vertices_index.insert({it_vertices.first, point_empty_vector});

      std::cout << it_vertices.first << " "
                << map_vertices_index[it_vertices.first].first << std::endl;
    }
}

template <int dim>
void
RPTNodalReconstruction<dim>::calculate_counts(
  std::map<unsigned int, std::pair<Point<dim>, std::vector<double>>>
    &index_count)
{
  std::vector<double> calculated_counts;
  // Create radioactive particle with positions
  particles.clear();
  unsigned int id = 0;
  for (auto &it_map : index_count)
    {
      RadioParticle<dim> particle(it_map.second.first, id);
      particles.push_back(particle);
      id++;
    }


  // Calculate count for every particle-detector pair
  for (unsigned int i_particle = 0; i_particle < particles.size(); i_particle++)
    {
      for (unsigned int i_detector = 0; i_detector < detectors.size();
           i_detector++)
        {
          // Create the particle-detector interaction object
          ParticleDetectorInteractions<dim> particle_detector_interactions(
            particles[i_particle], detectors[i_detector], parameters);

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
      it.second.second = counts;
      i += detectors.size();
    }
}

template <int dim>
std::vector<int>
RPTNodalReconstruction<dim>::find_cells(
  IteratorRange<TriaIterator<CellAccessor<dim, dim>>> cell_iterators,
  std::vector<int> parent_cell_indexes /* {-1} */)
{
  double            min, max;
  std::vector<bool> candidate_detector(detectors.size());
  std::vector<int>  candidates;

  // Calculate least squared error cell by cell for all counts per detector and
  // all vertices of the cell
  for (const auto &cell : cell_iterators)
    {
      bool parent_is_candidate = false;
      if (cell->level() > 1)
        {
          for (auto &id : parent_cell_indexes)
            {
              if (id == cell->parent_index())
                parent_is_candidate = true;
            }
        }
      else
        parent_is_candidate = true; // Level 0 : no parent

      if (parent_is_candidate)
        {
          for (unsigned int j = 0; j < detectors.size(); j++)
            {
              std::vector<double> counts_vertices(8);
              for (unsigned int i = 0; i < cell->n_vertices(); i++)
                {
                  unsigned int vertex_index = cell->vertex_index(i);
                  counts_vertices[i] =
                    map_vertices_index[vertex_index].second[j];
                }
              min = *std::min_element(counts_vertices.begin(),
                                      counts_vertices.end());
              max = *std::max_element(counts_vertices.begin(),
                                      counts_vertices.end());
              if (reconstruction_counts[j] >= min &&
                  reconstruction_counts[j] <= max)
                candidate_detector[j] = true;
              else
                candidate_detector[j] = false;
            }
          if (std::all_of(candidate_detector.begin(),
                          candidate_detector.end(),
                          [](bool v) { return v; }))
            candidates.push_back(cell->index());
        }
    }


  // If there's no children candidates, parent cells are keep
  if (candidates.empty())
    candidates = parent_cell_indexes;

  for (unsigned int i = 0; i < candidates.size(); i++)
    std::cout << candidates[i] << std::endl;

  return candidates;
}

template <int dim>
int
RPTNodalReconstruction<dim>::find_best_cell(
  IteratorRange<TriaIterator<CellAccessor<dim, dim>>> cell_iterators,
  std::vector<int>                                    candidate)
{
  std::vector<std::pair<int, double>> cellid_error;
  double                              error;
  int                                 best_cell;

  // Calculate least squared error cell by cell for all counts per detector and
  // all vertices of the cell
  for (const auto &cell : cell_iterators)
    {
      bool is_candidate = false;
      for (auto &id : candidate)
        {
          if (id == cell->index())
            is_candidate = true;
        }
      if (is_candidate)
        {
          error = 0;
          for (unsigned int j = 0; j < detectors.size(); j++)
            {
              for (unsigned int i = 0; i < cell->n_vertices(); i++)
                {
                  unsigned int vertex_index = cell->vertex_index(i);
                  error +=
                    std::pow(std::fabs(
                               reconstruction_counts[j] -
                               map_vertices_index[vertex_index].second[j]),
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

  return best_cell = cellid_error[0].first;
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



template class RPTNodalReconstruction<3>;