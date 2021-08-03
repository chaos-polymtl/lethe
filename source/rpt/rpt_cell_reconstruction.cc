
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
#include <rpt/rpt_cell_reconstruction.h>

#include <fstream>

template <int dim>
RPTCellReconstruction<dim>::RPTCellReconstruction(
  RPTCalculatingParameters &rpt_parameters)
  : parameters(rpt_parameters)
{
  // Seed the random number generator
  srand(parameters.rpt_param.seed);
}

template <int dim>
void
RPTCellReconstruction<dim>::execute_cell_reconstruction()
{
  assign_detector_positions(); // Assign detector position & create detectors
  read_counts(); // Read reconstruction counts of unknown particle positions
  create_grid(); // Read a grid for the reactor vessel
  set_coarse_mesh_counts(); // Calculate counts at vertices of the coarse mesh

  // Show positions in the terminal if verbose is enable
  for (unsigned int i = 0; i < reconstruction_counts.size(); i++)
    {
      find_unknown_position(reconstruction_counts[i]);
      if (parameters.rpt_param.verbose)
        std::cout << "Unknown particle : " << i
                  << " Cell volume : " << cells_volumes[i]
                  << " Position : " << reconstruction_positions[i] << std::endl;
    }

  // Export positions in file
  export_positions();
}

template <int dim>
void
RPTCellReconstruction<dim>::set_coarse_mesh_counts()
{
  find_vertices_positions(0);
  calculate_counts();
}

template <int dim>
void
RPTCellReconstruction<dim>::find_unknown_position(
  std::vector<double> &particle_reconstruction_counts)
{
  unsigned int     level = 0;
  std::vector<int> cells_indexes, new_cells_indexes;
  Point<dim>       reconstruction_position;
  double           volume;
  bool             still_new_candidates = true;

  // Find the first cell candidates from coarse mesh
  cells_indexes = find_cells(level, particle_reconstruction_counts);

  // Find the positions with the more refine candidate cell
  while (level < parameters.reconstruction_param.reactor_refinement &&
         still_new_candidates)
    {
      level++; // Start at level 1

      // Get vertices positions at current level and calculate their counts
      find_vertices_positions(level, cells_indexes);
      calculate_counts();

      // Find cell candidate for the position
      // (reconstruction counts are in the counts range of the cells for every
      // detector)
      new_cells_indexes =
        find_cells(level, particle_reconstruction_counts, cells_indexes);

      // If no refined cells are candidates at current level, research has to
      // stop and need to go back to parent level
      if (new_cells_indexes[0] == cells_indexes[0])
        {
          still_new_candidates = false;
          level--;
        }

      cells_indexes = new_cells_indexes;
    }

  // If many cells are candidate at the more refine level, the least squared
  // method is used to get the best cell
  int best_cell;
  if (cells_indexes.size() > 1)
    {
      best_cell =
        find_best_cell(level, particle_reconstruction_counts, cells_indexes);
    }
  else
    best_cell = cells_indexes[0];

  // Find the volume and the center position of the best cell
  for (const auto &cell : triangulation.cell_iterators_on_level(level))
    {
      if (cell->index() == best_cell)
        {
          volume                  = cell->measure();
          reconstruction_position = cell->center(true);
        }
    }

  // Store reconstruction positions and volume prior exportation
  reconstruction_positions.push_back(reconstruction_position);
  cells_volumes.push_back(volume);
}

template <int dim>
void
RPTCellReconstruction<dim>::create_grid()
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
  triangulation.refine_global(
    parameters.reconstruction_param.reactor_refinement);
}

template <int dim>
void
RPTCellReconstruction<dim>::find_vertices_positions(
  unsigned int     level,
  std::vector<int> parent_cell_indexes /* {-1} */)
{
  // Find the positions associated to the vertices of cell
  for (const auto &cell : triangulation.cell_iterators_on_level(level))
    {
      bool parent_is_candidate = false;

      // Check if parent cell of current level cell is a candidate
      if (cell->level() > 0)
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
          // Calculate counts at new vertices
          for (unsigned int i = 0; i < cell->n_vertices(); i++)
            {
              unsigned int vertex_index    = cell->vertex_index(i);
              Point<dim>   vertex_location = cell->vertex(i);

              // Check if vertex already exist and store the position with the
              // key as the vertex id
              if (map_vertices_index.find(vertex_index) ==
                  map_vertices_index.end())
                {
                  std::pair<Point<dim>, std::vector<double>> point_empty_vector(
                    vertex_location, {});
                  map_vertices_index[vertex_index] = point_empty_vector;
                }
            }
        }
    }
}

template <int dim>
void
RPTCellReconstruction<dim>::calculate_counts()
{
  // Vector fir calculated counts for one particle postiion
  std::vector<double> calculated_counts;

  unsigned int i = 0;
  for (auto &it : map_vertices_index)
    {
      // Check if counts are not already calculated at vertex
      if (it.second.second.empty())
        {
          calculated_counts.clear();
          RadioParticle<dim> particle(it.second.first, i);

          for (unsigned int i_detector = 0; i_detector < detectors.size();
               i_detector++)
            {
              // Create the particle-detector interaction object
              ParticleDetectorInteractions<dim> particle_detector_interactions(
                particle, detectors[i_detector], parameters);

              // Calculate count and store in vector for the position
              double count = particle_detector_interactions.calculate_count();
              calculated_counts.push_back(count);
            }
          it.second.second = calculated_counts;
          i++;
        }
    }
}

template <int dim>
std::vector<int>
RPTCellReconstruction<dim>::find_cells(
  unsigned int         level,
  std::vector<double> &particle_reconstruction_counts,
  std::vector<int>     parent_cell_indexes /* {-1} */)
{
  double            min, max;
  std::vector<bool> candidate_detector(detectors.size());
  std::vector<int>  candidates;

  // Check if parent cell of current level cell is a candidate
  for (const auto &cell : triangulation.cell_iterators_on_level(level))
    {
      bool parent_is_candidate = false;
      if (cell->level() > 0)
        {
          for (auto &id : parent_cell_indexes)
            {
              if (id == cell->parent_index())
                parent_is_candidate = true;
            }
        }
      else
        parent_is_candidate = true;

      if (parent_is_candidate)
        {
          for (unsigned int j = 0; j < detectors.size(); j++)
            {
              std::vector<double> counts_vertices(8);
              for (unsigned int i = 0; i < cell->n_vertices(); i++)
                {
                  // Get counts associate to the current vertex and detector
                  unsigned int vertex_index = cell->vertex_index(i);
                  counts_vertices[i] =
                    map_vertices_index[vertex_index].second[j];
                }
              // Get the min/max counts values
              min = *std::min_element(counts_vertices.begin(),
                                      counts_vertices.end());
              max = *std::max_element(counts_vertices.begin(),
                                      counts_vertices.end());

              // If counts is in the range, counts for the detector is true
              if (particle_reconstruction_counts[j] >= min &&
                  particle_reconstruction_counts[j] <= max)
                candidate_detector[j] = true;
              else
                candidate_detector[j] = false;
            }

          // If reconstruction counts are in ranges of cell for every detector,
          // cell in candidate
          if (std::all_of(candidate_detector.begin(),
                          candidate_detector.end(),
                          [](bool v) { return v; }))
            candidates.push_back(cell->index());
        }
    }

  // If there's no new candidates, parent cells are kept as candidate
  if (candidates.empty())
    candidates = parent_cell_indexes;

  return candidates;
}

template <int dim>
int
RPTCellReconstruction<dim>::find_best_cell(
  unsigned int         level,
  std::vector<double> &particle_reconstruction_counts,
  std::vector<int>     candidate)
{
  std::vector<std::pair<int, double>> cellid_error;
  double                              error;
  int                                 best_cell;

  // Use least squared error cell by cell for all counts per detector and
  // all vertices of the cell
  for (const auto &cell : triangulation.cell_iterators_on_level(level))
    {
      bool is_candidate = false;

      // Check if candidate cell
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
                               particle_reconstruction_counts[j] -
                               map_vertices_index[vertex_index].second[j]),
                             2);
                }
            }
          cellid_error.push_back(std::make_pair(cell->index(), error));
        }
    }

  // Sorted index of cell by error
  std::sort(cellid_error.begin(),
            cellid_error.end(),
            [](std::pair<int, double> &pair1, std::pair<int, double> &pair2) {
              return (pair1.second < pair2.second);
            });

  // Return the best cell index
  return best_cell = cellid_error[0].first;
}

template <int dim>
void
RPTCellReconstruction<dim>::read_counts()
{
  // Read counts for reconstruction
  const std::string filename =
    parameters.reconstruction_param.reconstruction_counts_file;
  std::ifstream counts_file(filename);

  std::vector<double> values;
  std::copy(std::istream_iterator<double>(counts_file),
            std::istream_iterator<double>(),
            std::back_inserter(values));

  // Get the number of detector
  unsigned int number_of_detectors = detectors.size();

  // Extract positions, create point objects and detectors
  for (unsigned int i = 0; i < values.size(); i += number_of_detectors)
    {
      std::vector<double>::const_iterator first = values.begin() + i;
      std::vector<double>::const_iterator last =
        values.begin() + i + number_of_detectors;

      std::vector<double> particle_counts(first, last);

      reconstruction_counts.push_back(particle_counts);
    }
}

template <int dim>
void
RPTCellReconstruction<dim>::export_positions()
{
  std::ofstream myfile;
  std::string   sep;
  std::string   filename =
    parameters.reconstruction_param.reconstruction_positions_file;
  myfile.open(filename);
  if (filename.substr(filename.find_last_of(".") + 1) == ".dat")
    {
      myfile
        << "unknown_positions_id volumes particle_positions_x particle_positions_y particle_positions_z"
        << std::endl;
      sep = " ";
    }
  else // .csv is default
    {
      myfile
        << "unknown_positions_id,volumes,particle_positions_x,particle_positions_y,particle_positions_z"
        << std::endl;
      sep = ",";
    }

  for (unsigned int i_particle = 0;
       i_particle < reconstruction_positions.size();
       i_particle++)
    {
      myfile << i_particle << sep << cells_volumes[i_particle] << sep
             << reconstruction_positions[i_particle][0] << sep
             << reconstruction_positions[i_particle][1] << sep
             << reconstruction_positions[i_particle][2] << std::endl;
    }
}

template <int dim>
void
RPTCellReconstruction<dim>::assign_detector_positions()
{
  // Read text file with detector positions and store it in vector
  std::ifstream detector_file(
    parameters.detector_param.detector_positions_file);

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

      Detector<dim> detector(parameters.detector_param,
                             i,
                             face_point,
                             middle_point);

      detectors.push_back(detector);
    }
}


template class RPTCellReconstruction<3>;