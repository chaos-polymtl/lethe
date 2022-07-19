#include <deal.II/grid/cell_id.h>

#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/data_out.h>

#include <rpt/particle_detector_interactions.h>
#include <rpt/particle_visualization.h>
#include <rpt/rpt_cell_reconstruction.h>
#include <rpt/rpt_utilities.h>

#include <cmath>
#include <fstream>

template <int dim>
RPTCellReconstruction<dim>::RPTCellReconstruction(
  Parameters::RPTParameters               &rpt_parameters,
  Parameters::RPTReconstructionParameters &rpt_reconstruction_parameters,
  Parameters::DetectorParameters          &rpt_detector_parameters)
  : parameters(rpt_parameters)
  , reconstruction_parameters(rpt_reconstruction_parameters)
  , detector_parameters(rpt_detector_parameters)
  , computing_timer(std::cout, TimerOutput::summary, TimerOutput::wall_times)
{
  // Seed the random number generator
  srand(parameters.seed);
}

template <int dim>
void
RPTCellReconstruction<dim>::execute_cell_reconstruction()
{
  // Assign detector position & create detectors
  detectors = assign_detector_positions<dim>(detector_parameters);

  // Read reconstruction (experimental) counts of unknown particle positions
  reconstruction_counts = read_detectors_counts<dim>(
    reconstruction_parameters.reconstruction_counts_file, detectors.size());

  // Read known positions if analyse positions if enable
  if (reconstruction_parameters.analyse_positions)
    {
      known_positions =
        read_positions<dim>(reconstruction_parameters.known_positions_file);
      AssertThrow(
        known_positions.size() == reconstruction_counts.size(),
        ExcMessage(
          "Reconstruction counts and known positions files do not have data for the same number of positions."))
    }

  // Create a grid for the reactor vessel
  attach_grid_to_triangulation<dim>(
    triangulation, parameters, reconstruction_parameters.reactor_refinement);

  // Calculate counts at vertices of the coarse mesh
  set_coarse_mesh_counts();
  find_all_positions();

  // Export positions in file
  export_positions();

  // Show information about cells and vertices
  if (parameters.verbosity == Parameters::Verbosity::verbose)
    {
      std::cout << std::endl;
      std::cout << "Number of active cells : " << triangulation.n_active_cells()
                << std::endl;
      std::cout << "Number of vertices in the mesh : "
                << triangulation.n_vertices() << std::endl;
      std::cout << "Number of vertices with calculated counts : "
                << map_vertices_index.size() << std::endl;
    }

  // Generate visualisation files
  visualize_positions();

  // Disable the output of time clock
  if (!reconstruction_parameters.verbose_clock)
    computing_timer.disable_output();
}

template <int dim>
void
RPTCellReconstruction<dim>::set_coarse_mesh_counts()
{
  find_vertices_positions(reconstruction_parameters.coarse_mesh_level);
  calculate_counts();
}

template <int dim>
void
RPTCellReconstruction<dim>::find_all_positions()
{
  // Show positions in the terminal if verbose is enable
  for (unsigned int i = 0; i < reconstruction_counts.size(); i++)
    {
      find_unknown_position(reconstruction_counts[i], i);
      if (parameters.verbosity == Parameters::Verbosity::verbose)
        std::cout << "Unknown particle : " << i
                  << " Cell volume : " << cells_volumes[i]
                  << " Position : " << reconstruction_positions[i] << std::endl;
    }
}

template <int dim>
void
RPTCellReconstruction<dim>::find_unknown_position(
  std::vector<double> &particle_reconstruction_counts,
  unsigned int         id)
{
  int              level = 0;
  std::vector<int> cells_indexes, new_cells_indexes;
  Point<dim>       reconstruction_position;
  double           volume;
  bool             still_new_candidates = true;
  std::string      status;

  // Find the first cell candidates from coarse mesh
  cells_indexes = find_cells(level, particle_reconstruction_counts);

  // Find the positions with the more refine candidate cell
  while (level < reconstruction_parameters.reactor_refinement &&
         still_new_candidates)
    {
      level++; // Start at coarse mesh level + 1

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
          status = "parent_cell"; // Cell isn't the highest refinement level
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

      // Update status
      if (status.empty())
        status = "cost_function";
      else
        status += "|cost_function";
    }
  else
    best_cell = cells_indexes[0];

  // Find the volume and the center position of the best cell and modify
  // status when analyse positions is enable
  for (const auto &cell : triangulation.cell_iterators_on_level(level))
    {
      if (cell->index() == best_cell)
        {
          volume                  = cell->measure();
          reconstruction_position = cell->center(true);

          // Verify if the best cell is the good cell and update status
          if (reconstruction_parameters.analyse_positions)
            {
              if (cell->point_inside(known_positions[id]))
                {
                  if (status.empty())
                    status = "right_cell";
                  else
                    status += "|right_cell";
                }
              else
                {
                  if (status.empty())
                    status = "wrong_cell";
                  else
                    status += "|wrong_cell";
                }
            }
        }
    }

  // Store information prior data exportation
  reconstruction_positions.push_back(reconstruction_position);
  cells_volumes.push_back(volume);
  final_cell_level.push_back(level);
  cell_status.push_back(status);
}

template <int dim>
void
RPTCellReconstruction<dim>::find_vertices_positions(
  unsigned int     level,
  std::vector<int> parent_cell_indexes /* {-1} */)
{
  TimerOutput::Scope t(computing_timer, "finding_vertices_positions");

  // Find the positions associated to the vertices of cell
  for (const auto &cell : triangulation.cell_iterators_on_level(level))
    {
      // Stay false if
      bool parent_is_candidate = false;

      // Check if parent cell of current level cell is a candidate
      // Stay false at coarsh mesh or if no candidate at refined level
      if (cell->level() > reconstruction_parameters.coarse_mesh_level)
        {
          for (auto &id : parent_cell_indexes)
            {
              if (id == cell->parent_index())
                parent_is_candidate = true;
            }
        }
      else
        parent_is_candidate = true; // Coarse mesh level : no parent

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
  TimerOutput::Scope t(computing_timer, "vertices_calculation_counts");

  // Vector fir calculated counts for one particle position
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
  TimerOutput::Scope t(computing_timer, "finding_cell");

  double            min, max;
  std::vector<bool> candidate_detector(detectors.size());
  std::vector<int>  candidates;

  // Check if parent cell of current level cell is a candidate
  for (const auto &cell : triangulation.cell_iterators_on_level(level))
    {
      bool parent_is_candidate = false;
      if (cell->level() > reconstruction_parameters.coarse_mesh_level)
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
  TimerOutput::Scope t(computing_timer, "finding_cell");

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
RPTCellReconstruction<dim>::export_positions()
{
  std::ofstream myfile;
  std::string   sep;
  std::string   filename =
    reconstruction_parameters.reconstruction_positions_file;
  myfile.open(filename);

  // Headers with space or comma
  if (filename.substr(filename.find_last_of(".") + 1) == ".dat")
    {
      myfile
        << "reconstructed_positions_id volumes particle_positions_x particle_positions_y particle_positions_z cell_level";
      sep = " ";
    }
  else // .csv is default
    {
      myfile
        << "reconstructed_positions_id,volumes,particle_positions_x,particle_positions_y,particle_positions_z,cell_level";
      sep = ",";
    }

  // Add status columns and distance if analyse positions is enable
  if (reconstruction_parameters.analyse_positions)
    myfile << sep << "status" << sep << "status_id" << sep << "distance"
           << std::endl;
  else
    myfile << std::endl;

  for (unsigned int i_particle = 0;
       i_particle < reconstruction_positions.size();
       i_particle++)
    {
      myfile << i_particle << sep << cells_volumes[i_particle] << sep
             << reconstruction_positions[i_particle][0] << sep
             << reconstruction_positions[i_particle][1] << sep
             << reconstruction_positions[i_particle][2] << sep
             << final_cell_level[i_particle];

      if (reconstruction_parameters.analyse_positions)
        {
          // Calculate distance between known and found positions
          double distance =
            std::sqrt(std::pow(reconstruction_positions[i_particle][0] -
                                 known_positions[i_particle][0],
                               2) +
                      std::pow(reconstruction_positions[i_particle][1] -
                                 known_positions[i_particle][1],
                               2) +
                      std::pow(reconstruction_positions[i_particle][2] -
                                 known_positions[i_particle][2],
                               2));

          // Associate status to a status id
          unsigned int status_id = UINT_MAX;
          if (cell_status[i_particle] == "right_cell")
            status_id = 0;
          else if (cell_status[i_particle] == "cost_function|right_cell")
            status_id = 1;
          else if (cell_status[i_particle] == "parent_cell|right_cell")
            status_id = 2;
          else if (cell_status[i_particle] ==
                   "parent_cell|cost_function|right_cell")
            status_id = 3;
          else if (cell_status[i_particle] == "wrong_cell")
            status_id = 4;
          else if (cell_status[i_particle] == "cost_function|wrong_cell")
            status_id = 5;
          else if (cell_status[i_particle] == "parent_cell|wrong_cell")
            status_id = 6;
          else if (cell_status[i_particle] ==
                   "parent_cell|cost_function|wrong_cell")
            status_id = 7;


          myfile << sep << cell_status[i_particle] << sep << status_id << sep
                 << distance << std::endl;
        }
      else
        myfile << std::endl;
    }
}


template <int dim>
void
RPTCellReconstruction<dim>::visualize_positions()
{
  TimerOutput::Scope t(computing_timer, "visualization_files_creation");

  ParticleVisualization<dim> particle_visualization(
    triangulation,
    reconstruction_parameters.reconstruction_positions_file,
    reconstruction_positions,
    reconstruction_counts);

  particle_visualization.build_visualization_files();
}


template class RPTCellReconstruction<3>;
