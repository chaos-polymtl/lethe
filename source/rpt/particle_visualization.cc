#include <core/solutions_output.h>

#include <deal.II/grid/grid_out.h>

#include <deal.II/numerics/data_out_faces.h>

#include <deal.II/particles/data_out.h>

#include <rpt/particle_visualization.h>


template <int dim>
ParticleVisualization<dim>::ParticleVisualization(
  Triangulation<dim>               &background_triangulation,
  std::string                      &filename,
  std::vector<Point<dim>>          &positions,
  std::vector<std::vector<double>> &counts)
  : particle_positions(positions)
  , particle_counts(counts)
  , visualization_filename(filename)
  , mapping(1)
{
  // Initialize particle handler and empty dof handler with reactor grid
  particle_handler.initialize(background_triangulation,
                              mapping,
                              particle_counts[0].size());
  empty_dof_handler.reinit(background_triangulation);
}

template <int dim>
void
ParticleVisualization<dim>::build_visualization_files()
{
  // Export grid at wall of the reactor grid
  DataOutFaces<dim> grid_out_faces;
  grid_out_faces.attach_dof_handler(empty_dof_handler);
  grid_out_faces.build_patches();
  std::string s        = visualization_filename;
  std::string filename = s.substr(0, s.find(".")) + "_grid_";
  write_boundaries_vtu<dim>(
    grid_out_faces, "./", 0, 0, MPI_COMM_WORLD, filename, 0);

  // Clear memory consuming used objects
  grid_out_faces.clear();
  empty_dof_handler.clear();

  // Insert and move the particle positions prior generate vtu/pvtu files
  std::vector<Point<dim>> position(1);
  for (unsigned int i = 0; i < particle_counts.size(); i++)
    {
      // Put one position in a vector (mandatory to insert particle without
      // using triangulation iterator)
      position[0] = particle_positions[i];

      // If first position, insert it in particle handler, otherwise the
      // position is moved
      if (i == 0)
        particle_handler.insert_particles(position);
      else
        particle_handler.set_particle_positions(position, false);

      // Set properties (counts for all detectors) to the particle position
      particle_handler.begin()->set_properties(particle_counts[i]);

      // Generate files
      output_particles(i);
    }
}

template <int dim>
void
ParticleVisualization<dim>::output_particles(unsigned int it)
{
  // Create DataOut for particles
  Particles::DataOut<dim> particle_output;

  // Set detector ids
  unsigned int             n_detectors  = particle_counts[0].size();
  std::vector<std::string> detector_ids = {"detector_counts_0"};
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      n_detectors, DataComponentInterpretation::component_is_scalar);
  if (n_detectors > 1)
    {
      for (unsigned int i = 1; i < n_detectors; i++)
        detector_ids.push_back("detector_counts_" + std::to_string(i));
    }

  particle_output.build_patches(particle_handler,
                                detector_ids,
                                data_component_interpretation);

  std::string s        = visualization_filename;
  std::string filename = s.substr(0, s.find(".")) + "_";
  particle_output.write_vtu_with_pvtu_record("./",
                                             filename,
                                             it,
                                             MPI_COMM_WORLD);
}

template class ParticleVisualization<3>;
