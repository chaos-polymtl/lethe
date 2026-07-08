// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>

#include <dem/insertion_packed.h>

#include <deal.II/grid/grid_tools.h>

using namespace DEM;

template <int dim, typename PropertiesIndex>
InsertionPacked<dim, PropertiesIndex>::InsertionPacked(
  const std::vector<std::shared_ptr<Distribution>>
    &size_distribution_object_container,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
  : Insertion<dim, PropertiesIndex>(size_distribution_object_container,
                                    triangulation,
                                    dem_parameters)
  , particles_of_each_type_remaining(
      dem_parameters.lagrangian_physical_properties.number.at(0))
  , acceptance_fct(dem_parameters.insertion_info.insertion_acceptance_fct)
{
  // Initializing the current inserting particle type
  // Getting properties as local parameters
  const InsertionInfo<dim> insertion_information =
    dem_parameters.insertion_info;

  std::vector<unsigned int> axis_list;
  axis_list.reserve(dim);
  axis_list = {insertion_information.direction_sequence[0],
               insertion_information.direction_sequence[1]};
  if constexpr (dim == 3)
    {
      axis_list.push_back(insertion_information.direction_sequence[2]);
    }

  axis_min.resize(3);
  axis_max.resize(3);
  for (unsigned int axis : axis_list)
    {
      AssertThrow(axis < dim,
                  ExcMessage("Insertion direction must be 0, 1 or 2"));

      axis_min[axis] = insertion_information.insertion_box_point_1(axis);
      axis_max[axis] = insertion_information.insertion_box_point_2(axis);
    }
}


template <int dim, typename PropertiesIndex>
void
InsertionPacked<dim, PropertiesIndex>::insert(
  Particles::ParticleHandler<dim>                 &particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
{
  // We only insert once with this method. Thus, if the particle handler already
  // possesses particles, we should not insert again.
  if (particle_handler.n_global_particles() != 0)
    return;

  // Message passing interface requirements.
  MPI_Comm           communicator = triangulation.get_mpi_communicator();
  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  auto this_mpi_process = Utilities::MPI::this_mpi_process(communicator);
  auto n_mpi_process    = Utilities::MPI::n_mpi_processes(communicator);

  // Obtaining global bounding boxes
  const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    triangulation, IteratorFilters::LocallyOwnedCell());
  const auto global_bounding_boxes =
    Utilities::MPI::all_gather(communicator, my_bounding_box);

  // Vector that will get filled with the insertion locations.
  std::vector<Point<dim>> insertion_points_on_proc;

  // Initialize the Random Number Generator using the configured seed
  std::mt19937 rng(dem_parameters.insertion_info.seed_for_insertion +
                   this_mpi_process);

  // We insert all the particle for every type all at once, thus we loop on each
  // particle type.
  for (unsigned int particle_type = 0;
       particle_type <
       dem_parameters.lagrangian_physical_properties.particle_type_number;
       particle_type++)
    {
      const unsigned int n_particle_current_type =
        dem_parameters.lagrangian_physical_properties.number.at(particle_type);

      if (n_particle_current_type == 0)
        continue;

      // Distributing particles between processors
      unsigned int n_particle_current_type_this_proc =
        static_cast<unsigned int>(
          floor(n_particle_current_type / n_mpi_process));
      if (this_mpi_process == (n_mpi_process - 1))
        n_particle_current_type_this_proc = static_cast<unsigned int>(
          n_particle_current_type -
          (n_mpi_process - 1) * floor(n_particle_current_type / n_mpi_process));

      // Clear and reserve space for insertion points for this particle type.
      // We clear, since if this is the second particle type, this vector is
      // already filled.
      insertion_points_on_proc.clear();
      insertion_points_on_proc.reserve(n_particle_current_type_this_proc);

      // Find the first and the last particle id for each process
      unsigned int first_id;
      unsigned int last_id;
      if (this_mpi_process == (n_mpi_process - 1))
        { // Last proc
          first_id =
            n_particle_current_type - n_particle_current_type_this_proc;
          last_id = n_particle_current_type;
        }
      else
        { // Every other proc
          first_id = this_mpi_process * n_particle_current_type_this_proc;
          last_id  = (this_mpi_process + 1) * n_particle_current_type_this_proc;
        }

      // Loop and generate valid insertion points.
      Point<dim>   insertion_location;
      unsigned int particle_counter = 0;
      for (unsigned int id = first_id; id < last_id; ++id, ++particle_counter)
        {
          // Call updated rejection sampling function passing the RNG state
          generate_insertion_location(insertion_location,
                                      rng,
                                      dem_parameters.insertion_info);

          insertion_points_on_proc.push_back(insertion_location);
        }

      // Assigning particle properties
      std::vector<std::vector<double>> particle_properties;
      this->assign_particle_properties(dem_parameters,
                                       n_particle_current_type_this_proc,
                                       particle_type,
                                       insertion_points_on_proc,
                                       particle_properties);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(insertion_points_on_proc,
                                               global_bounding_boxes,
                                               particle_properties);

      this->print_insertion_info(n_particle_current_type,
                                 0,
                                 particle_type,
                                 pcout);
    }
}

template <int dim, typename PropertiesIndex>
void
InsertionPacked<dim, PropertiesIndex>::update_previous_position(
  const parallel::distributed::Triangulation<dim> &triangulation,
  Particles::ParticleHandler<dim>                 &particle_handler)
{
  for (auto cell : triangulation.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell = particle_handler.particles_in_cell(cell);

      if (particles_in_cell.empty())
        continue;

      for (auto &particle : particles_in_cell)
        {
          auto particle_properties = particle.get_properties();

          Point<dim> particle_previous_position = particle.get_location();

          particle_properties[PropertiesIndex::v_x] =
            particle_previous_position[0];
          particle_properties[PropertiesIndex::v_y] =
            particle_previous_position[1];
          if constexpr (dim == 3)
            particle_properties[PropertiesIndex::v_z] =
              particle_previous_position[2];
        }
    }
}

template <int dim, typename PropertiesIndex>
void
InsertionPacked<dim, PropertiesIndex>::clamp_displacement(
  Particles::ParticleHandler<dim> &particle_handler,
  const double                     max_disp,
  std::vector<double>             &displacement)
{
  for (auto &particle : particle_handler)
    {
      auto particle_properties = particle.get_properties();

      Point<dim> previous_position;
      previous_position[0] = particle_properties[PropertiesIndex::v_x];
      previous_position[1] = particle_properties[PropertiesIndex::v_y];
      if constexpr (dim == 3)
        previous_position[2] = particle_properties[PropertiesIndex::v_z];

      const Tensor<1, dim> displacement_tensor =
        particle.get_location() - previous_position;

      const double disp_norm = displacement_tensor.norm();

      // No movement
      if (std::isnan(disp_norm))
        continue;

      const unsigned int particle_id = particle.get_local_index();
      if (disp_norm > max_disp)
        {
          particle.set_location(previous_position +
                                (max_disp / disp_norm) * displacement_tensor);
          displacement[particle_id] += max_disp;
        }
      else
        {
          displacement[particle_id] += disp_norm;
        }
    }
}


template <int dim, typename PropertiesIndex>
void
InsertionPacked<dim, PropertiesIndex>::generate_insertion_location(
  Point<dim>               &insertion_location,
  std::mt19937             &rng,
  const InsertionInfo<dim> &insertion_information)
{
  // A continuous random distribution on the range [min, max) with equal
  // probability throughout the range. Thus the generated points are equally
  // likely to be located in the middle of the box or on its edges.
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  bool accepted = false;
  while (!accepted)
    {
      // Generate a new candidate point within the bounding box
      for (int d = 0; d < dim; ++d)
        {
          // Get the actual coordinate axis (e.g., 0 for X, 1 for Y, 2 for Z)
          // based on the configured insertion sequence
          const unsigned int axis =
            insertion_information.direction_sequence.at(d);

          // Pull a fresh 0.0 to 1.0 random value
          double random_value = distribution(rng);

          // Map it uniformly across the specific axis span
          insertion_location[axis] =
            this->axis_min[axis] +
            random_value * (this->axis_max[axis] - this->axis_min[axis]);
        }
      // Reject the point if the value is negative.
      accepted = acceptance_fct->value(insertion_location) > 0.;
    }
}


template class InsertionPacked<2, DEM::DEMProperties::PropertiesIndex>;
template class InsertionPacked<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionPacked<2, DEM::DEMMPProperties::PropertiesIndex>;
template class InsertionPacked<3, DEM::DEMProperties::PropertiesIndex>;
template class InsertionPacked<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionPacked<3, DEM::DEMMPProperties::PropertiesIndex>;
