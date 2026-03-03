// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/insertion.h>

#include <cmath>
#include <sstream>

template <int dim, typename PropertiesIndex>
Insertion<dim, PropertiesIndex>::Insertion(
  const std::vector<std::shared_ptr<Distribution>>
    &size_distribution_object_container,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
  : removing_particles_in_region(
      dem_parameters.insertion_info.removing_particles_in_region)
{
  distributions_objects = size_distribution_object_container;

  // Boost signal for load balancing
  this->mark_for_update         = true;
  this->change_to_triangulation = triangulation.signals.any_change.connect(
    [&] { this->mark_for_update = true; });

  for (unsigned int i = 0; i < dim; ++i)
    {
      if (dem_parameters.insertion_info.clear_box_point_1[i] <=
          dem_parameters.insertion_info.clear_box_point_2[i])
        {
          p_min[i] = dem_parameters.insertion_info.clear_box_point_1[i];
          p_max[i] = dem_parameters.insertion_info.clear_box_point_2[i];
        }
      else
        {
          p_min[i] = dem_parameters.insertion_info.clear_box_point_2[i];
          p_max[i] = dem_parameters.insertion_info.clear_box_point_1[i];
        }
    }
}

// Prints the insertion information
template <int dim, typename PropertiesIndex>
void
Insertion<dim, PropertiesIndex>::print_insertion_info(
  const unsigned int       &inserted_this_step,
  const unsigned int       &remained_particles,
  const unsigned int       &particle_type,
  const ConditionalOStream &pcout)
{
  std::stringstream ss;

  ss << inserted_this_step << " particles of type " << particle_type
     << " were inserted, " << remained_particles << " particles of type "
     << particle_type << " remaining";

  announce_string(pcout, ss.str(), '*');
}

// Carries out assigning the properties of inserted particles. The output vector
// is used in insert_global_particles as input argument
template <int dim, typename PropertiesIndex>
void
Insertion<dim, PropertiesIndex>::assign_particle_properties(
  const DEMSolverParameters<dim>   &dem_parameters,
  const unsigned int               &inserted_this_step_this_proc,
  const unsigned int               &current_inserting_particle_type,
  const std::vector<Point<dim>>    &insertion_points,
  std::vector<std::vector<double>> &particle_properties)
{
  // Clearing and resizing particle_properties
  particle_properties.reserve(inserted_this_step_this_proc);

  // Getting properties as local parameters
  // TODO: MAYBE CHANGE THE INPUT TO PHYSICAL PROPERTIES DIRECTLY
  auto physical_properties = dem_parameters.lagrangian_physical_properties;

  distributions_objects[current_inserting_particle_type]
    ->particle_size_sampling(inserted_this_step_this_proc);

  // A loop is defined over the number of particles which are going to be
  // inserted at this step
  for (unsigned int particle_counter = 0;
       particle_counter < inserted_this_step_this_proc;
       ++particle_counter)
    {
      double type = current_inserting_particle_type;
      // We make sure that the diameter is positive
      double diameter =
        std::abs(distributions_objects[current_inserting_particle_type]
                   ->particle_sizes[particle_counter]);
      double density =
        physical_properties.density_particle[current_inserting_particle_type];
      double vel_x   = dem_parameters.insertion_info.initial_vel[0];
      double vel_y   = dem_parameters.insertion_info.initial_vel[1];
      double vel_z   = dem_parameters.insertion_info.initial_vel[2];
      double omega_x = dem_parameters.insertion_info.initial_omega[0];
      double omega_y = dem_parameters.insertion_info.initial_omega[1];
      double omega_z = dem_parameters.insertion_info.initial_omega[2];
      double mass    = density * 4. / 3. * M_PI *
                    Utilities::fixed_power<3, double>(diameter * 0.5);

      std::vector<double> properties_of_one_particle{
        type, diameter, mass, vel_x, vel_y, vel_z, omega_x, omega_y, omega_z};

      if constexpr (std::is_same_v<PropertiesIndex,
                                   DEM::DEMMPProperties::PropertiesIndex>)
        {
          double T =
            dem_parameters.insertion_info.initial_temperature_function->value(
              insertion_points[particle_counter]);
          double specific_heat =
            physical_properties
              .specific_heat_particle[current_inserting_particle_type];
          properties_of_one_particle.push_back(T);
          properties_of_one_particle.push_back(specific_heat);
        }

      particle_properties.push_back(properties_of_one_particle);
      properties_of_one_particle.clear();
    }
}

template <int dim, typename PropertiesIndex>
void
Insertion<dim, PropertiesIndex>::
  calculate_insertion_domain_maximum_particle_number(
    const DEMSolverParameters<dim> &dem_parameters,
    const ConditionalOStream       &pcout)
{
  // Getting properties as local parameters
  const auto insertion_information = dem_parameters.insertion_info;

  // Checking if the insertion directions are valid (no repetition)
  int axis_sum = 0;
  if constexpr (dim == 2)
    {
      axis_sum = insertion_information.direction_sequence[0] +
                 insertion_information.direction_sequence[1];

      AssertThrow(
        axis_sum == 1,
        ExcMessage("Invalid insertion directions: 2 directions are the same "));
    }
  if constexpr (dim == 3)
    {
      axis_sum = insertion_information.direction_sequence[0] +
                 insertion_information.direction_sequence[1] +
                 insertion_information.direction_sequence[2];

      AssertThrow(
        axis_sum == 3,
        ExcMessage(
          "Invalid insertion directions: at least 2 directions are the same "));
    }


  // This variable is used for calculation of the maximum number of particles
  // that can fit in the chosen insertion box
  int maximum_particle_number = 1;

  number_of_particles_directions.resize(dim);
  axis_min.resize(dim);
  axis_max.resize(dim);

  std::vector<unsigned int> axis_list = {
    insertion_information.direction_sequence[0],
    insertion_information.direction_sequence[1]};

  if constexpr (dim == 3)
    {
      axis_list.push_back(insertion_information.direction_sequence[2]);
    }

  // Assigning the minimum and maximum positions of the insertion box in respect
  // to the axis order
  for (unsigned int axis : axis_list)
    {
      switch (axis)
        {
          case 0:
            axis_min[0] = insertion_information.insertion_box_point_1(0);
            axis_max[0] = insertion_information.insertion_box_point_2(0);
            break;
          case 1:
            axis_min[1] = insertion_information.insertion_box_point_1(1);
            axis_max[1] = insertion_information.insertion_box_point_2(1);
            break;
          case 2:
            axis_min[2] = insertion_information.insertion_box_point_1(2);
            axis_max[2] = insertion_information.insertion_box_point_2(2);
            break;
          default:
            AssertThrow(false,
                        ExcMessage("Insertion direction must be 0, 1 or 2"));
        }

      // Assign max number of particles according to the direction and calculate
      // the total max number (maximum_particle_number = max_x * max_y * max_z)
      int number_of_particles = static_cast<int>(
        (axis_max[axis] - axis_min[axis]) /
        (insertion_information.distance_threshold * this->maximum_diameter));
      number_of_particles_directions[axis] = number_of_particles;
      maximum_particle_number *= number_of_particles;
    }

  // If the inserted number of particles at this step exceeds the maximum
  // number, a warning is printed
  if (insertion_information.inserted_this_step > maximum_particle_number)
    {
      pcout << "Warning: the requested number of particles for insertion ("
            << insertion_information.inserted_this_step
            << ") is higher than maximum expected number of particles ("
            << maximum_particle_number << ")" << std::endl;

      // Updating the number of inserted particles at each step
      inserted_this_step = maximum_particle_number;
    }
  else
    {
      inserted_this_step = insertion_information.inserted_this_step;
    }
}


template <int dim, typename PropertiesIndex>
void
Insertion<dim, PropertiesIndex>::find_cells_in_removing_box(
  const parallel::distributed::Triangulation<dim> &triangulation)
{
  // Clearing the containers
  in_removal_box.clear();
  edge_of_removal_box.clear();
  bool partially_inside, completely_inside;
  // Looping through cells
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          partially_inside  = false;
          completely_inside = true;
          // Loop through its vertices
          for (unsigned int vertex_id = 0; vertex_id < cell->n_vertices();
               ++vertex_id)
            {
              Point<3> vertex = point_nd_to_3d(cell->vertex(vertex_id));

              // Check if the n_th vertex is in the box
              if (p_min[0] <= vertex[0] && vertex[0] <= p_max[0] &&
                  p_min[1] <= vertex[1] && vertex[1] <= p_max[1] &&
                  p_min[2] <= vertex[2] && vertex[2] <= p_max[2])
                {
                  partially_inside = true;
                }
              else
                {
                  completely_inside = false;
                }
            }
          // If the cell is completely inside the clearing box, both bool will
          // be true. We need to add the cell iterator to the
          // in_removal_box container.

          // If the cell is partially inside, at least one of the "if" will have
          // failed, thus the second bool will be false but the first one true.
          // We need to add the cell iterator to the edge_of_removal_box
          // container.

          // If the cell is outside, both bool will be false, and we do nothing.
          if (completely_inside)
            in_removal_box.insert(cell);
          else if (!completely_inside && partially_inside)
            edge_of_removal_box.insert(cell);
        }
    }
}

template <int dim, typename PropertiesIndex>
void
Insertion<dim, PropertiesIndex>::remove_particles_in_box(
  Particles::ParticleHandler<dim> &particle_handler)
{
  std::vector<
    typename dealii::Particles::ParticleHandler<dim>::particle_iterator>
    to_remove_iterators;

  // Reserve to the maximum number of particle on this processor.
  to_remove_iterators.reserve(particle_handler.n_locally_owned_particles());

  // Loop over the first container
  for (const auto &cell_in_box : in_removal_box)
    {
      // Check if this cell has particles
      auto particles_in_cell = particle_handler.particles_in_cell(cell_in_box);

      if (!particles_in_cell.empty())
        {
          // Loop over the particle in the cell
          for (auto particle_in_cell = particles_in_cell.begin();
               particle_in_cell != particles_in_cell.end();
               ++particle_in_cell)
            {
              // Since we know the cell is fully inside the box, we can
              // remove every particles in it.

              to_remove_iterators.push_back(particle_in_cell);
            }
        }
    }

  // Loop over the second container
  for (auto cell_edge_of_removal_box = edge_of_removal_box.begin();
       cell_edge_of_removal_box != edge_of_removal_box.end();
       ++cell_edge_of_removal_box)
    {
      // Check if this cell has particle
      auto particles_in_cell =
        particle_handler.particles_in_cell(*cell_edge_of_removal_box);
      const bool particles_exist_in_cell = !particles_in_cell.empty();

      if (particles_exist_in_cell)
        {
          // Loop over the particles in the cell
          for (auto particle_in_cell = particles_in_cell.begin();
               particle_in_cell != particles_in_cell.end();
               ++particle_in_cell)
            {
              // We need to check if the particle is inside since the cell
              // is at the edge.
              Point<3> particle_position =
                point_nd_to_3d(particle_in_cell->get_location());
              if (p_min[0] <= particle_position[0] &&
                  p_max[0] >= particle_position[0] &&
                  p_min[1] <= particle_position[1] &&
                  p_max[1] >= particle_position[1] &&
                  p_min[2] <= particle_position[2] &&
                  p_max[2] >= particle_position[2])
                {
                  to_remove_iterators.push_back(particle_in_cell);
                }
            }
        }
    }
  particle_handler.remove_particles(to_remove_iterators);
}

template class Insertion<2, DEM::DEMProperties::PropertiesIndex>;
template class Insertion<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class Insertion<2, DEM::DEMMPProperties::PropertiesIndex>;
template class Insertion<3, DEM::DEMProperties::PropertiesIndex>;
template class Insertion<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class Insertion<3, DEM::DEMMPProperties::PropertiesIndex>;
