// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>

#include <dem/insertion_packed.h>

#include <deal.II/grid/grid_tools.h>

#include <random>

using namespace DEM;

template <int dim, typename PropertiesIndex>
InsertionPacked<dim, PropertiesIndex>::InsertionPacked(
  const std::vector<std::shared_ptr<Distribution>>
    &size_distribution_object_container,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters,
  const double                                     maximum_particle_diameter)
  : Insertion<dim, PropertiesIndex>(size_distribution_object_container,
                                    triangulation,
                                    dem_parameters)
  , particles_of_each_type_remaining(
      dem_parameters.lagrangian_physical_properties.number.at(0))
{
  // Initializing the current inserting particle type
  current_inserting_particle_type = 0;
  this->inserted_this_step        = 0;
  this->maximum_diameter          = maximum_particle_diameter;


  // Getting properties as local parameters
  const auto insertion_information = dem_parameters.insertion_info;

  std::vector<unsigned int> axis_list = {
    insertion_information.direction_sequence[0],
    insertion_information.direction_sequence[1]};

  if constexpr (dim == 3)
    {
      axis_list.push_back(insertion_information.direction_sequence[2]);
    }

  axis_min.resize(3);
  axis_max.resize(3);
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
    }
}

// The main insertion function.
template <int dim, typename PropertiesIndex>
void
InsertionPacked<dim, PropertiesIndex>::insert(
  Particles::ParticleHandler<dim>                 &particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
{
  // We only insert once with this method.
  if (particle_handler.n_global_particles() != 0)
    return;

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

  std::vector<Point<dim>> insertion_points_on_proc;

  // 1. Initialize the Random Number Generator using the configured seed
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  // We add the particle type, otherwise the position of the Nth particle
  // for type 0 and 1 will be the exact same.
  std::mt19937 rng(dem_parameters.insertion_info.seed_for_insertion +
                   this_mpi_process);


  // We insert all the particle for every type all at once.
  for (unsigned int particle_type = 0;
       particle_type <
       dem_parameters.lagrangian_physical_properties.particle_type_number;
       particle_type++)
    {
      this->inserted_this_step =
        dem_parameters.lagrangian_physical_properties.number.at(particle_type);

      if (this->inserted_this_step == 0)
        continue;

      // Distributing particles between processors
      this->inserted_this_step_this_proc = static_cast<unsigned int>(
        floor(this->inserted_this_step / n_mpi_process));
      if (this_mpi_process == (n_mpi_process - 1))
        this->inserted_this_step_this_proc = static_cast<unsigned int>(
          this->inserted_this_step -
          (n_mpi_process - 1) *
            floor(this->inserted_this_step / n_mpi_process));

      // Clear and reserve space for insertion points for this particle type
      insertion_points_on_proc.clear();
      insertion_points_on_proc.reserve(this->inserted_this_step_this_proc);

      Point<dim> insertion_location;

      // Find the first and the last particle id for each process
      unsigned int first_id;
      unsigned int last_id;
      if (this_mpi_process == (n_mpi_process - 1))
        {
          first_id =
            this->inserted_this_step - this->inserted_this_step_this_proc;
          last_id = this->inserted_this_step;
        }
      else
        {
          first_id = this_mpi_process * this->inserted_this_step_this_proc;
          last_id = (this_mpi_process + 1) * this->inserted_this_step_this_proc;
        }

      // 2. Loop and generate valid insertion points using rejection sampling
      unsigned int particle_counter = 0;
      for (unsigned int id = first_id; id < last_id; ++id, ++particle_counter)
        {
          // Call updated rejection sampling function passing the RNG state
          find_insertion_location(insertion_location,
                                  rng,
                                  dem_parameters.insertion_info,
                                  distribution);

          insertion_points_on_proc.push_back(insertion_location);
        }

      // Assigning inserted particles properties
      std::vector<std::vector<double>> particle_properties;
      this->assign_particle_properties(dem_parameters,
                                       this->inserted_this_step_this_proc,
                                       particle_type,
                                       insertion_points_on_proc,
                                       particle_properties);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(insertion_points_on_proc,
                                               global_bounding_boxes,
                                               particle_properties);

      this->print_insertion_info(this->inserted_this_step,
                                 0,
                                 particle_type,
                                 pcout);
    }
}


template <int dim, typename PropertiesIndex>
void
InsertionPacked<dim, PropertiesIndex>::find_insertion_location(
  Point<dim>                                       &insertion_location,
  std::mt19937                                     &rng,
  const Parameters::Lagrangian::InsertionInfo<dim> &insertion_information,
  std::uniform_real_distribution<double>           &distribution)
{
  // Generate a new candidate point within the bounding box


  for (int d = 0; d < dim; ++d)
    {
      // Get the actual coordinate axis (e.g., 0 for X, 1 for Y, 2 for Z)
      // based on the configured insertion sequence
      const unsigned int axis = insertion_information.direction_sequence.at(d);

      // Pull a fresh 0.0 to 1.0 random value
      double random_value = distribution(rng);

      // Map it uniformly across the specific axis span
      insertion_location[axis] =
        this->axis_min[axis] +
        random_value * (this->axis_max[axis] - this->axis_min[axis]);
    }
}


template class InsertionPacked<2, DEM::DEMProperties::PropertiesIndex>;
template class InsertionPacked<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionPacked<2, DEM::DEMMPProperties::PropertiesIndex>;
template class InsertionPacked<3, DEM::DEMProperties::PropertiesIndex>;
template class InsertionPacked<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionPacked<3, DEM::DEMMPProperties::PropertiesIndex>;
