// SPDX-FileCopyrightText: Copyright (c) 2021, 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/insertion_list.h>

using namespace DEM;

DeclException2(DiameterSizeCoherence,
               int,
               int,
               << "Incoherent particle diameter lists: list x has " << arg1
               << " element(s) and list d has " << arg2 << " element(s).");

/* @brief The constructor of list insertion class does not accomplish anything other
 * than check if the position list are of the coherent size and to
 * create the insertion_points member
 */
template <int dim, typename PropertiesIndex>
InsertionList<dim, PropertiesIndex>::InsertionList(
  const std::vector<std::shared_ptr<Distribution>>
    &size_distribution_object_container,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
  : Insertion<dim, PropertiesIndex>(size_distribution_object_container,
                                    triangulation,
                                    dem_parameters)
  , remaining_particles_of_each_type(
      dem_parameters.lagrangian_physical_properties.number.at(0))
{
  // Initializing current inserting particle type
  current_inserting_particle_type = 0;

  const std::vector<double> &list_x  = dem_parameters.insertion_info.list_x;
  const std::vector<double> &list_y  = dem_parameters.insertion_info.list_y;
  const std::vector<double> &list_z  = dem_parameters.insertion_info.list_z;
  const std::vector<double> &list_vx = dem_parameters.insertion_info.list_vx;
  const std::vector<double> &list_vy = dem_parameters.insertion_info.list_vy;
  const std::vector<double> &list_vz = dem_parameters.insertion_info.list_vz;
  const std::vector<double> &list_wx = dem_parameters.insertion_info.list_wx;
  const std::vector<double> &list_wy = dem_parameters.insertion_info.list_wy;
  const std::vector<double> &list_wz = dem_parameters.insertion_info.list_wz;
  std::vector<double>        list_d  = dem_parameters.insertion_info.list_d;
  std::vector<double>        list_T  = dem_parameters.insertion_info.list_T;

  // If the default diameter is negative, the diameter list is resized to the
  // number of particles and assigned the average particle diameter to all
  // particles
  if (list_d[0] < 0)
    {
      double particle_average_diameter =
        dem_parameters.lagrangian_physical_properties.particle_average_diameter
          .at(current_inserting_particle_type);
      list_d.assign(list_x.size(), particle_average_diameter);
    }

  // Throw an error if the number diameter provided is not coherent with the
  // number of particles.
  AssertThrow(list_x.size() == list_d.size(),
              DiameterSizeCoherence(list_x.size(), list_d.size()));

  // Generate vector of insertion position
  for (unsigned int i = 0; i < list_x.size(); ++i)
    {
      if constexpr (dim == 2)
        {
          insertion_points.emplace_back(Point<dim>({list_x[i], list_y[i]}));
          velocities.emplace_back(Tensor<1, 3>({list_vx[i], list_vy[i], 0.}));
          angular_velocities.emplace_back(Tensor<1, 3>({0., 0., list_wz[i]}));
        }
      if constexpr (dim == 3)
        {
          insertion_points.emplace_back(
            Point<dim>({list_x[i], list_y[i], list_z[i]}));
          velocities.emplace_back(
            Tensor<1, 3>({list_vx[i], list_vy[i], list_vz[i]}));
          angular_velocities.emplace_back(
            Tensor<1, 3>({list_wx[i], list_wy[i], list_wz[i]}));
        }
    }
  diameters    = list_d;
  temperatures = list_T;
}

// The main insertion function. Insert_global_function is used to insert the
// particles
template <int dim, typename PropertiesIndex>
void
InsertionList<dim, PropertiesIndex>::insert(
  Particles::ParticleHandler<dim>                 &particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
{
  // TODO refactor into a function call
  if (remaining_particles_of_each_type == 0 &&
      current_inserting_particle_type !=
        dem_parameters.lagrangian_physical_properties.particle_type_number - 1)
    {
      remaining_particles_of_each_type =
        dem_parameters.lagrangian_physical_properties.number.at(
          ++current_inserting_particle_type);
    }

  if (remaining_particles_of_each_type > 0)
    {
      if (this->removing_particles_in_region)
        {
          if (this->mark_for_update)
            {
              this->find_cells_in_removing_box(triangulation);
              this->mark_for_update = false;
            }
          this->remove_particles_in_box(particle_handler);
        }

      unsigned int n_total_particles_to_insert = insertion_points.size();
      n_total_particles_to_insert =
        std::min(remaining_particles_of_each_type, n_total_particles_to_insert);

      // All processors except 0 will not insert particles
      MPI_Comm communicator = triangulation.get_mpi_communicator();
      auto this_mpi_process = Utilities::MPI::this_mpi_process(communicator);
      const unsigned int n_particles_to_insert_this_proc =
        this_mpi_process == 0 ? n_total_particles_to_insert : 0;

      std::vector<Point<dim>> insertion_points_on_proc_this_step;
      insertion_points_on_proc_this_step.reserve(
        n_particles_to_insert_this_proc);

      // Because the list insertion is made to insert only a few particles
      // only processor 0 manages the insertion of these particles
      if (this_mpi_process == 0)
        {
          for (unsigned int p = 0; p < n_particles_to_insert_this_proc; ++p)
            {
              insertion_points_on_proc_this_step.emplace_back(
                insertion_points[p]);
            }
        }

      // Obtain global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);

      // A vector of vectors, which contains all the properties of all inserted
      // particles at each insertion step
      std::vector<std::vector<double>> particle_properties;

      // Assign inserted particles properties
      this->assign_particle_properties_for_list_insertion(
        dem_parameters,
        n_particles_to_insert_this_proc,
        current_inserting_particle_type,
        particle_properties);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(
        insertion_points_on_proc_this_step,
        global_bounding_boxes,
        particle_properties);

      // Update number of particles remaining to be inserted
      remaining_particles_of_each_type -= n_total_particles_to_insert;


      ConditionalOStream pcout(
        std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);
      this->print_insertion_info(n_total_particles_to_insert,
                                 remaining_particles_of_each_type,
                                 current_inserting_particle_type,
                                 pcout);
    }
}

template <int dim, typename PropertiesIndex>
void
InsertionList<dim, PropertiesIndex>::
  assign_particle_properties_for_list_insertion(
    const DEMSolverParameters<dim>   &dem_parameters,
    const unsigned int               &inserted_this_step_this_proc,
    const unsigned int               &current_inserting_particle_type,
    std::vector<std::vector<double>> &particle_properties)
{
  // Clearing and resizing particle_properties
  particle_properties.reserve(inserted_this_step_this_proc);

  // Getting properties as local parameters
  // TODO: MAYBE CHANGE THE INPUT TO PHYSICAL PROPERTIES DIRECTLY
  auto physical_properties = dem_parameters.lagrangian_physical_properties;

  // A loop is defined over the number of particles which are going to be
  // inserted at this step
  for (unsigned int particle_counter = 0;
       particle_counter < inserted_this_step_this_proc;
       ++particle_counter)
    {
      double type     = current_inserting_particle_type;
      double diameter = this->diameters[particle_counter];
      double density =
        physical_properties.density_particle[current_inserting_particle_type];
      double vel_x   = this->velocities[particle_counter][0];
      double vel_y   = this->velocities[particle_counter][1];
      double vel_z   = this->velocities[particle_counter][2];
      double omega_x = this->angular_velocities[particle_counter][0];
      double omega_y = this->angular_velocities[particle_counter][1];
      double omega_z = this->angular_velocities[particle_counter][2];
      double mass    = density * 4. / 3. * M_PI *
                    Utilities::fixed_power<3, double>(diameter * 0.5);

      std::vector<double> properties_of_one_particle{
        type, diameter, mass, vel_x, vel_y, vel_z, omega_x, omega_y, omega_z};

      if constexpr (std::is_same_v<PropertiesIndex,
                                   DEM::DEMMPProperties::PropertiesIndex>)
        {
          double T = this->temperatures[particle_counter];
          double specific_heat =
            physical_properties
              .specific_heat_particle[current_inserting_particle_type];
          properties_of_one_particle.push_back(T);
          properties_of_one_particle.push_back(specific_heat);
        }

      if constexpr (std::is_same_v<PropertiesIndex,
                                   DEM::CFDDEMProperties::PropertiesIndex>)
        {
          // Push back all zero variables for the CFD-DEM coupling properties
          for (unsigned int i = properties_of_one_particle.size() - 1;
               i < PropertiesIndex::n_properties;
               ++i)
            properties_of_one_particle.push_back(0.);
        }

      particle_properties.push_back(properties_of_one_particle);
      properties_of_one_particle.clear();
    }
}

template class InsertionList<2, DEM::DEMProperties::PropertiesIndex>;
template class InsertionList<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionList<2, DEM::DEMMPProperties::PropertiesIndex>;
template class InsertionList<3, DEM::DEMProperties::PropertiesIndex>;
template class InsertionList<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionList<3, DEM::DEMMPProperties::PropertiesIndex>;
