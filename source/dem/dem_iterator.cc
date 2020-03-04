/*
 * DEMiterator.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */
#include <dem/dem_iterator.h>

using namespace dealii;

template <int dim>
DEM_iterator<dim>::DEM_iterator()
{}

template <int dim>
void
DEM_iterator<dim>::reinitialize_force(
  Particles::ParticleHandler<dim> &particle_handler)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Getting properties of particle as local variable
      auto particle_properties = particle->get_properties();

      // Reinitializing forces and momentums of particles in the system
      particle_properties[DEM::PropertiesIndex::force_x] = 0;
      particle_properties[DEM::PropertiesIndex::force_y] = 0;
      particle_properties[DEM::PropertiesIndex::force_z] = 0;

      particle_properties[DEM::PropertiesIndex::M_x] = 0;
      particle_properties[DEM::PropertiesIndex::M_y] = 0;
      particle_properties[DEM::PropertiesIndex::M_z] = 0;
    }
}

template <int dim>
void
DEM_iterator<dim>::engine(
  Particles::ParticleHandler<dim> &particle_handler,
  const Triangulation<dim> &       triangulation,
  int &                            DEM_step,
  double &                         DEM_time,
  std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
    &cell_neighbor_list,
  std::vector<std::map<int, pp_contact_info_struct<dim>>>
    &                                           pairs_in_contact_info,
  std::vector<boundary_cells_info_struct<dim>> &boundary_cells_information,
  std::vector<std::map<int, pw_contact_info_struct<dim>>> &pw_pairs_in_contact,
  DEMSolverParameters<dim> &                               dem_parameters,
  Tensor<1, dim> &                                         g,
  std::vector<std::pair<std::string, int>>                 properties,
  Particles::PropertyPool &                                property_pool,
  PPContactForce<dim> *                                    pp_force_object,
  PWContactForce<dim> *                                    pw_force_object,
  Integrator<dim> *                                        integrator_object,
  PPBroadSearch<dim> *pp_broad_search_object,
  PPFineSearch<dim> * pp_fine_search_object,
  PWBroadSearch<dim> *pw_broad_search_object,
  PWFineSearch<dim> * pw_fine_search_object,
  TimerOutput &       computing_timer)
{
  // Defining parameters as a local variable
  auto local_parameter = dem_parameters;

  // Moving walls

  // Insertion
  computing_timer.enter_subsection("insertion");
  if (fmod(DEM_step, local_parameter.insertionInfo.insertion_frequency) == 1)
    {
      if (DEM_step < local_parameter.insertionInfo.insertion_steps_number)
        {
          // put this if inside the insertion class or use a local variable
          // instead of n_global_particles
          if (particle_handler.n_global_particles() <
              local_parameter.simulationControl
                .total_particle_number) // number < total number
            {
              NonUniformInsertion<dim> ins2(local_parameter);
              // UniformInsertion<dim> ins2(local_parameter);

              ins2.insert(particle_handler,
                          triangulation,
                          property_pool,
                          local_parameter);
            }
        }
    }
  computing_timer.leave_subsection();

  // Sort particles in cells
  computing_timer.enter_subsection("sort_particles_in_cells");
  particle_handler.sort_particles_into_subdomains_and_cells();
  computing_timer.leave_subsection();

  // Force reinitilization
  computing_timer.enter_subsection("reinitialize_forces");
  reinitialize_force(particle_handler);
  computing_timer.leave_subsection();

  // PP contact search
  // PP broad search
  computing_timer.enter_subsection("pp_broad_search");
  // if (step % model_parameters_struct.pp_broad_search_frequency == 0) {
  std::vector<std::pair<Particles::ParticleIterator<dim>,
                        Particles::ParticleIterator<dim>>>
    contact_pair_candidates;

  contact_pair_candidates =
    pp_broad_search_object->find_PP_Contact_Pairs(particle_handler,
                                                  cell_neighbor_list);

  // }
  computing_timer.leave_subsection();

  // PP fine search
  computing_timer.enter_subsection("pp_fine_search");
  pp_fine_search_object->pp_Fine_Search(contact_pair_candidates,
                                        pairs_in_contact_info,
                                        local_parameter.simulationControl.dt);
  computing_timer.leave_subsection();

  // PP contact force
  computing_timer.enter_subsection("pp_contact_force");
  pp_force_object->calculate_pp_contact_force(pairs_in_contact_info,
                                              local_parameter);
  computing_timer.leave_subsection();

  // PW contact search
  // PW broad search
  computing_timer.enter_subsection("pw_broad_search");
  std::vector<std::tuple<std::pair<Particles::ParticleIterator<dim>, int>,
                         Tensor<1, dim>,
                         Point<dim>>>
    pwContactList;
  if (DEM_step % local_parameter.model_parmeters.pw_broad_search_frequency == 0)
    {
      pwContactList = pw_broad_search_object->find_PW_Contact_Pairs(
        boundary_cells_information, particle_handler);
    }
  computing_timer.leave_subsection();

  // PW fine search
  computing_timer.enter_subsection("pw_fine_search");
  pw_fine_search_object->pw_Fine_Search(pwContactList,
                                        pw_pairs_in_contact,
                                        local_parameter.simulationControl.dt);
  computing_timer.leave_subsection();

  // PW contact force:
  computing_timer.enter_subsection("pw_contact_force");
  pw_force_object->calculate_pw_contact_force(pw_pairs_in_contact,
                                              local_parameter);
  computing_timer.leave_subsection();

  // Integration
  computing_timer.enter_subsection("integration");
  integrator_object->integrate(particle_handler,
                               g,
                               local_parameter.simulationControl.dt);
  computing_timer.leave_subsection();

  // Visualization
  computing_timer.enter_subsection("visualization");
  if (DEM_step % local_parameter.simulationControl.write_frequency == 1)
    {
      Visualization<dim> visObj;
      visObj.build_patches(particle_handler, properties);
      WriteVTU<dim> writObj;
      writObj.write_master_files(visObj, local_parameter);
      writObj.writeVTUFiles(visObj, DEM_step, DEM_time, local_parameter);
    }
  computing_timer.leave_subsection();

  // Print iteration

  if (fmod(DEM_step, 1000) == 1)
    {
      std::cout << "Step " << DEM_step << std::endl;
      computing_timer.print_summary();
      std::cout
        << "-------------------------------------------------------------"
           "----------"
        << std::endl;
    }

  // Update:
  DEM_step = DEM_step + 1;
  DEM_time = DEM_step * local_parameter.simulationControl.dt;
}

template class DEM_iterator<2>;
template class DEM_iterator<3>;
