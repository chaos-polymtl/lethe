/*
 * DEMiterator.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */
#include <dem/dem_iterator.h>

using namespace dealii;

template <int dim> DEM_iterator<dim>::DEM_iterator() {}

template <int dim>
void DEM_iterator<dim>::reinitialize_force(
    Particles::ParticleHandler<dim> &particle_handler) {
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end(); ++particle) {

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
void DEM_iterator<dim>::engine(
    Particles::ParticleHandler<dim> &particle_handler,
    const Triangulation<dim> &triangulation, int &DEM_step, double &DEM_time,
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
        &cell_neighbor_list,
    std::vector<std::map<int, pp_contact_info_struct<dim>>>
        &pairs_in_contact_info,
    std::vector<boundary_cells_info_struct<dim>> &boundary_cells_information,
    std::vector<std::map<int, pw_contact_info_struct<dim>>>
        &pw_pairs_in_contact,
    DEMSolverParameters<dim> &dem_parameters, Tensor<1, dim> &g,
    std::vector<std::pair<std::string, int>> properties,
    Particles::PropertyPool &property_pool,
    PPContactForce<dim> *pp_force_object, PWContactForce<dim> *pw_force_object,
    Integrator<dim> *integrator_object,
    PPBroadSearch<dim> *pp_broad_search_object,
    PPFineSearch<dim> *pp_fine_search_object,
    PWBroadSearch<dim> *pw_broad_search_object,
    PWFineSearch<dim> *pw_fine_search_object) {

  // Defining parameters as a local variable
  auto local_parameter = dem_parameters;

  // moving walls

  auto t1 = std::chrono::high_resolution_clock::now();
  // insertion
  if (fmod(DEM_step, local_parameter.insertionInfo.insertion_frequency) == 1) {
    if (DEM_step < local_parameter.insertionInfo.insertion_steps_number) {
      // put this if inside the insertion class or use a local variable
      // instead of n_global_particles
      if (particle_handler.n_global_particles() <
          local_parameter.simulationControl
              .total_particle_number) // number < total number
      {
        NonUniformInsertion<dim> ins2(dem_parameters);
        // UniformInsertion<dim> ins2(dem_parameters);

        ins2.insert(particle_handler, triangulation, property_pool,
                    dem_parameters);
      }
    }
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration_Insertion =
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

  // sort particles in cells
  particle_handler.sort_particles_into_subdomains_and_cells();

  // force reinitilization
  reinitialize_force(particle_handler);

  // contact search

  auto t3 = std::chrono::high_resolution_clock::now();
  // if (step % model_parameters_struct.pp_broad_search_frequency == 0) {
  std::vector<std::pair<Particles::ParticleIterator<dim>,
                        Particles::ParticleIterator<dim>>>
      contact_pair_candidates;

  contact_pair_candidates = pp_broad_search_object->find_PP_Contact_Pairs(
      particle_handler, cell_neighbor_list);

  // }
  auto t4 = std::chrono::high_resolution_clock::now();
  auto duration_PPContactPairs =
      std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();

  auto t5 = std::chrono::high_resolution_clock::now();
  pp_fine_search_object->pp_Fine_Search(contact_pair_candidates,
                                        pairs_in_contact_info,
                                        local_parameter.simulationControl.dt);

  auto t6 = std::chrono::high_resolution_clock::now();
  auto duration_PPFineSearch =
      std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count();

  auto t7 = std::chrono::high_resolution_clock::now();
  // contact force
  pp_force_object->calculate_pp_contact_force(pairs_in_contact_info,
                                              dem_parameters);
  auto t8 = std::chrono::high_resolution_clock::now();
  auto duration_PPContactForce =
      std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count();

  // p-w contact detection:
  std::vector<std::tuple<std::pair<Particles::ParticleIterator<dim>, int>,
                         Tensor<1, dim>, Point<dim>>>
      pwContactList;

  auto t9 = std::chrono::high_resolution_clock::now();
  if (DEM_step % local_parameter.model_parmeters.pw_broad_search_frequency ==
      0) {
    pwContactList = pw_broad_search_object->find_PW_Contact_Pairs(
        boundary_cells_information, particle_handler);
  }
  auto t10 = std::chrono::high_resolution_clock::now();
  auto duration_PWContactPairs =
      std::chrono::duration_cast<std::chrono::microseconds>(t10 - t9).count();

  auto t11 = std::chrono::high_resolution_clock::now();
  pw_fine_search_object->pw_Fine_Search(pwContactList, pw_pairs_in_contact,
                                        local_parameter.simulationControl.dt);

  auto t12 = std::chrono::high_resolution_clock::now();
  auto duration_PWFineSearch =
      std::chrono::duration_cast<std::chrono::microseconds>(t12 - t11).count();

  auto t13 = std::chrono::high_resolution_clock::now();
  // p-w contact force:
  pw_force_object->calculate_pw_contact_force(pw_pairs_in_contact,
                                              dem_parameters);
  auto t14 = std::chrono::high_resolution_clock::now();
  auto duration_PWContactForce =
      std::chrono::duration_cast<std::chrono::microseconds>(t14 - t13).count();

  auto t15 = std::chrono::high_resolution_clock::now();
  // Integration
  integrator_object->integrate(particle_handler, g,
                               local_parameter.simulationControl.dt);
  auto t16 = std::chrono::high_resolution_clock::now();
  auto duration_Integration =
      std::chrono::duration_cast<std::chrono::microseconds>(t16 - t15).count();
  auto t17 = std::chrono::high_resolution_clock::now();

  // visualization
  if (DEM_step % local_parameter.simulationControl.write_frequency == 1) {
    Visualization<dim> visObj;
    visObj.build_patches(particle_handler, properties);
    WriteVTU<dim> writObj;
    writObj.write_master_files(visObj, dem_parameters);
    writObj.writeVTUFiles(visObj, DEM_step, DEM_time, dem_parameters);
  }
  auto t18 = std::chrono::high_resolution_clock::now();
  auto duration_Visualization =
      std::chrono::duration_cast<std::chrono::microseconds>(t18 - t17).count();

  // print iteration
  if (fmod(DEM_step, 1000) == 1) {
    std::cout << "Step " << DEM_step << std::endl;
    std::cout << "CPU time of insertion is: " << duration_Insertion << " micros"
              << std::endl;
    std::cout << "CPU time of P-P borad search is: " << duration_PPContactPairs
              << " micros" << std::endl;
    std::cout << "CPU time of P-P fine search is: " << duration_PPFineSearch
              << " micros" << std::endl;
    std::cout << "CPU time of P-P contact force is: " << duration_PPContactForce
              << " micros" << std::endl;
    std::cout << "CPU time of P-W borad search is: " << duration_PWContactPairs
              << " micros" << std::endl;
    std::cout << "CPU time of P-W fine search is: " << duration_PWFineSearch
              << " micros" << std::endl;
    std::cout << "CPU time of P-W contact force is: " << duration_PWContactForce
              << " micros" << std::endl;
    std::cout << "CPU time of integration is: " << duration_Integration
              << " micros" << std::endl;
    std::cout << "CPU time of visualization is: " << duration_Visualization
              << " micros" << std::endl;
    std::cout << "-------------------------------------------------------------"
                 "----------"
              << std::endl;
  }

  // update:
  DEM_step = DEM_step + 1;
  DEM_time = DEM_step * local_parameter.simulationControl.dt;
}

template class DEM_iterator<2>;
template class DEM_iterator<3>;
