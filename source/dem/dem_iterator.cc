/*
 * DEMiterator.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */
#include "dem/dem_iterator.h"

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/property_pool.h>

#include <math.h>

#include <chrono>

#include "dem/contact_info_struct.h"
#include "dem/dem_solver_parameters.h"
#include "dem/nonuniform_insertion.h"
#include "dem/uniform_insertion.h"
#include "dem/visualization.h"
#include "dem/write_vtu.h"

using namespace dealii;

template <int dim, int spacedim> DEM_iterator<dim, spacedim>::DEM_iterator() {}

template <int dim, int spacedim>
void DEM_iterator<dim, spacedim>::forceReinit(
    Particles::ParticleHandler<dim, spacedim> &particle_handler) {
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end(); ++particle) {
    particle->get_properties()[13] = 0;
    particle->get_properties()[14] = 0;
    particle->get_properties()[15] = 0;

    particle->get_properties()[21] = 0;
    particle->get_properties()[22] = 0;
    particle->get_properties()[23] = 0;
  }
}

/*
void DEM_iterator::checkSimBound(
  Particles::ParticleHandler<3, 3> &particle_handler,
  ReadInputScript                   readInput)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      if (particle->get_properties()[4] < readInput.x_min ||
          particle->get_properties()[4] > readInput.x_max ||
          particle->get_properties()[5] < readInput.y_min ||
          particle->get_properties()[5] > readInput.y_max ||
          particle->get_properties()[6] < readInput.z_min ||
          particle->get_properties()[6] > readInput.z_max)
        {
          particle_handler.remove_particle(particle);
        }
    }
}
*/

template <int dim, int spacedim>
void DEM_iterator<dim, spacedim>::engine(
    Particles::ParticleHandler<dim, spacedim> &particle_handler,
    const Triangulation<dim, spacedim> &tr, int &step, float &time,
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
        &cellNeighbor,
    std::vector<std::map<int, Particles::ParticleIterator<dim, spacedim>>>
        &inContactPairs,
    std::vector<std::map<int, contact_info_struct<dim, spacedim>>>
        &inContactInfo,
    std::vector<boundary_cells_info_struct<dim>> boundary_cells_information,
    std::vector<std::tuple<Particles::ParticleIterator<dim, spacedim>,
                           Tensor<1, dim>, Point<dim>, double, double, double,
                           Point<dim>, double>> &pwContactInfo,
    std::vector<std::tuple<std::string, int>> properties,
    Particles::PropertyPool &property_pool,
    ParticleWallContactDetection<dim, spacedim> pw,
    PPContactForce<dim, spacedim> *pplf,
    ParticleWallContactForce<dim, spacedim> pwcf,
    Integrator<dim, spacedim> *Integ1, double dt, int nTotal, int writeFreq,
    physical_info_struct<dim> physical_info_struct,
    insertion_info_struct<dim, spacedim> insertion_info_struct,
    Tensor<1, dim> g, int numFields, int numProperties,
    PPBroadSearch<dim, spacedim> ppbs, PPFineSearch<dim, spacedim> ppfs,
    PWBroadSearch<dim, spacedim> pwbs) {
  // moving walls

  auto t1 = std::chrono::high_resolution_clock::now();
  // insertion
  if (fmod(step, insertion_info_struct.insertion_frequency) == 1) {
    if (step < insertion_info_struct.insertion_steps_number) {
      // put this if inside the insertion class or use a local variable
      // instead of n_global_particles
      if (particle_handler.n_global_particles() <
          nTotal) // number < total number
      {
        NonUniformInsertion<dim, spacedim> ins2(physical_info_struct,
                                                insertion_info_struct);

        ins2.insert(particle_handler, tr, property_pool, physical_info_struct,
                    insertion_info_struct);
      }
    }
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration_Insertion =
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

  // force reinitilization
  forceReinit(particle_handler);

  // contact search

  auto t3 = std::chrono::high_resolution_clock::now();
  // if (fmod(step, 10) == 1) {
  std::vector<std::pair<Particles::ParticleIterator<dim, spacedim>,
                        Particles::ParticleIterator<dim, spacedim>>>
      contact_pair_candidates;

  contact_pair_candidates =
      ppbs.find_PP_Contact_Pairs(particle_handler, cellNeighbor);
  //}
  auto t4 = std::chrono::high_resolution_clock::now();
  auto duration_PPContactPairs =
      std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();

  auto t5 = std::chrono::high_resolution_clock::now();
  ppfs.pp_Fine_Search(contact_pair_candidates, inContactPairs, inContactInfo,
                      dt, particle_handler);

  auto t6 = std::chrono::high_resolution_clock::now();
  auto duration_PPFineSearch =
      std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count();

  auto t7 = std::chrono::high_resolution_clock::now();
  // contact force
  // cf.nonLinearCF(inContactInfo, physical_info_struct);
  pplf->calculate_pp_contact_force(inContactInfo, physical_info_struct);
  auto t8 = std::chrono::high_resolution_clock::now();
  auto duration_PPContactForce =
      std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count();

  // p-w contact detection:
  std::vector<std::tuple<Particles::ParticleIterator<dim, spacedim>,
                         Tensor<1, dim>, Point<dim>>>
      pwContactList;

  auto t9 = std::chrono::high_resolution_clock::now();
  // if (fmod(step, 10) == 1) {

  pwContactList =
      pwbs.find_PW_Contact_Pairs(boundary_cells_information, particle_handler);
  //}
  auto t10 = std::chrono::high_resolution_clock::now();
  auto duration_PWContactPairs =
      std::chrono::duration_cast<std::chrono::microseconds>(t10 - t9).count();

  auto t11 = std::chrono::high_resolution_clock::now();
  pw.pwFineSearch(pwContactList, pwContactInfo, dt);
  auto t12 = std::chrono::high_resolution_clock::now();
  auto duration_PWFineSearch =
      std::chrono::duration_cast<std::chrono::microseconds>(t12 - t11).count();

  auto t13 = std::chrono::high_resolution_clock::now();
  // p-w contact force:
  pwcf.pwNonLinearCF(pwContactInfo, physical_info_struct);
  // pwcf.pwLinearCF(pwContactInfo, physical_info_struct);
  auto t14 = std::chrono::high_resolution_clock::now();
  auto duration_PWContactForce =
      std::chrono::duration_cast<std::chrono::microseconds>(t14 - t13).count();

  auto t15 = std::chrono::high_resolution_clock::now();
  // Integration
  // Integ1.eulerIntegration(particle_handler, g, dt);
  // Integ1.rk2Integration(particle_handler, g, dt);
  Integ1->integrate(particle_handler, g, dt);
  auto t16 = std::chrono::high_resolution_clock::now();
  auto duration_Integration =
      std::chrono::duration_cast<std::chrono::microseconds>(t16 - t15).count();

  auto t17 = std::chrono::high_resolution_clock::now();
  // visualization
  if (fmod(step, writeFreq) == 1) {
    Visualization<dim, spacedim> visObj;
    visObj.build_patches(particle_handler, numFields, numProperties,
                         properties);
    WriteVTU<dim, spacedim> writObj;
    writObj.write_master_files(visObj);
    writObj.writeVTUFiles(visObj, step, time);
  }
  auto t18 = std::chrono::high_resolution_clock::now();
  auto duration_Visualization =
      std::chrono::duration_cast<std::chrono::microseconds>(t18 - t17).count();

  // print iteration
  if (fmod(step, 1000) == 1) {
    std::cout << "Step " << step << std::endl;
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
  particle_handler.sort_particles_into_subdomains_and_cells();
  step = step + 1;
  time = step * dt;
}

template class DEM_iterator<3, 3>;
