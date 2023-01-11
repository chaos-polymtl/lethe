/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2022 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 *
 */

#include <core/dem_properties.h>

#include <dem/data_containers.h>

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_vector.templates.h>

#include <deal.II/particles/particle_handler.h>

#include <vector>

using namespace dealii;

#ifndef lethe_disable_particle_contact_h
#  define lethe_disable_particle_contact_h

// Special template instance for this class.
// Unsigned integer would have been a better choice, but is not working
// with the Vector class (error because of calling of abs() function)
template class LinearAlgebra::distributed::Vector<int>;

template <int dim>
class DisableParticleContact
{
public:
  DisableParticleContact<dim>();

  enum mobility_status
  {
    inactive,
    active,
    mobile,
    empty,
    n_mobility_status
  };

  void
  identify_mobility_status(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const DoFHandler<dim> &                          background_dh,
    const Particles::ParticleHandler<dim> &          particle_handler,
    MPI_Comm                                         mpi_communicator);

  /**
   * Carries out the calculation of the granular temperature in each local cell.
   * These values are summed up during the post-processing steps (imposed by
   * post-processing frequency), and finally divided by the sampling counter to
   * calculate time-averaged granular temperature distribution.
   *
   * @param triangulation Triangulation
   * @param particle_handler Particle handler
   */
  void
  calculate_average_granular_temperature(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const Particles::ParticleHandler<dim> &          particle_handler);


  inline unsigned int
  check_cell_mobility(
    const typename Triangulation<dim>::active_cell_iterator &cell) const
  {
    unsigned int status = mobility_status::mobile;
    auto         cell_iterator_from_mobile_container =
      status_to_cell[status].find(cell);

    if (cell_iterator_from_mobile_container != status_to_cell[status].end())
      {
        return status;
      }

    status = mobility_status::active;
    auto cell_iterator_from_active_container =
      status_to_cell[status].find(cell);

    if (cell_iterator_from_active_container != status_to_cell[status].end())
      {
        return status;
      }

    return mobility_status::inactive;
  }

  void
  print_cell_mobility_info(
    const Particles::ParticleHandler<dim> &particle_handler,
    double                                 time)
  {
    int num   = 0;
    int total = 0;

    if (!status_to_cell[mobility_status::inactive].empty())
      {
        for (auto &cell : status_to_cell[mobility_status::inactive])
          {
            const unsigned int n_particles_in_cell =
              particle_handler.n_particles_in_cell(cell);

            if (n_particles_in_cell > 0)
              {
                num++;
                total++;
              }
          }
      }

    for (auto &status_set : {status_to_cell[mobility_status::active],
                             status_to_cell[mobility_status::mobile],
                             status_to_cell[mobility_status::empty]})
      {
        if (!status_set.empty())
          {
            for (auto &cell : status_set)
              {
                const unsigned int n_particles_in_cell =
                  particle_handler.n_particles_in_cell(cell);

                if (n_particles_in_cell > 0)
                  {
                    total++;
                  }
              }
          }
      }

    unsigned int all_num   = Utilities::MPI::sum(num, MPI_COMM_WORLD);
    unsigned int all_total = Utilities::MPI::sum(total, MPI_COMM_WORLD);

    if (all_total > 0)
      {
        if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
          {
            std::cout << "info;" << time << ";" << all_num << ";" << all_total
                      << std::endl;
          }
      }
  };

  std::vector<unsigned int>
  get_mobility_status_vector(unsigned int n_active_cells)
  {
    std::vector<unsigned int> status(n_active_cells);

    // Loop over all set with different mobility status
    for (unsigned int i_set = 0; i_set < mobility_status::n_mobility_status;
         i_set++)
      {
        for (auto &cell : status_to_cell[i_set])
          {
            unsigned int cell_id = cell->active_cell_index();
            status[cell_id]      = i_set;
          }
      }

    return status;
  };

  std::vector<typename DEM::dem_data_structures<dim>::cell_set>
  get_mobility_status()
  {
    return status_to_cell;
  }

  LinearAlgebra::distributed::Vector<int>
  get_mobility_at_nodes()
  {
    return mobility_at_nodes;
  }

  unsigned int
  get_n_mobility_status()
  {
    return mobility_status::n_mobility_status;
  }

  void
  set_limit_value(const double granular_temperature,
                  const double solid_fraction)
  {
    granular_temperature_limit = granular_temperature;
    solid_fraction_limit       = solid_fraction;
  }

private:
  void
  check_granular_temperature(
    const typename DEM::dem_data_structures<dim>::cells_neighbor_list
      &                                    cells_neighbor_list,
    const Particles::ParticleHandler<dim> &particle_handler);

  std::vector<typename DEM::dem_data_structures<dim>::cell_set> status_to_cell;
  Vector<double>                                                solid_fractions;

  Vector<double>                          granular_temperature_average;
  LinearAlgebra::distributed::Vector<int> mobility_at_nodes;

  double granular_temperature_limit;
  double solid_fraction_limit;
};



#endif // lethe_disable_particle_contact_h
