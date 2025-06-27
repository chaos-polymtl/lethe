// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_insertion_plane_h
#define lethe_insertion_plane_h

#include <fem-dem/insertion.h>

/**
 * @brief Insertion of particles using cells cut by a plane.
 * Particle insertion using cells cut by a plane. Locally owned cells that are
 * cut by the plane are flagged. From those flagged cells, we insert a particle
 * at their centroids if they are individually empty (contain no particle). This
 * way, no significant overlap is occurring on the insertion of new particle
 * which can occur with other insertion method.
 */
template <int dim, typename PropertiesIndex>
class InsertionPlane : public Insertion<dim, PropertiesIndex>
{
public:
  /**
   * The plane insertion class insert particles using a plane
   * define by a point an a normal vector. This method of insertion can be
   * useful when dealing with a domain close to be fully filled with particle.
   * In this situation, other insertion method have a high risk to create a big
   * overlap between particles on insertion. The plane insertion method mitigate
   * this risk by insertion at the center of empty cells.
   *
   * @param size_distribution_object_container Contains all distribution for each
   * particle type
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  InsertionPlane(const std::vector<std::shared_ptr<Distribution>>
                   &size_distribution_object_container,
                 const parallel::distributed::Triangulation<dim> &triangulation,
                 const DEMSolverParameters<dim> &dem_parameters);

  /**
   * Carries out the insertion of particles.
   * @param particle_handler The particle handler of particles which are being
   * inserted.
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted.
   * @param dem_parameters DEM parameters declared in the .prm file.
   */
  virtual void
  insert(Particles::ParticleHandler<dim>                 &particle_handler,
         const parallel::distributed::Triangulation<dim> &triangulation,
         const DEMSolverParameters<dim> &dem_parameters) override;


  /**
   * @brief Serialize the plane insertion object to an output archive.
   *
   * @param ar Output archive where the attributes are stored.
   */
  virtual void
  serialize(boost::archive::text_oarchive &ar) override
  {
    ar &particles_of_each_type_remaining &current_inserting_particle_type;
  }

  /**
   * @brief Deserialize an input archive to the plane insertion object.
   *
   * @param ar Input archive where the attributes are stored.
   *
   */
  virtual void
  deserialize(boost::archive::text_iarchive &ar) override
  {
    ar &particles_of_each_type_remaining &current_inserting_particle_type;
  }

private:
  /**
   * @brief Store the cells that are cut by the plane.
   *
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param plane_point Point which define the plane
   * @param plane_normal_vector Vector which define the normal direction of
   * the plane.
   */
  void
  find_inplane_cells(
    const parallel::distributed::Triangulation<dim> &triangulation,
    Point<3>                                         plane_point,
    Tensor<1, 3>                                     plane_normal_vector);

  /**
   * @brief Store the location of the centers of all the cells that are in the plane
   */
  void
  find_centers_of_inplane_cells();


  std::set<typename Triangulation<dim>::active_cell_iterator>
                                               plane_cells_for_insertion;
  double                                       maximum_range_for_randomness;
  int                                          particles_of_each_type_remaining;
  unsigned int                                 current_inserting_particle_type;
  std::unordered_map<unsigned int, Point<dim>> cells_centers;
};

#endif
