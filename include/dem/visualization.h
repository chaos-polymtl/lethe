/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include <deal.II/base/data_out_base.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>

#include <tuple>
#include <vector>

using namespace dealii;

#ifndef visualization_h
#  define visualization_h

/**
 * Building patches of particle properties for visualization
 *
 * @note This function is taken from Aspect and dealii and implemented here
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class Visualization : public dealii::DataOutInterface<0, dim>
{
public:
  Visualization<dim>();

  /**
   * Carries out building the patches of properties of particles for
   * visualization
   *
   * @param particle_handler The particle handler of active particles for
   * visulization
   * @param properties Properties of particles for visulization. This is a
   * vector of pairs and each pair contains the property name as the first
   * element and the size of the property as the second element. For vectors
   * only the size of the first element of the vector is defined equal to the
   * dimension
   * @param g Body force for calculation of the exerted force on particles
   */
  void
  build_patches(Particles::ParticleHandler<dim> &        particle_handler,
                std::vector<std::pair<std::string, int>> properties,
                const Tensor<1, dim> &                   g);

  /**
   * Prints the data of particles in the xyz format
   *
   * @param particle_handler The particle handler of active particles
   * @param g Body force for calculation of the exerted force on particles
   * @param pcout Printing in parallel
   */
  void
  print_xyz(dealii::Particles::ParticleHandler<dim> &particle_handler,
            const Tensor<1, dim> &                   g,
            const ConditionalOStream &               pcout);

  ~Visualization();

private:
  /**
   * Implementation of the corresponding function of the base class.
   */
  virtual const std::vector<DataOutBase::Patch<0, dim>> &
  get_patches() const;

  /**
   * Implementation of the corresponding function of the base class.
   */
  virtual std::vector<std::string>
  get_dataset_names() const;

#  if DEAL_II_VERSION_GTE(9, 1, 0)
  virtual std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
  get_nonscalar_data_ranges() const;
#  else
  virtual std::vector<std::tuple<unsigned int, unsigned int, std::string>>
  get_vector_data_ranges() const;
#  endif

  /**
   * Output information that is filled by build_patches() and
   * written by the write function of the base class.
   */
  std::vector<DataOutBase::Patch<0, dim>> patches;

  /**
   * A list of field names for all data components stored in patches.
   */
  std::vector<std::string> dataset_names;

  /**
   * Store which of the data fields are vectors.
   */
#  if DEAL_II_VERSION_GTE(9, 1, 0)
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vector_datasets;
#  else
  std::vector<std::tuple<unsigned int, unsigned int, std::string>>
    vector_datasets;
#  endif

  /**
   * Number of properties that are written in output files
   */
  const unsigned int number_of_properties_to_write = 16;

  /**
   * Particle properties that are written in output files
   */
  std::vector<std::pair<std::string, int>> properties_to_write;
};
#endif
