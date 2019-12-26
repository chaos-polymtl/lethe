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

#include <tuple>
#include <vector>

using namespace dealii;

#ifndef VISUALIZATION_H_
#  define VISUALIZATION_H_

template <int dim, int spacedim>
class Visualization : public dealii::DataOutInterface<0, dim>
{
public:
  Visualization<dim,spacedim>();

  void
  build_patches(const Particles::ParticleHandler<dim,spacedim> &,
                const unsigned int,
                const unsigned int,
                std::vector<std::tuple<std::string, int>>);

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
};
#endif
