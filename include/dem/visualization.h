// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_visualization_h
#define lethe_visualization_h

#include <dem/dem_solver_parameters.h>

#include <deal.II/base/data_out_base.h>

#include <deal.II/particles/particle_handler.h>

#include <tuple>
#include <vector>

using namespace dealii;

/**
 * @brief Building patches of particle properties for visualization.
 * This function is taken from Aspect and deal.ii and implemented here.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 */
template <int dim, typename PropertiesIndex>
class Visualization : public dealii::DataOutInterface<0, dim>
{
public:
  Visualization();

  /**
   * @brief Build the patches of properties of particles for visualization.
   *
   * @param particle_handler The particle handler of active particles for
   * visualization.
   * @param properties Properties of particles for visualization. This is a
   * vector of pairs and each pair contains the property name as the first
   * element and the size of the property as the second element. For vectors
   * only the size of the first element of the vector is defined equal to the
   * dimension.
   */
  void
  build_patches(Particles::ParticleHandler<dim> &particle_handler,
                const std::vector<std::pair<std::string, int>> &properties);

  /**
   * @brief Print the data of particles in the xyz format.
   *
   * @param particle_handler The particle handler of active particles.
   * @param pcout Printing in parallel.
   */
  void
  print_xyz(dealii::Particles::ParticleHandler<dim> &particle_handler,
            const MPI_Comm                          &mpi_communicator,
            const ConditionalOStream                &pcout);

  /**
   * @brief Print the data of particles in the deal.II intermediate format.
   *
   * @param data_to_print The vector of data to be printed.
   * @param background_dh The DoFHandler of the background grid.
   * @param mpi_communicator The MPI communicator.
   */
  void
  print_intermediate_format(const Vector<float>   &data_to_print,
                            const DoFHandler<dim> &background_dh,
                            const MPI_Comm        &mpi_communicator);

  ~Visualization();

private:
  /**
   * @brief Implementation of the corresponding function of the base class.
   */
  virtual const std::vector<DataOutBase::Patch<0, dim>> &
  get_patches() const override;

  /**
   * @brief Implementation of the corresponding function of the base class.
   */
  virtual std::vector<std::string>
  get_dataset_names() const override;

  virtual std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
  get_nonscalar_data_ranges() const override;

  /**
   * @brief Output information that is filled by build_patches() and
   * written by the write function of the base class.
   */
  std::vector<DataOutBase::Patch<0, dim>> patches;

  /**
   * @brief A list of field names for all data components stored in patches.
   */
  std::vector<std::string> dataset_names;

  /**
   * @brief Store which of the data fields are vectors.
   */
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vector_datasets;

  /**
   * @brief Particle properties that are written in output files.
   */
  std::vector<std::pair<std::string, int>> properties_to_write;
};
#endif
