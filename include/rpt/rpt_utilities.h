/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 - by the Lethe authors
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
 */


#ifndef LETHE_RPT_UTILITIES_H
#define LETHE_RPT_UTILITIES_H

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <rpt/detector.h>
#include <rpt/radioactive_particle.h>

#include <fstream>
#include <iostream>
#include <iterator>

using namespace dealii;

/**
 * @brief Create grid of the reactor vessel (cylinder). This function will be refactored to be deprecated and use the regular mesh functionnalities of Lethe.
 */
template <int dim>
void
attach_grid_to_triangulation_temporary(
  Triangulation<dim>              &triangulation,
  const Parameters::RPTParameters &parameters,
  const unsigned int               n_refinement)
{
  // Generate cylinder (needs rotation and shift to get origin at the bottom
  // with z towards top)
  GridGenerator::subdivided_cylinder(triangulation,
                                     5,
                                     parameters.reactor_radius,
                                     parameters.reactor_height / 2.);
  Tensor<1, dim> axis({0, 1, 0});
  GridTools::rotate(axis, M_PI_2, triangulation);
  const Tensor<1, dim> shift_vector({0, 0, parameters.reactor_height / 2.});
  GridTools::shift(shift_vector, triangulation);

  // Add cylindrical manifold
  const Tensor<1, dim>                direction({0, 0, 1});
  const CylindricalManifold<dim, dim> manifold(direction, {0, 0, 1});
  triangulation.set_manifold(0, manifold);

  // Refine all cells
  triangulation.prepare_coarsening_and_refinement();
  triangulation.refine_global(n_refinement);
}


/**
 * @brief Read text file for positions
 *
 *  @param positions_filename Data filename for particle positions
 */
template <int dim>
std::vector<Point<dim>>
read_positions(std::string &positions_filename)
{
  // Read text file with positions and store it in vector
  std::ifstream particle_file(positions_filename);

  std::string skip;
  std::getline(particle_file, skip); // Skip header line
  std::vector<double> values;
  std::copy(std::istream_iterator<double>(particle_file),
            std::istream_iterator<double>(),
            std::back_inserter(values));

  unsigned int number_of_positions = values.size() / dim;

  // Extract positions
  std::vector<Point<dim>> positions;
  for (unsigned int i = 0; i < number_of_positions; i++)
    {
      Point<dim> point(values[dim * i],
                       values[dim * i + 1],
                       values[dim * i + 2]);

      positions.push_back(point);
    }

  return positions;
}

/**
 * @brief Read text file for particle positions and assign them to particle
 * objects
 *
 *  @param particle_positions_filename Data filename for particle positions
 */
template <int dim>
std::vector<RadioParticle<dim>>
assign_particle_positions(std::string &particle_positions_filename)
{
  std::vector<Point<dim>> positions =
    read_positions<dim>(particle_positions_filename);

  // Extract positions, create point objects and radioactive particles
  std::vector<RadioParticle<dim>> particle_positions;
  for (unsigned int i = 0; i < positions.size(); i++)
    {
      RadioParticle<dim> position(positions[i], i);
      particle_positions.push_back(position);
    }

  return particle_positions;
}


/**
 * @brief Read text file for detector positions and assign them to detector
 * objects
 *
 * @param detector_parameters All parameters related to detector
 */

template <int dim>
std::vector<Detector<dim>>
assign_detector_positions(Parameters::DetectorParameters &detector_parameters)
{
  // Read text file with detector positions and store it in vector
  std::ifstream detector_file(detector_parameters.detector_positions_file);

  std::string skip;
  std::getline(detector_file, skip);
  // Skip header line if the header is present
  if (!isalpha(skip[0]))
    {
      detector_file.seekg(0, std::ios::beg);
    }
  std::vector<double> values;
  std::copy(std::istream_iterator<double>(detector_file),
            std::istream_iterator<double>(),
            std::back_inserter(values));

  // Get the number of detector (2 positions for 1 detector, face and middle)
  int number_of_detector = values.size() / (2 * dim);

  // Extract positions, create point objects and detectors
  std::vector<Detector<dim>> detectors;
  for (int i = 0; i < number_of_detector; i++)
    {
      Point<dim> face_point(values[2 * dim * i],
                            values[2 * dim * i + 1],
                            values[2 * dim * i + 2]);
      Point<dim> middle_point(values[2 * dim * i + dim],
                              values[2 * dim * i + dim + 1],
                              values[2 * dim * i + dim + 2]);

      Detector<dim> detector(detector_parameters, i, face_point, middle_point);

      detectors.push_back(detector);
    }
  return detectors;
}

/**
 * @brief Read text file for count data and store it in a vector
 *
 * @param counts_filename Data filename for counts
 */
template <int dim>
std::vector<double>
read_counts(std::string &counts_filename)
{
  // Read text file with counts
  std::ifstream counts_file(counts_filename);

  std::string skip;
  std::getline(counts_file, skip);
  // Skip header line if the header is present
  if (!isalpha(skip[0]))
    {
      counts_file.seekg(0, std::ios::beg);
    }
  std::vector<double> counts;
  std::copy(std::istream_iterator<double>(counts_file),
            std::istream_iterator<double>(),
            std::back_inserter(counts));

  return counts;
}

/**
 * @brief Read text file for count data and sort them to match all counts
 * to a detector
 * @param counts_filename Data filename for counts
 */

template <int dim>
std::vector<std::vector<double>>
read_detectors_counts(std::string &counts_filename, unsigned int n_detectors)
{
  std::vector<double> counts = read_counts<dim>(counts_filename);

  std::vector<std::vector<double>> detectors_counts;
  for (unsigned int i = 0; i < counts.size(); i += n_detectors)
    {
      std::vector<double>::const_iterator first = counts.begin() + i;
      std::vector<double>::const_iterator last =
        counts.begin() + i + n_detectors;

      std::vector<double> particle_counts(first, last);

      detectors_counts.push_back(particle_counts);
    }

  return detectors_counts;
}



#endif // LETHE_RPT_UTILITIES_H
