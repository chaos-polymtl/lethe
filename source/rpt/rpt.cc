/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

*
* Author: Bruno Blais, Ghazaleh Mirakhori Polytechnique Montreal, 2020-
*/

#include <deal.II/base/point.h>

#include <rpt/rpt.h>

#include <fstream>
#include <iostream>
#include <iterator>

template <int dim>
RPT<dim>::RPT(RPTCalculatingParameters &RPTparameters)
  : rpt_parameters(RPTparameters)
{}

template <int dim>
void
RPT<dim>::calculate()
{
  // Reading and assigning positions to particles and detectors
  assign_particle_positions();
  assign_detector_positions();

  // Open a .csv file if exporting results in enable
  std::ofstream myfile;
  if (rpt_parameters.rpt_param.export_counts)
    {
      std::string filename = rpt_parameters.rpt_param.particle_positions_file;
      myfile.open(filename.substr(0, filename.find(".")) + ".csv");
      myfile
        << "particle_positions_x particle_positions_y particle_positions_z detector_id counts"
        << std::endl;
    }

  // Seed the random number generator
  srand(rpt_parameters.rpt_param.seed);

  // Calculate count for every particle-detector pair
  for (unsigned int i_particle = 0; i_particle < particle_positions.size();
       i_particle++)
    {
      for (unsigned int i_detector = 0; i_detector < detectors.size();
           i_detector++)
        {
          // Create the particle-detector interaction object
          ParticleDetectorInteractions<dim> particle_detector_interactions(
            particle_positions[i_particle],
            detectors[i_detector],
            rpt_parameters);

          // Calculate count and print it in terminal
          double count = particle_detector_interactions.calculate_count();
          std::cout << "Count for particle position " << i_particle
                    << " and detector " << i_detector << " : " << count
                    << std::endl;

          // Export results in .csv if enable
          if (myfile.is_open())
            myfile << particle_positions[i_particle].get_position() << " "
                   << detectors[i_detector].get_id() << " " << count
                   << std::endl;
        }
    }

  if (myfile.is_open())
    myfile.close();
}

template <int dim>
void
RPT<dim>::assign_particle_positions()
{
  // Read text file with particle positions and store it in vector
  std::string   line;
  std::ifstream particle_file(rpt_parameters.rpt_param.particle_positions_file);

  std::vector<double> values;
  std::copy(std::istream_iterator<double>(particle_file),
            std::istream_iterator<double>(),
            std::back_inserter(values));

  int number_of_positions = values.size() / dim;

  // Extract positions, create point objects and radioactive particles
  for (int i = 0; i < number_of_positions; i++)
    {
      Point<dim>         point(values[dim * i],
                       values[dim * i + 1],
                       values[dim * i + 2]);
      RadioParticle<dim> position(point, i);
      particle_positions.push_back(position);
    }
}

template <int dim>
void
RPT<dim>::assign_detector_positions()
{
  // Read text file with detector positions and store it in vector
  std::string   line;
  std::ifstream detector_file(
    rpt_parameters.detector_param.detector_positions_file);

  std::vector<double> values;
  std::copy(std::istream_iterator<double>(detector_file),
            std::istream_iterator<double>(),
            std::back_inserter(values));

  // Get the number of detector (2 positions for 1 detector, face and middle)
  int number_of_detector = values.size() / (2 * dim);

  // Extract positions, create point objects and detectors
  for (int i = 0; i < number_of_detector; i++)
    {
      Point<dim> face_point(values[2 * dim * i],
                            values[2 * dim * i + 1],
                            values[2 * dim * i + 2]);
      Point<dim> middle_point(values[2 * dim * i + dim],
                              values[2 * dim * i + dim + 1],
                              values[2 * dim * i + dim + 2]);

      Detector<dim> detector(rpt_parameters.detector_param,
                             i,
                             face_point,
                             middle_point);

      detectors.push_back(detector);
    }
}

template class RPT<3>;
