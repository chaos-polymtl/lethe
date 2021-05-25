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

#include <rpt/parameters_rpt.h>
#include <rpt/radioactive_particle.h>
#include <rpt/rpt.h>
#include <rpt/rpt_calculating_parameters.h>

#include <fstream>
#include <iostream>
#include <iterator>

template <int dim>
RPT<dim>::RPT(RPTCalculatingParameters<dim> &RPTparameters)
  : rpt_parameters(RPTparameters)
{}

template <int dim>
void
RPT<dim>::calculate()
{
  /*
  Point<dim> pt_1;
  pt_1[0] = 1.;
  pt_1[1] = 0.;
  pt_1[2] = 3.;

  Point<dim> pt_2;
  pt_2[0] = 1.;
  pt_2[1] = 0.;
  pt_2[2] = 1.;

  std::cout << "Point 1 is : " << pt_1 << std::endl;
  std::cout << "Point 2 is : " << pt_2 << std::endl;
  std::cout << "Distance between two points : " << pt_1.distance(pt_2)
            << std::endl;
  std::cout << "Vector from 1 to 2 : " << pt_2 - pt_1 << std::endl;


  Tensor<1, dim> v_1({1, 0, 0});
  Tensor<1, dim> v_2({0, 1, 0});

  std::cout << " v_1 : " << v_1 << " v_2 : " << v_2 << std::endl;

  std::cout << " cross-product : " << cross_product_3d(v_1, v_2) << std::endl;
  */

  assign_particle_positions();
  std::cout << particle_positions.size() << std::endl;
  std::cout << particle_positions[0].position[0] << std::endl;

  assign_detector_positions();
  std::cout << detector_positions.size() << std::endl;
  std::cout << detector_positions[0].face_position[0] << std::endl;
}

template <int dim>
void
RPT<dim>::assign_particle_positions()
{
  // Read text file with particle positions and store it in vector
  std::string   line;
  std::ifstream particle_file(rpt_parameters.rpt_param.particle_positions_file);
  std::getline(particle_file, line);
  std::vector<double> values;
  std::copy(std::istream_iterator<double>(particle_file),
            std::istream_iterator<double>(),
            std::back_inserter(values));

  int number_of_positions = (values.size() + 1) / dim;

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

  std::getline(detector_file, line);
  std::vector<double> values;
  std::copy(std::istream_iterator<double>(detector_file),
            std::istream_iterator<double>(),
            std::back_inserter(values));

  // Get the number of detector (2 positions for 1 detector, face and middle)
  int number_of_detector = (values.size() + 1) / (2 * dim);

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
      detector_positions.push_back(detector);
    }
}



template class RPT<3>;
