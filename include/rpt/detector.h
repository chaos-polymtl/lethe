/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2021 by the Lethe authors
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
 * Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2021-
 */

#ifndef lethe_detector_h
#define lethe_detector_h

/**
 * Contains all properties of a detector
 */

#include <deal.II/base/point.h>

#include <rpt/parameters_rpt.h>

using namespace dealii;

template <int dim>
class Detector
{
public:
  /**
   * @brief Constructor for the Detector.
   *
   * @param detector_parameters All parameters related to the detector
   *
   * @param n ID number of the detector
   *
   * @param face_point Position on the center of the face of the detector
   *
   * @param middle_point Position of the middle of the detector
   *
   */
  Detector(Parameters::DetectorParameters &detector_parameters,
           int                             n,
           Point<dim> &                    face_point,
           Point<dim> &                    middle_point)
    : radius(detector_parameters.radius)
    , length(detector_parameters.length)
    , id(n)
    , face_position(face_point)
    , middle_position(middle_point)
  {}

  double
  get_radius()
  {
    return radius;
  }

  double
  get_length()
  {
    return length;
  }

  int
  get_id()
  {
    return id;
  }

  Point<dim>
  get_face_position()
  {
    return face_position;
  }

  Point<dim>
  get_middle_position()
  {
    return middle_position;
  }

private:
  double     radius;
  double     length;
  int        id;
  Point<dim> face_position;
  Point<dim> middle_position;
};


#endif // lethe_detector_h
