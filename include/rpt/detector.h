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
  Detector(Parameters::DetectorParameters detector_parameters,
           int &                          n,
           Point<dim>                     face_point,
           Point<dim>                     middle_point)
    : radius(detector_parameters.radius)
    , length(detector_parameters.length)
    , id(n)
    , face_position(face_point)
    , middle_position(middle_point)
  {}

  double     radius;
  double     length;
  int        id;
  Point<dim> face_position;
  Point<dim> middle_position;
};


#endif // lethe_detector_h
