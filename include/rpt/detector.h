// SPDX-FileCopyrightText: Copyright (c) 2021-2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
           Point<dim>                     &face_point,
           Point<dim>                     &middle_point)
    : radius(detector_parameters.radius)
    , length(detector_parameters.length)
    , id(n)
    , face_position(face_point)
    , middle_position(middle_point)
    , dead_time(detector_parameters.dead_time[n])
    , activity(detector_parameters.activity[n])
    , attenuation_coefficient_reactor(
        detector_parameters.attenuation_coefficient_reactor[n])
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

  double
  get_dead_time()
  {
    return dead_time;
  }

  double
  get_activity()
  {
    return activity;
  }

  double
  get_attenuation_coefficient_reactor()
  {
    return attenuation_coefficient_reactor;
  }



private:
  double     radius;
  double     length;
  int        id;
  Point<dim> face_position;
  Point<dim> middle_position;
  double     dead_time; // Dead time of the detector per accepted pulse
  double     activity;  // Activity of the tracer
  double     attenuation_coefficient_reactor; // Total linear attenuation
  // coefficient of the medium
};


#endif // lethe_detector_h
