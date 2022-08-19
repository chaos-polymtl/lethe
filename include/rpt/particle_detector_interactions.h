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
 */

#ifndef lethe_particle_detector_interactions_h
#define lethe_particle_detector_interactions_h

/**
 * This class allows to calculate the photon count from a particle received by
 * a detector with the Monte Carlo method.
 */

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <rpt/detector.h>
#include <rpt/parameters_rpt.h>
#include <rpt/radioactive_particle.h>
#include <rpt/rpt_calculating_parameters.h>

template <int dim>
class ParticleDetectorInteractions
{
public:
  /**
   * @brief Constructor for the ParticleDetectorInteractions.
   *
   * @param particle Particle which contains information about its position
   *
   * @param detector Detector which contains information about its positions,
   * radius and length
   *
   * @param rpt_parameters All other parameters needed for the count calculation
   *
   * @param dead_time Dead_time defines the parameter dead-time for the detector
   *
   * @param activity Activity defines the activity of the radioactive source with respect to the detector
   *
   * @param attenuation_coefficient_reactor Attenuation_coefficient_reactor defines the homogeneous attenuation coefficient of the reactor with respect to the reactor
   */
  ParticleDetectorInteractions(RadioParticle<dim> &       particle,
                               Detector<dim> &            detector,
                               Parameters::RPTParameters &rpt_parameters)
    : particle_position(particle.get_position())
    , detector_face_position(detector.get_face_position())
    , detector_middle_position(detector.get_middle_position())
    , detector_radius(detector.get_radius())
    , detector_length(detector.get_length())
    , dead_time(detector.get_dead_time())
    , activity(detector.get_activity())
    , attenuation_coefficient_reactor(
        detector.get_attenuation_coefficient_reactor())
    , parameters(rpt_parameters)
  {}

  /**
   * @brief Calculate photon count of a detector with the Monte Carlo method.
   */
  double
  calculate_count();

  double
  get_h();

  double
  get_rho();

  double
  get_alpha(double n_alpha, double n_theta);

  double
  get_theta(double n_alpha, double n_theta);

  double
  get_detector_path_length(double n_alpha, double n_theta);

  double
  get_reactor_path_length(double n_alpha, double n_theta);

private:
  /**
   * @brief Calculate position parameters (h & rho) of the particle with the
   * detector.
   */
  void
  calculate_position_parameters();

  /**
   * @brief Calculate related angles (theta and alpha) to the solid angle.
   *
   * @param n_alpha The random value [0; 1] generated in the Monte Carlo for
   * alpha
   *
   * @param n_theta The random value [0; 1] generated in the Monte Carlo for
   * theta
   */
  void
  calculate_solid_angle(double n_alpha, double n_theta);

  /**
   * @brief Calculate the length of the photon path through the detector.
   */
  double
  calculate_detector_path_length();

  /**
   * @brief Calculate the length of the photon path through the reactor/tank.
   */
  double
  calculate_reactor_path_length();

  /**
   * @brief Calculate the probability functions of the gamma-rays interaction
   * with the detector.
   *
   * @param detector_path_length Length of the photon path through the detector
   */
  double
  calculate_detector_interaction_probability(double &detector_path_length);

  /**
   * @brief Calculate the probability of non-interaction between the gamma-rays
   * emitted whithin the solid angle and the material inside the reactor/tank
   * its body.
   *
   * @param reactor_path_length Length of the photon path through the reactor
   */
  double
  calculate_non_interaction_probability(double &reactor_path_length);

  /**
   * @brief Calculate the efficiency of the detector with the Monte Carlo method.
   */
  void
  calculate_efficiency();

  /**
   * @brief Solve the t variable of the equation of a straigth line in parametric
   * form. It used the Newton's method with a numerical derivative with initial
   * values of -1 and 1.
   *
   * @param e_inverse Inverse of the detector direction matrix
   *
   * @param detector_particle_origin Origin of the detector-particle coordinate
   *
   * @param particle_position_rotation Transfer of the lab coordinate particle
   * position to the detector-particle coordinate
   */
  std::vector<double> solve_t(Tensor<2, dim> e_inverse,
                              Tensor<1, dim> detector_particle_origin,
                              Tensor<1, dim> particle_position_rotation);

  double         efficiency;
  double         alpha;
  double         theta;
  double         theta_cri;
  double         weighting_factor_alpha;
  double         weighting_factor_theta;
  double         h;
  double         rho;
  double         OA_distance;
  double         OB_distance;
  Tensor<1, dim> detector_orientation_x;
  Tensor<1, dim> detector_orientation_y;
  Tensor<1, dim> detector_orientation_z;

  Point<dim> particle_position;
  Point<dim> detector_face_position;
  Point<dim> detector_middle_position;

  double                    detector_radius;
  double                    detector_length;
  double                    dead_time;
  double                    activity;
  double                    attenuation_coefficient_reactor;
  Parameters::RPTParameters parameters;
};



#endif // lethe_particle_detector_interactions_h
