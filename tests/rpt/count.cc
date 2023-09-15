/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 -  by the Lethe authors
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
 * Author: Audrey Collard-Daigneault Polytechnique Montreal, 2021
 */

/**
 * @brief Assure that the count is always the same with a fixed seed
 * number for random number generator and iteration number.
 */

#include <deal.II/base/point.h>

#include <../tests/tests.h>
#include <rpt/detector.h>
#include <rpt/parameters_rpt.h>
#include <rpt/particle_detector_interactions.h>
#include <rpt/radioactive_particle.h>
#include <rpt/rpt_calculating_parameters.h>

#include <vector>

void
test()
{
  // Particle positions
  Point<3> P0 = {0, 0, 0};       // rho > r | h > 0
  Point<3> P1 = {0.5, 0, 0};     // rho > r | h = 0
  Point<3> P2 = {0, 0, 1};       // rho < r | aligned with detector
  Point<3> P3 = {0.25, 0, 1.25}; // rho < r

  RadioParticle<3> position0(P0, 0);
  RadioParticle<3> position1(P1, 1);
  RadioParticle<3> position2(P2, 2);
  RadioParticle<3> position3(P3, 2);

  const unsigned int n_particle = 4;
  RadioParticle<3>   particle_positions[n_particle]{position0,
                                                    position1,
                                                    position2,
                                                    position3};

  // Detector positions
  Parameters::DetectorParameters detector_param;
  detector_param.radius = 0.5;
  detector_param.length = 1;
  detector_param.dead_time.push_back(1);
  detector_param.activity.push_back(1);
  detector_param.attenuation_coefficient_reactor.push_back(1);
  Point<3>    FP = {0.5, 0, 1};
  Point<3>    MP = {1, 0, 1};
  Detector<3> detector(detector_param, 0, FP, MP);

  // Other parameters
  Parameters::RPTParameters rpt_param;
  rpt_param.reactor_radius                   = 0.5;
  rpt_param.peak_to_total_ratio              = 1;
  rpt_param.sampling_time                    = 1;
  rpt_param.n_monte_carlo_iteration          = 10000;
  rpt_param.seed                             = 1;
  rpt_param.attenuation_coefficient_detector = 1;
  rpt_param.gamma_rays_emitted               = 1;

  // Counts for every particle positions with Monte Carlo
  for (unsigned int i_particle = 0; i_particle < n_particle; i_particle++)
    {
      ParticleDetectorInteractions<3> particle_detector_interactions(
        particle_positions[i_particle], detector, rpt_param);

      double count = particle_detector_interactions.calculate_count();

      deallog << "Particle position " << i_particle << " : count = " << count
              << std::endl;
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      test();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
