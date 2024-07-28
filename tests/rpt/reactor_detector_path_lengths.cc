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
 */

/**
 * @brief Assure path lengths calculation are proceed correctly.
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
  rpt_param.reactor_radius = 0.5;

  // Two sets of n values to check when rho < r => theta < theta_cri & theta >
  // theta_cri
  double n_alpha[3] = {0.1, 0.75, 0.5};
  double n_theta[3] = {0.1, 0.75, 0};

  for (unsigned int i_particle = 0; i_particle < n_particle; i_particle++)
    {
      for (unsigned int i_n = 0; i_n < 3; i_n++)
        {
          ParticleDetectorInteractions<3> particle_detector_interactions(
            particle_positions[i_particle], detector, rpt_param);

          double detector_path_length =
            particle_detector_interactions.get_detector_path_length(
              n_alpha[i_n], n_theta[i_n]);
          double reactor_path_length =
            particle_detector_interactions.get_reactor_path_length(
              n_alpha[i_n], n_theta[i_n]);

          deallog << "Particle position " << i_particle
                  << " : detector path length = " << detector_path_length
                  << std::endl;
          deallog << "                      reactor path length =  "
                  << reactor_path_length << std::endl;
          deallog << "                      n_alpha = " << n_alpha[i_n]
                  << std::endl;
          deallog << "                      n_theta = " << n_theta[i_n]
                  << std::endl;
        }
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
