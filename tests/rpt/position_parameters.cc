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
 * @brief Assure position parameters are proceed correctly in any different
 * cases related to the particle position regarding the detector position.
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

  // Detector position
  Parameters::DetectorParameters detector_param;
  detector_param.radius = 0.5;
  detector_param.length = 1;
  detector_param.dead_time.push_back(1);
  detector_param.activity.push_back(1);
  detector_param.attenuation_coefficient_reactor.push_back(1);
  Point<3> FP0 = {0.5, 0, 1};
  Point<3> MP0 = {1, 0, 1};
  Point<3> FP1 = {0.5, 0, 0};
  Point<3> MP1 = {1, 0, 0};

  Detector<3> detector0(detector_param, 0, FP0, MP0);
  Detector<3> detector1(detector_param, 0, FP1, MP1);

  const unsigned int n_detector = 2;
  Detector<3>        detector_positions[n_detector]{detector0, detector1};

  Parameters::RPTParameters rpt_parameters;

  for (unsigned int i_particle = 0; i_particle < n_particle; i_particle++)
    {
      for (unsigned int i_detector = 0; i_detector < n_detector; i_detector++)
        {
          ParticleDetectorInteractions<3> particle_detector_interactions(
            particle_positions[i_particle],
            detector_positions[i_detector],
            rpt_parameters);
          double h   = particle_detector_interactions.get_h();
          double rho = particle_detector_interactions.get_rho();
          deallog << "Particle position " << i_particle << " | Detector "
                  << i_detector << std::endl;
          deallog << " h = " << h << std::endl;
          deallog << " rho = " << rho << std::endl;
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
