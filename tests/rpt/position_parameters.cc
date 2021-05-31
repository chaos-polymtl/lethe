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
 * @brief
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
  Point<3> P1 = {0, 0, 0};
  Point<3> P2 = {1, 0, 0};
  Point<3> P3 = {0, 0, 1};
  Point<3> P4 = {0, 0, 1.25};

  RadioParticle<3> position0(P1, 0);
  RadioParticle<3> position1(P2, 1);
  RadioParticle<3> position2(P3, 2);
  RadioParticle<3> position3(P4, 2);

  std::vector<RadioParticle<3>> particle_positions;

  particle_positions.push_back(position0);
  particle_positions.push_back(position1);
  particle_positions.push_back(position2);
  particle_positions.push_back(position3);

  // Detector position
  Parameters::DetectorParameters detector_param;
  detector_param.radius = 0.5;
  detector_param.length = 1;
  Point<3>    FP        = {0.5, 0, 1};
  Point<3>    MP        = {1, 0, 1};
  Detector<3> detector(detector_param, 0, FP, MP);

  RPTCalculatingParameters rpt_parameters;

  for (unsigned int i_particle = 0; i_particle < 4; i_particle++)
    {
      ParticleDetectorInteractions<3> particle_detector_interactions(
        particle_positions[i_particle], detector, rpt_parameters);
      double h   = particle_detector_interactions.get_h();
      double rho = particle_detector_interactions.get_rho();
      deallog << "Particle position " << i_particle << " : h = " << h
              << ", rho = " << rho << std::endl;
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
