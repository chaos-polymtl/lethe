// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_postprocessing_probes_h
#define lethe_postprocessing_probes_h

#include <core/parameters.h>
#include <core/simulation_control.h>
#include <core/utilities.h>


template <int dim>
class MultiphysicsInterface;


template <int dim>
class PostprocessingProbes
{
public:
  /**
   * TODO AA
   * @brief
   *
   * @param multiphysics_interface
   * @param triangulation
   * @param simulation_control
   * @param[in] probing_points_parameters
   * @param[in] output_folder
   */
  PostprocessingProbes(
    MultiphysicsInterface<dim> *multiphysics_interface,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
    std::shared_ptr<SimulationControl> simulation_control, // for time
    const Parameters::PostProcessing<dim>::ProbingPoints
                      &probing_points_parameters,
    const std::string &output_folder);

  void
  postprocess_probes(Variable variable);
  // TODO AA AssertThrow not implemented for void fraction assertion

  // write, read checkpoint gather tables


private:
  /**
   * Multiphysics interface through which the mapping and the DoF handler is
   * taken.
   */
  MultiphysicsInterface<dim> *multiphysics_interface;

  /// TODO AA
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;

  ///
  std::shared_ptr<SimulationControl>
    simulation_control; // TODO AA check if all of it is needed

  ///
  const Parameters::PostProcessing<dim>::ProbingPoints
    probing_points_parameters;

  ///
  const std::string output_folder;

  /// Tables of
  std::vector<TableHandler> probe_tables;
};


#endif
