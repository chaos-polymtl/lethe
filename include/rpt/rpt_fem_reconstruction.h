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

/**
 * This class calculates an L2 projection of the counts on a given
 * triangulation then uses this projection to reconstruct
 * the position of the particles
 */

#ifndef lethe_rpt_fem_reconstruction_h
#define lethe_rpt_fem_reconstruction_h

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <rpt/detector.h>
#include <rpt/parameters_rpt.h>
#include <rpt/particle_detector_interactions.h>
#include <rpt/radioactive_particle.h>
#include <rpt/rpt_calculating_parameters.h>


using namespace dealii;

template <int dim>
class RPTFEMReconstruction
{
public:
  /**
   * @brief Constructor for the RPTFEMReconstructon
   *
   * @param RPTparameters All parameters and positions needed for the count
   * calculation and the reconstruction
   *
   */
  RPTFEMReconstruction(RPTCalculatingParameters &RPTparameters)
    : fe(1)
    , dof_handler(triangulation)
    , rpt_parameters(RPTparameters)
  {}


  void
  L2_project();


private:
  void
  setup_system();
  void
  assemble_system(unsigned detector_no);
  void
  assign_detector_positions();
  void
  solve_linear_system(unsigned detector_no);
  void
  output_results();
  void
  Loop_over_cells();
  double
  Calculate_Jacobian_1();
  double
  Calculate_Jacobian_2();
  double
  Calculate_Jacobian_3();

  double
  Calculate_Jacobian_4();

  double
  Calculate_Jacobian_5();

  double
  Calculate_Jacobian_6();

  double
  Calculate_Jacobian_7();

  double
  Calculate_Jacobian_8();

  double
  Calculate_Jacobian_9();
  void
  solve();
  double
  f1();

  double
  f2();

  double
  f3();

  void
  assemble_jacobian_for_Newton_method();

  void
  assemble_rhs();

  Tensor<2, dim>jacobian_matrix ;
  Tensor<1, dim>rhs_matrix ;
  Tensor<1,dim>unknown;
  Tensor<1,dim>experimental_count;
  std::vector<std::vector<double>> c;



  Triangulation<dim> triangulation;
  FE_Q<dim>          fe;
  DoFHandler<dim>    dof_handler;

  AffineConstraints<double> constraints;
  SparseMatrix<double>      system_matrix;
  SparsityPattern           sparsity_pattern;
  Vector<double>            system_rhs;

  std::vector<Vector<double>> nodal_counts;

  RPTCalculatingParameters  rpt_parameters;
  std::vector<Detector<dim>> detectors;
};



#endif
