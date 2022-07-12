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

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <rpt/detector.h>
#include <rpt/parameters_rpt.h>
#include <rpt/particle_detector_interactions.h>
#include <rpt/radioactive_particle.h>
#include <rpt/rpt_calculating_parameters.h>
#include <deal.II/base/timer.h>


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
            , mapping(FE_SimplexP<dim>(1))
    {}


    void
    L2_project();

    void
    test();

    void
    rpt_fem_reconstruct();


private:
    void
    setup_triangulation();
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
    output_raw_results_per_level();
    void
    output_counts_on_level(
            unsigned int                                  level,
            std::map<types::global_dof_index, Point<dim>> &dof_index_and_location);
    void
    find_cell(std::vector<double> experimental_count);
    void
    trajectory();
    void
    checkpoint();
    void
    load_from_checkpoint();


    void
    cylinder_shell();

    void
    post_process_L2_projection();


    Triangulation<dim> triangulation;
    FE_SimplexP<dim>   fe;
    DoFHandler<dim>    dof_handler;
    RPTCalculatingParameters    rpt_parameters;
    MappingFE<dim>     mapping;
    unsigned int n_detector;

    AffineConstraints<double>   constraints;
    SparseMatrix<double>        system_matrix;
    SparsityPattern             sparsity_pattern;
    Vector<double>              system_rhs;
    std::vector<Vector<double>> nodal_counts;
    std::vector<Detector<dim>>  detectors;
};

template <int dim>
Vector<double>
assemble_matrix_and_rhs(std::vector<std::vector<double> > &vertex_count,
                        std::vector<double> &experimental_count
);

#endif
