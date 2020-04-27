/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
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
 * Author: Bruno Blais, Polytechnique Montreal, 2020 -
 */

#ifndef lethe_solutions_output_h
#define lethe_solutions_output_h

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>


// Lethe includes
#include <core/pvd_handler.h>

using namespace dealii;


/**
 * @brief Calculate the CFL condition on the simulation domain
 * @return CFL maximal value in the domain
 * Post-processing function
 * This function calculates the maximal CFL value in the domain
 *
 * @param mpi_communicator The mpi communicator. It is used to reduce the CFL calculation.
 */
template <int dim>
void
write_vtu_and_pvd(PVDHandler &        pvd_handler,
                  const DataOut<dim> &data_out,
                  const std::string   folder,
                  const std::string   file_prefix,
                  const double        time,
                  const unsigned int  iter,
                  const unsigned int  group_files,
                  const MPI_Comm &    mpi_communicator,
                  const unsigned int  digits = 4);

template <int dim>
void
write_boundaries_vtu(const DataOutFaces<dim> &data_out,
                     const std::string        folder,
                     const double             time,
                     const unsigned int       iter,
                     const MPI_Comm &         mpi_communicator,
                     const std::string  file_prefix = std::string("boundaries"),
                     const unsigned int digits      = 4);
#endif
