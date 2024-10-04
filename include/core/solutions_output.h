// SPDX-FileCopyrightText: Copyright (c) 2020, 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_solutions_output_h
#define lethe_solutions_output_h

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>


// Lethe includes
#include <core/pvd_handler.h>

using namespace dealii;


/**
 * @brief Output the data out to "group_files" vtu files, with a pvtu file and a pvd to store the timing
 * This function outputs the data out to "group_files" vtu file that are
 * orchestrated by a .pvtu file. An additional pvd file takes care of storing
 * the time associated with each .pvtu file.
 *
 * @param pvd_handler a PVDHandler to store the information about the file name and time associated with it
 *
 * @param data_out the DataOut class to which the data has been attached
 *
 * @param folder a string that contains the path where the results are to be saved
 *
 * @param file_prefix a string that stores the name of the file without the iteration number and the extension
 *
 * @param time the time associated with the file
 *
 * @param iter the iteration number associated with the file
 *
 * @param group_files the number of vtu files that will be generated.
 *
 * @param mpi_communicator The mpi communicator
 *
 * @param digits An optional parameter that specifies the amount of digit used to store iteration number in the file name
 */
template <int dim, int spacedim = dim>
void
write_vtu_and_pvd(PVDHandler                            &pvd_handler,
                  const DataOutInterface<dim, spacedim> &data_out,
                  const std::string                     &folder,
                  const std::string                     &file_prefix,
                  const double                           time,
                  const unsigned int                     iter,
                  const unsigned int                     group_files,
                  const MPI_Comm                        &mpi_communicator,
                  const unsigned int                     digits = 5);

/**
 * @brief Output the Data Out Faces to a single vtu file
 * This function outputs the DataOutFaces to a vtu file.
 *
 * @param data_out the DataOutFaces class to which the data has been attached
 *
 * @param time the time associated with the file
 *
 * @param iter the iteration number associated with the file
 *
 * @param mpi_communicator The mpi communicator
 *
 * @param digits An optional parameter that specifies the amount of digit used to store iteration number in the file name
 */
template <int dim>
void
write_boundaries_vtu(const DataOutFaces<dim> &data_out,
                     const std::string       &folder,
                     const double             time,
                     const unsigned int       iter,
                     const MPI_Comm          &mpi_communicator,
                     const std::string &file_prefix = std::string("boundaries"),
                     const unsigned int digits      = 5);
#endif
