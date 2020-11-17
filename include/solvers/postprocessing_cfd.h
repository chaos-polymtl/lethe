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
 * Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2020 -
 */

#include <solvers/flow_control.h>

#ifndef lethe_postprocessing_cfd_h

/**
 * @brief calculate_flow_rate. This function calculates the volumetric flow
 * rate at the selected boundary. It actually calculates the flow rate through
 * the summation of the value at each cell surface with the normal vector,
 * the velocity value and the area.
 *
 * @param dof_handler. The argument used to get finite elements
 *
 * @param present_solution. The vector which contains all the values to
 *                          calculate the flow rate
 *
 * @param boundary_id. The inlet boundary
 *
 * @param fem_parameters. The FEM parameters
 *
 * @param mpi_communicator. The mpi communicator
 */
template <int dim, typename VectorType>
std::pair<double, double>
calculate_flow_rate(const DoFHandler<dim> &dof_handler,
                    const VectorType &     present_solution,
                    const unsigned int &   boundary_id,
                    const Parameters::FEM &fem_parameters,
                    const MPI_Comm &       mpi_communicator);



#  define lethe_postprocessing_cfd_h

#endif
