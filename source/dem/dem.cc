/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#include <dem/dem.h>
#include <dem/dem_properties.h>

template <int dim>
DEMSolver<dim>::DEMSolver(ParametersDEM<dim> dem_parameters):
  mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , parameters(dem_parameters)
  , property_pool(DEM::get_number_properties())
  , mapping(1)
{

}

template <int dim>
void DEMSolver<dim>::solve()
{

}


template class DEMSolver<3>;
