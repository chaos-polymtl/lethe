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

#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor.h>

#ifndef lethe_utilities_h
#define lethe_utilities_h

using namespace dealii;

template <int dim>
TableHandler make_table_from_vector_of_tensors(std::vector<Tensor<1, dim>> vector,
                                               std::vector<std::string> column_names,
                                               const unsigned int display_precision);



#endif
