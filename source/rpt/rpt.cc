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

*
* Author: Bruno Blais, Ghazaleh Mirakhori Polytechnique Montreal, 2020-
*/

#include "rpt/rpt.h"

#include <iostream>


using namespace dealii;

template <int dim>
void
RPT<dim>::dummy_function()
{
  Point<dim> pt_1;
  pt_1[0] = 1.;
  pt_1[1] = 0.;
  pt_1[2] = 3.;

  Point<dim> pt_2;
  pt_2[0] = 1.;
  pt_2[1] = 0.;
  pt_2[2] = 1.;

  std::cout << "Point 1 is : " << pt_1 << std::endl;
  std::cout << "Point 2 is : " << pt_2 << std::endl;
  std::cout << "Distance between two points : " << pt_1.distance(pt_2)
            << std::endl;
  std::cout << "Vector from 1 to 2 : " << pt_2 - pt_1 << std::endl;


  Tensor<1, dim> v_1({1, 0, 0});
  Tensor<1, dim> v_2({0, 1, 0});

  std::cout << " v_1 : " << v_1 << " v_2 : " << v_2 << std::endl;

  std::cout << " cross-product : " << cross_product_3d(v_1, v_2) << std::endl;
}



template class RPT<3>;
