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
* Author: Lucka, Polytechnique Montreal, 2023 -
                           */


#ifndef lethe_mesh_controler_h
#define lethe_mesh_controler_h


using namespace dealii;


template <int dim>
class MeshControler
{
public:
  double
  calculate_corsen_fraction();

};

#endif
