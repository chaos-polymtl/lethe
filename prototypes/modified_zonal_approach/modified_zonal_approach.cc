/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019- by the Lethe authors
 *
 * This file is part of the Lethe library.
 *
 * The LEthe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019
 */


// This is a template folder for prototype executables to be created when
// developing totally new Lethe functionnalities

#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>


using namespace dealii;

int
main()
{
  try
    {
      Tensor<2, 2> mat;
      mat[0][0] = 1.;
      mat[0][1] = 2;
      mat[1][0] = 3;
      mat[1][1] = 4;
      std::cout << "The tensor initiated is : " << mat << std::endl;
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
