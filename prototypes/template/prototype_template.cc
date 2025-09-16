// SPDX-FileCopyrightText: Copyright (c) 2019 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// This is a template folder for prototype executables to be created when
// developing totally new Lethe functionnalities

#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>

#include <iostream>


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
