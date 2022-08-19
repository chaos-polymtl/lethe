#include <rpt/rpt_calculating_parameters.h>
#include <rpt/rpt_fem_reconstruction.h>

#include <fstream>
#include <iostream>

using namespace dealii;

int
main(int argc, char *argv[])
{
  try
    {
      if (argc != 2)
        {
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
          std::exit(1);
        }

      ParameterHandler         prm;
      RPTCalculatingParameters rpt_parameters;
      rpt_parameters.declare(prm);

      // Parsing of the file
      prm.parse_input(argv[1]);
      rpt_parameters.parse(prm);

      RPTFEMReconstruction<3> rpt_reconstruct(
        rpt_parameters.rpt_param,
        rpt_parameters.fem_reconstruction_param,
        rpt_parameters.detector_param);
      if (rpt_parameters.fem_reconstruction_param.l2_project_and_reconstruct)
        rpt_reconstruct.L2_project();
      rpt_reconstruct.rpt_fem_reconstruct();
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
