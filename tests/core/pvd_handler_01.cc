// check the read and write of simulationcontrol

#include "../tests.h"
#include "core/pvd_handler.h"

int
main()
{
  try
    {
      initlog();

      deallog << "Beggining" << std::endl;

      PVDHandler pvdhandlerMaster;

      pvdhandlerMaster.append(0.03, "time_at_step_1");
      pvdhandlerMaster.append(2.00, "time_at_step_2");
      pvdhandlerMaster.append(7.00, "time_at_step_3");
      pvdhandlerMaster.save("restart");

      PVDHandler pvdhandlerWorker;
      pvdhandlerWorker.read("restart");

      if (pvdhandlerWorker.size() != pvdhandlerMaster.size())
        throw std::runtime_error("Size are not equal");
      unsigned int size = pvdhandlerWorker.size();
      for (unsigned int i = 0; i < size; ++i)
        {
          deallog << pvdhandlerWorker.times_and_names_[i].first << " "
                  << pvdhandlerWorker.times_and_names_[i].second << std::endl;
          if (!approximatelyEqual(pvdhandlerMaster.times_and_names_[i].first,
                                  pvdhandlerWorker.times_and_names_[i].first,
                                  1e-8))
            throw std::runtime_error("Time not equal");
          if (pvdhandlerMaster.times_and_names_[i].second !=
              pvdhandlerWorker.times_and_names_[i].second)
            throw std::runtime_error("File not equal");
        }
      deallog << "OK" << std::endl;
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
}
