#include <dem/print_initial_information.h>

void
print_initial_information(const ConditionalOStream &pcout,
                          const unsigned int &      n_mpi_processes)
{
  pcout
    << "***************************************************************** \n";
  pcout << "Starting simulation with Lethe/DEM on " << n_mpi_processes
        << " processors" << std::endl;
  pcout << "***************************************************************** "
           "\n\n";
}
