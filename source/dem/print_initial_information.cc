#include <dem/print_initial_information.h>

void
print_initial_information(const ConditionalOStream &pcout,
                          const unsigned int &      n_mpi_processes)
{
  pcout << "Running on " << n_mpi_processes << " rank(s)" << std::endl;
  pcout << "***************************************************************** "
           "\n\n";
}
