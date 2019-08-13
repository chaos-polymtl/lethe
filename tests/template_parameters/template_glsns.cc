// check the read and write of simulationcontrol

#include "../tests.h"
#include "navierstokessolverparameters.h"

int main() {
  try {
    initlog();

    deallog << "Beggining" << std::endl;
    ParameterHandler prm;
    NavierStokesSolverParameters<2> nsparam;
    nsparam.declare(prm);
    std::ofstream output_prm("template.prm");
    prm.print_parameters(output_prm, prm.Text);

    std::ofstream output_xml("template.xml");
    prm.print_parameters(output_xml, prm.XML);

  } catch (std::exception &exc) {
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
  } catch (...) {
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
