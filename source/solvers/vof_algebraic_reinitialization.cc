#include <solvers/vof_algebraic_reinitialization.h>

template <int dim>
void
VOFAlgebraicReinitialization<dim>::setup_assemblers()
{
  this->assemblers.clear();

  // Time-stepping schemes
  if (is_bdf(this->simulation_control->get_assembly_method()))
    {
      this->assemblers.push_back(
        std::make_shared<VOFAlgebraicReinitializationAssemblerBDF<dim>>(
          this->simulation_control));
    }

  // Core assembler
  this->assemblers.push_back(
    std::make_shared<VOFAlgebraicReinitializationAssemblerCore<dim>>(
      this->simulation_control,
      this->simulation_parameters.fem_parameters,
      this->simulation_parameters.multiphysics
        .vof_parameters // TODO AMISHGA revoir pour les parameters
      ));
}

template class VOFAlgebraicReinitialization<2>;
template class VOFAlgebraicReinitialization<3>;
