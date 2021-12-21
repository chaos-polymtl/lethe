#include <core/rheological_model.h>

template <int dim>
void
RheologicalModel<dim>::model_cast(const Parameters::PhysicalProperties & physical_properties)
{
  if (!physical_properties.non_newtonian_flow)
    *this =
      std::make_unique<Newtonian<dim>>(physical_properties.fluids[0].viscosity);
  else if (physical_properties.non_newtonian_parameters.model ==
          Parameters::NonNewtonian::Model::carreau)
    *this = std::make_unique<Carreau<dim>>(
      physical_properties.non_newtonian_parameters);
}
