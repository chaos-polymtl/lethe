#include <core/rheological_model.h>

template <int dim>
std::shared_ptr<RheologicalModel<dim>>
RheologicalModel<dim>::model_cast(
  const Parameters::PhysicalProperties &physical_properties)
{
  if (physical_properties.non_newtonian_parameters.model ==
      Parameters::NonNewtonian::Model::powerlaw)
    return std::make_shared<PowerLaw<dim>>(
      physical_properties.non_newtonian_parameters);
  else if (physical_properties.non_newtonian_parameters.model ==
           Parameters::NonNewtonian::Model::carreau)
    return std::make_shared<Carreau<dim>>(
      physical_properties.non_newtonian_parameters);
  else // if (!physical_properties.non_newtonian_flow)
    return std::make_shared<Newtonian<dim>>(
      physical_properties.fluids[0].viscosity);
}

template class RheologicalModel<2>;
template class RheologicalModel<3>;
template class Newtonian<2>;
template class Newtonian<3>;
template class PowerLaw<2>;
template class PowerLaw<3>;
template class Carreau<2>;
template class Carreau<3>;
