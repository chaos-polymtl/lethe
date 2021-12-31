#include <core/rheological_model.h>

template <int dim>
std::shared_ptr<RheologicalModel<dim>>
RheologicalModel<dim>::model_cast(
  const Parameters::PhysicalProperties &physical_properties)
{
  if (!physical_properties.non_newtonian_flow)
    return std::make_shared<Newtonian<dim>>(
      physical_properties.fluids[0].viscosity);
  else //  (physical_properties.non_newtonian_parameters.model ==
       // Parameters::NonNewtonian::Model::carreau)
    return std::make_shared<Carreau<dim>>(
      physical_properties.non_newtonian_parameters);
}

template class RheologicalModel<2>;
template class RheologicalModel<3>;
template class Newtonian<2>;
template class Newtonian<3>;
template class Carreau<2>;
template class Carreau<3>;
