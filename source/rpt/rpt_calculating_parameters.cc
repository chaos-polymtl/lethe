#include <rpt/rpt_calculating_parameters.h>

template <int dim>
void
RPTCalculatingParameters<dim>::declare(ParameterHandler &prm)
{
  Parameters::RPTParameters::declare_parameters(prm);
  Parameters::InitialRPTParameters::declare_parameters(prm);
  Parameters::DetectorParameters::declare_parameters(prm);
}

template <int dim>
void
RPTCalculatingParameters<dim>::parse(ParameterHandler &prm)
{
  rpt_param.parse_parameters(prm);
  initial_param.parse_parameters(prm);
  detector_param.parse_parameters(prm);
}

template class RPTCalculatingParameters<3>;