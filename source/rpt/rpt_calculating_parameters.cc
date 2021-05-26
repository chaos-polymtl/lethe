#include <rpt/rpt_calculating_parameters.h>

void
RPTCalculatingParameters::declare(ParameterHandler &prm)
{
  Parameters::RPTParameters::declare_parameters(prm);
  Parameters::InitialRPTParameters::declare_parameters(prm);
  Parameters::DetectorParameters::declare_parameters(prm);
}

void
RPTCalculatingParameters::parse(ParameterHandler &prm)
{
  rpt_param.parse_parameters(prm);
  initial_param.parse_parameters(prm);
  detector_param.parse_parameters(prm);
}
