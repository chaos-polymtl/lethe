#include <rpt/rpt_calculating_parameters.h>

void
RPTCalculatingParameters::declare(ParameterHandler &prm)
{
  Parameters::RPTParameters::declare_parameters(prm);
  Parameters::RPTTuningParameters::declare_parameters(prm);
  Parameters::DetectorParameters::declare_parameters(prm);
  Parameters::RPTReconstructionParameters::declare_parameters(prm);
  Parameters::RPTFEMReconstructionParameters::declare_parameters(prm);
  Parameters::Mesh::declare_parameters(prm);
}

void
RPTCalculatingParameters::parse(ParameterHandler &prm)
{
  rpt_param.parse_parameters(prm);
  tuning_param.parse_parameters(prm);
  detector_param.parse_parameters(prm);
  reconstruction_param.parse_parameters(prm);
  fem_reconstruction_param.parse_parameters(prm);
  mesh.parse_parameters(prm);
}
