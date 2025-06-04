// SPDX-FileCopyrightText: Copyright (c) 2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/signed_distance_transformation.h>

std::shared_ptr<SignedDistanceTransformationBase>
SignedDistanceTransformationBase::model_cast(
  const Parameters::VOF_GeometricInterfaceReinitialization
    &geometric_redistanciation_parameters)
{
  if (geometric_redistanciation_parameters.transformation_type ==
      Parameters::RedistanciationTransformationType::piecewise_polynomial)
    return std::make_shared<SignedDistanceTransformationPiecewisePolynomial>(
      geometric_redistanciation_parameters.max_reinitialization_distance);
  else
    return std::make_shared<SignedDistanceTransformationTanh>(
      geometric_redistanciation_parameters.tanh_thickness);
}
