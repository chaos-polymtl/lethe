// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/effective_medium_model.h>

template <typename ValueType>
std::shared_ptr<EffectiveMediumModel<ValueType>>
EffectiveMediumModel<ValueType>::model_cast(
  const Parameters::CFDDEM &cfd_dem_parameters)
{
  if (cfd_dem_parameters.effective_medium_model_type ==
      Parameters::CFDDEM::EffectiveMediumModelType::none)
    return std::make_shared<NoEffectiveMediumModel<ValueType>>();
  else if (cfd_dem_parameters.effective_medium_model_type ==
           Parameters::CFDDEM::EffectiveMediumModelType::bruggeman)
    return std::make_shared<BruggemanEffectiveMediumModel<ValueType>>();
  else
    AssertThrow(
      false,
      ExcMessage(
        "Invalid effective medium model type. The choices are <none|bruggeman>"));
}

template <>
double
BruggemanEffectiveMediumModel<double>::value(const double void_fraction,
                                             const double fluid_property,
                                             const double solid_property) const
{
  if (fluid_property < 0.0)
    AssertThrow(false, InvalidPhysicalPropertyReal(fluid_property));

  if (solid_property < 0.0)
    AssertThrow(false, InvalidPhysicalPropertyReal(solid_property));

  // We clamp the void fraction to be between 0 and 1, since the Bruggeman model
  // is only valid for void fractions in this range and slightly out of range
  // values can occur from the projected void fraction due to numerical errors.
  // The clamping ensures that the effective property is always computed for
  // valid void fractions.
  double fluid_fraction = std::clamp(void_fraction, 0., 1.);
  double solid_fraction = 1.0 - fluid_fraction;
  double factor         = (3.0 * solid_fraction - 1.0) * solid_property +
                  (3.0 * fluid_fraction - 1.0) * fluid_property;

  double effective_property =
    0.25 * (factor +
            std::sqrt(factor * factor + 8.0 * solid_property * fluid_property));

  return effective_property;
}

template <>
std::complex<double>
BruggemanEffectiveMediumModel<std::complex<double>>::value(
  const double               void_fraction,
  const std::complex<double> fluid_property,
  const std::complex<double> solid_property) const
{
  if (!(std::imag(fluid_property) >= 0.0 && std::real(fluid_property) > 0.0))
    AssertThrow(false, InvalidPhysicalPropertyComplex(fluid_property));

  if (!(std::imag(solid_property) >= 0.0 && std::real(solid_property) > 0.0))
    AssertThrow(false, InvalidPhysicalPropertyComplex(solid_property));

  // We clamp the void fraction to be between 0 and 1, since the Bruggeman model
  // is only valid for void fractions in this range and slightly out of range
  // values can occur from the projected void fraction due to numerical errors.
  // The clamping ensures that the effective property is always computed for
  // valid void fractions.
  double               fluid_fraction = std::clamp(void_fraction, 0., 1.);
  double               solid_fraction = 1.0 - fluid_fraction;
  std::complex<double> factor = (3.0 * solid_fraction - 1.0) * solid_property +
                                (3.0 * fluid_fraction - 1.0) * fluid_property;

  // The std::sqrt function gives one of the branch of the complex plane, to
  // ensure that the effective property is in the upper-half right quadrant of
  // the complex plane, we manually compute both branches and select the proper
  // one.
  std::complex<double> effective_property_branch_1 =
    0.25 * (factor +
            std::sqrt(factor * factor + 8.0 * solid_property * fluid_property));
  std::complex<double> effective_property_branch_2 =
    0.25 * (factor -
            std::sqrt(factor * factor + 8.0 * solid_property * fluid_property));

  // The effective property must be in the upper-half right quadrant of the
  // complex plane, so we select the branch that satisfies this condition.
  if (std::real(effective_property_branch_1) > 0.0 &&
      std::imag(effective_property_branch_1) >= 0.0)
    {
      return effective_property_branch_1;
    }
  else if (std::real(effective_property_branch_2) > 0.0 &&
           std::imag(effective_property_branch_2) >= 0.0)
    {
      return effective_property_branch_2;
    }
  else
    {
      AssertThrow(
        false,
        ExcMessage(
          "The computed effective property is not in the upper-half right quadrant of the complex plane."));
      return std::complex<double>(0.0, 0.0);
    }
}

template <typename ValueType>
void
BruggemanEffectiveMediumModel<ValueType>::vector_value(
  const std::vector<double>    &void_fraction_vector,
  const std::vector<ValueType> &fluid_property_vector,
  const std::vector<ValueType> &solid_property_vector,
  std::vector<ValueType>       &property_vector) const
{
  Assert(void_fraction_vector.size() == fluid_property_vector.size(),
         ExcMessage(
           "The size of the void fraction vector and the fluid property "
           "vector must be the same."));
  Assert(void_fraction_vector.size() == solid_property_vector.size(),
         ExcMessage(
           "The size of the void fraction vector and the solid property "
           "vector must be the same."));
  property_vector.resize(void_fraction_vector.size());

  for (unsigned int i = 0; i < void_fraction_vector.size(); ++i)
    property_vector[i] = this->value(void_fraction_vector[i],
                                     fluid_property_vector[i],
                                     solid_property_vector[i]);
}

template class EffectiveMediumModel<double>;
template class EffectiveMediumModel<std::complex<double>>;
template class NoEffectiveMediumModel<double>;
template class NoEffectiveMediumModel<std::complex<double>>;
template class BruggemanEffectiveMediumModel<double>;
template class BruggemanEffectiveMediumModel<std::complex<double>>;
