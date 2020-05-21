#ifndef lethe_time_integration_utilities_h
#define lethe_time_integration_utilities_h

#include <core/parameters.h>

inline bool
is_sdirk(Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2);
}

inline bool
is_sdirk2(Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2);
}

inline bool
is_sdirk3(Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3);
}

inline bool
is_sdirk_step1(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1);
}

inline bool
is_sdirk_step2(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2);
}

inline bool
is_sdirk_step3(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method ==
          Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3);
}

inline bool
is_bdf(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method == Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf3);
}

inline bool
time_stepping_method_has_two_stages(
  const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::bdf3 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2);
}

inline bool
time_stepping_method_has_three_stages(
  const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::bdf3 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3);
}

#endif
