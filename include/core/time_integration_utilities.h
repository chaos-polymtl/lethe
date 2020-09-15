#ifndef lethe_time_integration_utilities_h
#define lethe_time_integration_utilities_h

#include <core/parameters.h>


/**
 * @brief Determines if the time integration method is a steady-state method
 *
 * @param method A time integration method
 */
inline bool
is_steady(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method == Parameters::SimulationControl::TimeSteppingMethod::steady ||
          method ==
            Parameters::SimulationControl::TimeSteppingMethod::steady_bdf);
}
/**
 * @brief Determines if the time integration method is within the sdirk family
 *
 * @param method A time integration method
 */
inline bool
is_sdirk(const Parameters::SimulationControl::TimeSteppingMethod method)
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

/**
 * @brief Determines if the time integration method is within the sdirk22 family
 *
 * @param method A time integration method
 */
inline bool
is_sdirk2(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2);
}

/**
 * @brief Determines if the time integration method is within the sdirk33 family
 *
 * @param method A time integration method
 */
inline bool
is_sdirk3(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3);
}

/**
 * @brief Determines if this is the first step of an sdirk method
 *
 * @param method A time integration method
 */
inline bool
is_sdirk_step1(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1);
}

/**
 * @brief Determines if this is the second step of an sdirk method
 *
 * @param method A time integration method
 */
inline bool
is_sdirk_step2(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2);
}

/**
 * @brief Determines if this is the third step of an sdirk method
 *
 * @param method A time integration method
 */
inline bool
is_sdirk_step3(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method ==
          Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3);
}

/**
 * @brief Determines if the time integration method is within the bdf-1 family
 *
 * @param method A time integration method
 */
inline bool
is_bdf1(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method ==
            Parameters::SimulationControl::TimeSteppingMethod::steady_bdf ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf1);
}

/**
 * @brief Determines if the time integration method is within the bdf family
 *
 * @param method A time integration method
 */
inline bool
is_bdf(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method ==
            Parameters::SimulationControl::TimeSteppingMethod::steady_bdf ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf3);
}

/**
 * @brief Determines if the time integration method is within the high-order (>=2) bdf family
 *
 * @param method A time integration method
 */
inline bool
is_bdf_high_order(
  const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf3);
}

/**
 * @brief Determines if the time integration method requires an additional array
 *
 * @param method A time integration method
 */
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

/**
 * @brief Determines if the time integration method requires two additional arrays
 *
 * @param method A time integration method
 */
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
