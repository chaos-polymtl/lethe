#include <solvers/flow_control.h>

#include <fstream>

template <int dim>
FlowControl<dim>::FlowControl(
  const Parameters::DynamicFlowControl &flow_control)
  : beta_0(flow_control.beta_0)
  , beta_n1(0)
  , flow_rate_0(flow_control.flow_rate_0)
  , flow_rate_n(0)
  , flow_direction(flow_control.flow_direction)
  , no_force(true)
  , threshold_factor(1.01)
{}

template <int dim>
void
FlowControl<dim>::calculate_beta(const std::pair<double, double> &flow_rate,
                                 const double &                   dt,
                                 const unsigned int &             step_number)
{
  // Getting flow rate and area of the last time step.
  flow_rate_n = flow_rate.first;
  area        = flow_rate.second;

  // (Only after step time 1)
  // If flow is now reached with no force, the "no_force" variable is set to
  // false. It may means that slowing down the flow with the pressure drop
  // (no force applied) is ineffective (too big pressure drop) or the flow rate
  // now reached the intended value without force.
  // In case the flow rate is in the initial threshold but the pressure drop is
  // too small to slow it down, the threshold factor decreases by 2.
  // Otherwise, the flow rate may take a lot time to reached the flow rate value
  // only with the pressure drop.
  if (abs(beta_n - 0) < 1e-6 && step_number > 1)
    {
      if (abs(flow_rate_n) < abs(flow_rate_0))
        no_force = false;
      else if (abs(flow_rate_1n - flow_rate_n) <
                 abs(flow_rate_n - flow_rate_0) &&
               threshold_factor > 1.00125)
        threshold_factor = 1 + 0.5 * (threshold_factor - 1);
    }

  // (Only at step time 3)
  // If flow rate is over the desired value at time step 1 and beta applied at
  // time step 2 decreased it under the value, "no_force" is disable.
  // As the previous condition, if flow rate is below the value after the small
  // beta applied at time step 2, it means the pressure drop is too big. Then,
  // slowing down the flow rate with the pressure drop is ineffective.
  // This early disabling prevents to set a potential beta to 0 which could
  // ruins a good flow convergence after many step time.
  if (step_number == 3 && abs(flow_rate_n) < abs(flow_rate_0) &&
      abs(flow_rate_1n) > abs(flow_rate_0))
    no_force = false;

  if (step_number == 0)
    {
      // No force applied at first time step.
      // At this moment, calculate_beta is not called at this time step but
      // the condition prevents a beta calculation if the way it's called
      // changes.
      beta_n1 = 0.0;
    }
  else if (step_number == 1)
    {
      // Initial beta
      beta_n1 = beta_0;
    }
  else if (abs(flow_rate_n) > abs(flow_rate_0) &&
           abs(flow_rate_n) < threshold_factor * abs(flow_rate_0) &&
           no_force == true)
    {
      // If the flow rate is between targeted flow rate value and the threshold
      // and if it didn't reached it the value (no_force is enable), it
      // decreases by itself (pressure drop).
      beta_n1 = 0;
    }
  else if (step_number == 2)
    {
      // The calculated beta at time step 2 is small if the initial beta brings
      // the flow rate close to the fixed value but not enough to get in the
      // threshold.
      beta_n1 = 0.5 * (flow_rate_n - flow_rate_0) / (area * dt);
    }
  else
    {
      // Standard flow controller.
      // Calculate the new beta to control the flow.
      beta_n1 =
        0.75 * beta_n - (flow_rate_0 - 2 * flow_rate_n + flow_rate_1n) / (area * dt);

      // If desired flow rate is reached, new beta only maintains the force to
      // keep the flow at the desired value. Is so, if calculated beta is
      // negative it's set to 0 to avoided +/- pressure.
      if (flow_rate_0 * beta_n1 > 0 && no_force == false)
        beta_n1 = 0;
    }


  // Setting beta coefficient to the tensor according to the flow direction.
  if (flow_direction == 0)
    beta[0] = beta_n1; // beta = f_x
  else if (flow_direction == 1)
    beta[1] = beta_n1; // beta = f_y
  else if (flow_direction == 2)
    beta[2] = beta_n1; // beta = f_z

  // Assigning values of this time step as previous values prior next
  // calculation.
  beta_n       = beta_n1;
  flow_rate_1n = flow_rate_n;
}

template <int dim>
void
FlowControl<dim>::save(std::string prefix)
{
  std::string   filename = prefix + ".flowcontrol";
  std::ofstream output(filename.c_str());
  output << "Flow control" << std::endl;
  output << "Previous_beta " << beta_n << std::endl;
  output << "Previous_flow_rate " << flow_rate_1n << std::endl;
  output << "No_force " << no_force << std::endl;
  output << "Threshold_factor " << threshold_factor << std::endl;
}

template <int dim>
void
FlowControl<dim>::read(std::string prefix)
{
  std::string   filename = prefix + ".flowcontrol";
  std::ifstream input(filename.c_str());
  AssertThrow(input, ExcFileNotOpen(filename));

  std::string buffer;
  std::getline(input, buffer);
  input >> buffer >> beta_n;
  input >> buffer >> flow_rate_1n;
  input >> buffer >> no_force;
  input >> buffer >> threshold_factor;
}

template class FlowControl<2>;
template class FlowControl<3>;
