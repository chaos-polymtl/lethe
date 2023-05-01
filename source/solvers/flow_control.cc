#include <solvers/flow_control.h>

#include <fstream>

template <int dim>
FlowControl<dim>::FlowControl(
  const Parameters::DynamicFlowControl &flow_control)
  : average_velocity_0(flow_control.average_velocity_0)
  , beta_0(flow_control.beta_0)
  , alpha(flow_control.alpha)
  , beta_threshold(flow_control.beta_threshold)
  , beta_n(0.0)
  , flow_direction(flow_control.flow_direction)
  , no_force(true)
  , threshold_factor(1.01)
{}

template <int dim>
void
FlowControl<dim>::calculate_beta(const double &      average_velocity_n,
                                 const double &      dt,
                                 const unsigned int &step_number)
{
  double beta_n1;
  // Check if disabling no_force and modify threshold factor (> step time 1):
  // If average velocity target is reached without any beta force applied,
  // the "no_force" is
  // disabled. It means that the average velocity is now reached the
  // target value without force OR that slowing down the flow with the only
  // pressure drop (no force applied) is ineffective (too big pressure drop).
  // In case the average velocity is in the initial threshold but the
  // pressure drop is too small to slow down the flow, the threshold factor
  // decreases by 2. Otherwise, the average velocity may take a lot time to
  // reached the target value only with pressure drop.
  if (abs(beta_n) < 1e-6 && step_number > 1)
    {
      if (abs(average_velocity_n) < abs(average_velocity_0))
        no_force = false;
      else if (abs(average_velocity_1n - average_velocity_n) <
                 abs(average_velocity_n - average_velocity_0) &&
               threshold_factor > 1.00125)
        threshold_factor = 0.5 * (1 + threshold_factor);
    }

  // Check if disabling no_force is needed (step 3):
  // If average velocity is over the target value at time step 1 and beta
  // force applied at time step 2 decreased it under the value, "no_force" is
  // disabled.
  // If average velocity is below the target value after the small beta
  // applied at time step 2, it means the pressure drop is significative.
  // That means attempting to let the flow slowing down by itself is going to be
  // ineffective. This early disabling prevents setting beta to 0 that could
  // ruins a good flow convergence after many step time.
  if (step_number == 3 && abs(average_velocity_1n) > abs(average_velocity_0) &&
      abs(average_velocity_n) < abs(average_velocity_0))
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
      // Apply initial beta force if user defined (> than 0)
      // otherwise, it is calculated by the proportional controller
      beta_n1 = (abs(beta_0) > 1e-6 ?
                   beta_0 :
                   proportional_flow_controller(average_velocity_n, dt));
    }
  else if (abs(average_velocity_n) > abs(average_velocity_0) &&
           abs(average_velocity_n) <
             threshold_factor * abs(average_velocity_0) &&
           no_force == true)
    {
      // If the average velocity is between targeted value and the threshold
      // (meaning it is sightly over the value) and it has not reached once the
      // value during the simulation (no_force is enable), it decreases by
      // itself (pressure drop).
      beta_n1 = 0.0;
    }
  else if (step_number == 2)
    {
      // Since the main controller take the previous beta value, it is better
      // to use the proportional controller for this second time step since
      // it won't be affected by the first beta value that may cause overshoot
      // or undershoot.
      beta_n1 = proportional_flow_controller(average_velocity_n, dt);
    }
  else
    {
      // Main controller
      beta_n1 = main_flow_controller(average_velocity_n, dt);
    }

  // If the new beta is in the user defined threshold of the old beta, the new
  // beta is kept as the previous beta. This prevents the reassembly of the
  // matrix because of the force term when reuse matrix is used for the
  // non-linear solver
  if (abs(beta_n1 - beta_n) < abs(beta_threshold * beta_n))
    beta_n1 = beta_n;

  // Setting beta coefficient to the tensor according to the flow direction.
  beta[flow_direction] = beta_n1;

  // Assigning values of this time step as previous values prior next
  // calculation.
  beta_n              = beta_n1;
  average_velocity_1n = average_velocity_n;
}

template <int dim>
void
FlowControl<dim>::save(std::string prefix)
{
  std::string   filename = prefix + ".flowcontrol";
  std::ofstream output(filename.c_str());
  output << "Flow control" << std::endl;
  output << "Previous_beta " << beta_n << std::endl;
  output << "Previous_average_velocity " << average_velocity_1n << std::endl;
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
  input >> buffer >> average_velocity_1n;
  input >> buffer >> no_force;
  input >> buffer >> threshold_factor;
}

template class FlowControl<2>;
template class FlowControl<3>;
