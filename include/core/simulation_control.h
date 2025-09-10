// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_simulation_control_h
#define lethe_simulation_control_h

#include <core/bdf.h>
#include <core/parameters.h>
#include <core/sdirk_stage_data.h>

/**
 * @brief The SimulationControl class is responsible for the control of steady-state and transient
 * simulations carried out with Lethe. This base class is a pure virtual class
 * and cannot be instantiated. However, it stores the core variables which are
 * necessary for its serialization (write and read)
 **/

class SimulationControl
{
protected:
  // Time stepping method used for the simulation
  Parameters::SimulationControl::TimeSteppingMethod method;

  // Current method used for the present assembly
  // This is used to differentiate the substeps of SDIRK methods or
  // to start-up BDF simulations
  Parameters::SimulationControl::TimeSteppingMethod assembly_method;

  // Time of the current iteration being solved for
  double current_time;

  // Time of the previous iteration
  double previous_time;

  // Time step linking the previous iteration and the current time
  double time_step;

  // Initial time step given in the parameter
  double initial_time_step;

  // Simulation end time
  double end_time;

  // Boolean to keep the time step for the last iteration regardless of the end
  // time specify. Both for fixed time step and adaptive time step.
  bool time_step_independent_of_end_time;

  // Time step vector. This vector accumulates the time steps of the previous
  // iterations. This is required for multiple steps methods such as the bdfs.
  std::vector<double> time_step_vector;

  // Iteration. Iterations start at 0, but the first actual iteration is
  // iteration 1.
  unsigned int iteration_number;

  // Number of mesh adaptation iteration
  unsigned int number_mesh_adapt;

  // Courant-Friedrich-Levy (CFL) condition
  // Since the simulation control is unaware of the information propagation
  // mechanism (for instance the velocity), the current CFL must be set by the
  // solver itself
  double CFL;

  // Maximal CFL condition. This is used to control adaptative time stepping. In
  // the case of constant time stepping or steady-state simulations, this
  // parameter remains unused.
  double max_CFL;

  // Current value of the norm of the rhs
  double residual;

  // Residual at which the simulation stops when adjoint time-stepping is used
  double stop_tolerance;

  // Number of time steps stored
  // BDF methods require a number of previous time steps. This number is known a
  // priori and depends on the method used. We do not keep all the time steps to
  // prevent the accumulation within a large vector.
  static const unsigned int n_previous_time_steps = 4;

  // Output iteration frequency
  // Controls the output of the simulation results when the output is controlled
  // by the iteration number.
  unsigned int output_iteration_frequency;

  // Log iteration frequency
  // Controls the frequency at which status of the simulation is written to
  // the terminal
  unsigned int log_frequency;

  // Log precision
  // Controls the number of significant digits displayed on the standard outputs
  unsigned int log_precision;

  // Number of mesh subdivision to be used when outputting the results
  unsigned int subdivision;

  // Number of parallel file to generate
  unsigned int group_files;

  // Output name
  std::string output_name;

  // Output path
  std::string output_path;

  // Output boundaries
  // Control if the boundaries of the domain are outputted when writing results
  //
  bool output_boundaries;

  // Indicator to tell if this is the first assembly of a step
  bool first_assembly;

  // The method use to start high order bdf scheme
  Parameters::SimulationControl::BDFStartupMethods bdf_start_method;

  // The time scaling used to do small time-steps at the startup of the
  // simulation
  double startup_timestep_scaling;

  /// BDF coefficients used for time-stepping methods
  Vector<double> bdf_coefs;

  /**
   * @brief Update the BDF coefficients. It is necessary to update the
   * coefficients when there is a change in the time step values or the
   * time stepping scheme.
   **/
  void
  update_bdf_coefficients()
  {
    bdf_coefs = calculate_bdf_coefficients(assembly_method, time_step_vector);
  }

public:
  /**
   * @brief The simulation control class is constructed by a simple parameter structure
   * from which it draws it's arguments. This structure is not kept internally.
   * This means that require information is copied from the struct to the class.
   *
   * @param[in] param Structure of the parameters for the simulation control
   *
   **/

  SimulationControl(const Parameters::SimulationControl &param);

  /**
 * @brief Default destructor.
 **/

  virtual ~SimulationControl()=default;


  /**
   * @brief Return the number of stages associated with the time-stepping method.
   * For SDIRK methods, the number of stages is given by the last digit of the
   * method name. For other methods, the number of stages is 1.
   *
   * @param[in] method The time-stepping method for which the number of stages
   * is requested.
   *
   * @return The number of stages associated with the time-stepping method.
   */
  static unsigned int
  get_number_of_stages(
    const Parameters::SimulationControl::TimeSteppingMethod &method);


  /**
   * @brief Pure virtual function to control the progression of the simulation.
   * As long as integrate returns true, a simulation should proceed. The
   * criteria used to stop/continue the simulation in iterate is a property of
   * each individual time stepping control. For example, steady-state simulation
   * will proceed until the number of required mesh adaptation has been
   *performed.
   **/
  virtual bool
  integrate() = 0;

  /**
   * @brief Establishes if a simulation has reached it's end or not. The concrete
   * implementation of the class decides what is the stopping criteria
   * (iteration number, time_end reached, etc.)
   **/
  virtual bool
  is_at_end() = 0;

  /**
   * @brief Add a time step and stores the previous one in a list.
   *
   *  @param[in] p_timestep the new value of the time step for the present
   *iteration.
   **/
  void
  add_time_step(double p_timestep);


  /**
   * @brief Establish if the iteration is the first iteration or if the simulation
   * has not begun
   *
   */
  bool
  is_at_start()
  {
    return iteration_number <= 1;
  }

  /**
   * @brief Establish if the simulation is a steady-state simulation or no
   *
   */
  bool
  is_steady()
  {
    return method ==
             Parameters::SimulationControl::TimeSteppingMethod::steady ||
           method ==
             Parameters::SimulationControl::TimeSteppingMethod::steady_bdf;
  }

  /**
   * @brief Establish if the method is a bdf method
   *
   * @return true if the method is BDF1, BDF2, BDF3, or steady BDF, false otherwise.
   */
  bool
  is_bdf()
  {
    return method == Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
           method ==
             Parameters::SimulationControl::TimeSteppingMethod::steady_bdf ||
           method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 ||
           method == Parameters::SimulationControl::TimeSteppingMethod::bdf3;
  }

  /**
   * @brief Establish if the method is a sdirk method
   *
   * @return true if the method is sdirk, false otherwise.
   */
  bool
  is_sdirk()
  {
    return method ==
             Parameters::SimulationControl::TimeSteppingMethod::sdirk22 ||
           method ==
             Parameters::SimulationControl::TimeSteppingMethod::sdirk33 ||
           method == Parameters::SimulationControl::TimeSteppingMethod::sdirk43;
  }


  /**
   * @brief Calculate the next value of the time step. The base function returns
   * the value of the time step, but derived class may implement adaptative time
   * stepping
   */
  virtual double
  calculate_time_step()
  {
    return time_step;
  }


  /**
   * @brief print_progress Function that prints the current progress status of the simulation
   *
   * @param[in] pcout the ConditionalOSStream that is use to write
   */
  virtual void
  print_progression(const ConditionalOStream &pcout) = 0;

  /**
   * @brief Check if VTU/PVTU/PVD outputs are disabled
   */
  bool
  output_enabled() const
  {
    return output_iteration_frequency != 0;
  }

  /**
   * @brief Check if the present iteration is an output iteration depending on the
   * output control chosen.
   */
  virtual bool
  is_output_iteration();

  /**
   * @brief Check if the boundaries of the domain should be outputted
   * when writing results
   */
  bool
  get_output_boundaries()
  {
    return output_boundaries;
  }

  /**
   * @brief Check if the present iteration is a verbose iteration where
   * output should be written to the terminal
   */
  virtual bool
  is_verbose_iteration();

  /**
   * @brief Check if this is the first assembly of the present iteration.
   * If it indeed is the first assembly, then the first_assembly is set to
   * false.
   */
  virtual bool
  is_first_assembly()
  {
    bool return_value = first_assembly;
    first_assembly    = false;
    return return_value;
  }

  /**
   * @brief Define the assembly method to use while integrating. Use to start bdf 2 and bdf 3 scheme.
   */
  void
  update_assembly_method();

  /**
   * @brief Set the value of the CFL condition
   *
   * @param[in] p_CFL Value of the CFL condition calculated by the solver.
   */

  void
  set_CFL(const double p_CFL)
  {
    CFL = p_CFL;
  }


  /**
   * @brief Manually force the value of the time step for the present iteration
   *
   * @param[in] new_time_step The new value of the time step.
   * This time step is appended to the time step history
   */
  void
  set_current_time_step(const double new_time_step)
  {
    time_step = new_time_step;
    add_time_step(new_time_step);
  }

  /**
   * @brief Suggest the value of the time step for the next iteration. Note that
   * for adaptative simulations this time step may be altered
   *
   * @param[in] new_time_step The new value of the time step.
   * This time step is not added to the time step vector
   */
  void
  set_suggested_time_step(const double new_time_step)
  {
    time_step = new_time_step;
  }


  /**
   * @brief Provide the value of the residual at the beginning
   * of the iteration to the simulation controller
   *
   * @param[in] new_residual Value of the residual at the beginning
   * of an adjoint time-stepping time step
   */
  void
  provide_residual(const double new_residual)
  {
    residual = new_residual;
  }

  // Relatively trivial getters.

  double
  get_time_step() const
  {
    return time_step;
  }

  /**
   * @brief Get current iteration number.
   *
   * @return Current iteration number.
   */
  unsigned int
  get_iteration_number() const
  {
    return iteration_number;
  }

  std::string
  get_output_name() const
  {
    return output_name;
  }

  std::string
  get_output_path() const
  {
    return output_path;
  }

  unsigned int
  get_group_files() const
  {
    return group_files;
  }

  unsigned int
  get_log_precision() const
  {
    return log_precision;
  }

  double
  get_current_time() const
  {
    return current_time;
  }

  double
  get_previous_time() const
  {
    return previous_time;
  }


  std::vector<double>
  get_time_steps_vector() const
  {
    return time_step_vector;
  }

  double
  get_CFL() const
  {
    return CFL;
  }

  unsigned int
  get_step_number() const
  {
    return iteration_number;
  }

  unsigned int
  get_number_subdivision() const
  {
    return subdivision;
  }

  unsigned int
  get_log_frequency() const
  {
    return log_frequency;
  }

  Parameters::SimulationControl::TimeSteppingMethod
  get_assembly_method() const
  {
    return assembly_method;
  }

  void
  set_assembly_method(Parameters::SimulationControl::TimeSteppingMethod method)
  {
    assembly_method = method;
  }

  unsigned int
  get_number_of_previous_solution_in_assembly() const
  {
    return number_of_previous_solutions(method);
  }

  std::vector<double>
  get_simulation_times() const
  {
    // Create a vector of the previous times
    std::vector<double> times(n_previous_time_steps + 1);
    for (unsigned int i = 0; i < n_previous_time_steps + 1; ++i)
      {
        times[i] = current_time;
        for (unsigned int j = 0; j < i; ++j)
          times[i] -= time_step_vector[j];
      }
    return times;
  }

  Vector<double> const &
  get_bdf_coefficients()
  {
    return bdf_coefs;
  }

  /**
   * @brief Save the simulation control information from the checkpoint file and updates the time step vector, the CFL value, the time and the iteration number.
   *
   * @param[in] prefix The prefix of the checkpoint of the simulation
   */
  virtual void
  save(const std::string &prefix);

  /**
   * @brief Reads the simulation control information from the checkpoint file and updates the time step vector, the CFL value, the time and the iteration number.
   *
   * @param[in] prefix The prefix of the checkpoint of the simulation
   */
  virtual void
  read(const std::string &prefix);

  /**
   * @brief Reads and returns the simulation control information from the checkpoint file filename without updating the simulation control information.
   *
   * @param[in] prefix The prefix of the checkpoint of the simulation
   *
   * @return A vector containing the last checkpointed file and time step.
   */
  std::vector<double>
  get_checkpointed_simulation_control_info(const std::string &prefix);
};


class SimulationControlTransient : public SimulationControl
{
protected:
  // Enable adaptative time stepping
  bool adapt;

  // Time step scaling for adaptative time stepping
  double adaptative_time_step_scaling;

  // Max time step
  double max_dt;

  // Time last output
  double time_last_output;

  // Specific output frequency for time output control
  double output_time_frequency;

  // Specific output times for time output control
  std::vector<double> output_times_vector;

  // Counter to move between output times given in previous vector
  unsigned int output_times_counter;

  // Variable to check whether we still have specific times to check in output
  // times vector
  bool no_more_output_times;

  // Time interval for output of transient iterations either with time output
  // control or iterations control
  std::vector<double> output_time_interval;

  // Output control type: iteration or type
  Parameters::SimulationControl::OutputControl output_control;

  /**
   * @brief Calculates the next value of the time step. If adaptation
   * is enabled, the time step is calculated in order to ensure
   * that the CFL condition is bound by the maximal CFL value.
   * The new time step is equal to adaptative_time_step_scaling * the previous
   * time step. If the time_step_independent_of_end_time is set to false, the
   * time step is asjusted to meet exactly the end time; the default is to not
   * modify the time step.
   */
  virtual double
  calculate_time_step() override;

public:
  SimulationControlTransient(const Parameters::SimulationControl &param);

  virtual void
  print_progression(const ConditionalOStream &pcout) override;

  /**
   * @brief Proceeds with the simulation until the end condition is reached
   */
  virtual bool
  integrate() override;

  /**
   * @brief Ends the simulation when the end time is reached
   */
  virtual bool
  is_at_end() override;

  /**
   * @brief Output iterations are calculated based on the value of the output time
   * or the output interval and the output frequency within the interval.
   */
  virtual bool
  is_output_iteration() override;

  /**
   * @brief Save the simulation control information from the checkpoint file and
   * updates the time step vector, the CFL value, the time, the iteration number
   * and the index of the output times vector if time control is used.
   *
   * @param prefix The prefix of the checkpoint of the simulation
   */
  void
  save(const std::string &prefix) override;

  /**
   * @brief Reads the simulation control information from the checkpoint file and
   * updates the time step vector, the CFL value, the time, the iteration number
   * and the index of the output times vector if time control is used. Allows to
   * change time step if adaptive time stepping is disabled.
   *
   * @param prefix The prefix of the checkpoint of the simulation
   */
  void
  read(const std::string &prefix) override;
};

/**
 * @brief Transient simulation control tailored around the discrete element method
 */
class SimulationControlTransientDEM : public SimulationControlTransient
{
public:
  SimulationControlTransientDEM(const Parameters::SimulationControl &param);

  virtual void
  print_progression(const ConditionalOStream &pcout) override;
};

class SimulationControlSteady : public SimulationControl
{
public:
  SimulationControlSteady(const Parameters::SimulationControl &param);

  virtual void
  print_progression(const ConditionalOStream &pcout) override;

  /**
   * @brief Proceeds with the simulation until the end condition is reached
   */
  virtual bool
  integrate() override;

  /**
   * @brief Ends the simulation when the number of mesh adaptation is reached
   */
  virtual bool
  is_at_end() override;
};

class SimulationControlAdjointSteady : public SimulationControlTransient
{
public:
  SimulationControlAdjointSteady(const Parameters::SimulationControl &param);

  virtual void
  print_progression(const ConditionalOStream &pcout) override;

  /**
   * @brief Ends the simulation when the desired residual is reached
   */
  virtual bool
  is_at_end() override;
};

#endif
