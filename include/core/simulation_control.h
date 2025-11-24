// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file simulation_control.h
 * @brief Simulation control classes for managing time-stepping and iteration control in Lethe
 *
 * This file defines the simulation control hierarchy used in Lethe to manage
 * steady-state and transient simulations. The SimulationControl base class and
 * its derived classes handle time-stepping methods (including BDF and SDIRK),
 * adaptive time stepping, output control, and simulation progression logic.
 */

#ifndef lethe_simulation_control_h
#define lethe_simulation_control_h

#include <core/bdf.h>
#include <core/parameters.h>

#include <deal.II/particles/particle_handler.h>

/**
 * @brief Base class for controlling steady-state and transient simulations in Lethe
 *
 * The SimulationControl class is responsible for managing the progression of
 * simulations carried out with Lethe. This is a pure virtual base class that
 * cannot be instantiated directly. It stores core variables necessary for
 * time-stepping control, output management, and serialization.
 *
 * The class supports multiple time-stepping methods including:
 * - Backward Differentiation Formula (BDF) methods: BDF1, BDF2, BDF3
 * - Singly Diagonally Implicit Runge-Kutta (SDIRK) methods
 * - Steady-state methods
 *
 * Derived classes implement specific simulation control strategies for
 * transient, steady-state, adjoint, and specialized simulations (e.g., DEM, ray
 * tracing).
 *
 * @note This class handles serialization through save() and read() methods
 * to enable checkpointing and restart capabilities.
 */
class SimulationControl
{
protected:
  /// Time stepping method used for the simulation (BDF, SDIRK, steady, etc.)
  Parameters::SimulationControl::TimeSteppingMethod method;

  /**
   * @brief Current method used for the present assembly
   *
   * This is used to differentiate the substeps of SDIRK methods or
   * to start-up BDF simulations. For example, BDF2 and BDF3 schemes
   * may start with lower-order methods before transitioning to the
   * full-order scheme.
   */
  Parameters::SimulationControl::TimeSteppingMethod assembly_method;

  /// Time of the current iteration being solved for
  double current_time;

  /// Time of the previous iteration
  double previous_time;

  /// Time step linking the previous iteration and the current time
  double time_step;

  /// Initial time step given in the parameter file
  double initial_time_step;

  /// Simulation end time
  double end_time;

  /**
   * @brief Boolean to keep the time step for the last iteration regardless of the end time
   *
   * When true, the last time step is not adjusted to match the end time
   * exactly. Applies to both fixed time step and adaptive time step
   * simulations.
   */
  bool time_step_independent_of_end_time;

  /**
   * @brief Time step vector storing previous time steps
   *
   * This vector accumulates the time steps of the previous iterations.
   * This is required for multiple-step methods such as BDF schemes.
   * The vector size is limited by n_previous_time_steps to prevent
   * accumulation within a large vector.
   */
  std::vector<double> time_step_vector;

  /**
   * @brief Iteration number
   *
   * Iterations start at 0, but the first actual iteration is iteration 1.
   */
  unsigned int iteration_number;

  /// Number of mesh adaptation iterations performed
  unsigned int number_mesh_adapt;

  /**
   * @brief Courant-Friedrich-Levy (CFL) condition value
   *
   * Since the simulation control is unaware of the information propagation
   * mechanism (for instance the velocity), the current CFL must be set by the
   * solver itself using set_CFL().
   */
  double CFL;

  /**
   * @brief Maximal CFL condition for adaptive time stepping
   *
   * This is used to control adaptive time stepping. In the case of constant
   * time stepping or steady-state simulations, this parameter remains unused.
   */
  double max_CFL;

  /// Current value of the norm of the right-hand side residual
  double residual;

  /**
   * @brief Residual tolerance at which the simulation stops
   *
   * This is used when adjoint time-stepping is employed.
   */
  double stop_tolerance;

  /**
   * @brief Number of time steps stored
   *
   * BDF methods require a number of previous time steps. This number is known a
   * priori and depends on the method used. We do not keep all the time steps to
   * prevent the accumulation within a large vector.
   *
   * @note Currently set to 4 to support up to BDF3 methods.
   */
  static constexpr unsigned int n_previous_time_steps = 4;

  /**
   * @brief Output iteration frequency
   *
   * Controls the output of the simulation results when the output is controlled
   * by the iteration number. A value of 0 disables output.
   */
  unsigned int output_iteration_frequency;

  /**
   * @brief Log iteration frequency
   *
   * Controls the frequency at which the status of the simulation is written to
   * the terminal.
   */
  unsigned int log_frequency;

  /**
   * @brief Log precision
   *
   * Controls the number of significant digits displayed on the standard
   * outputs.
   */
  unsigned int log_precision;

  /// Number of mesh subdivisions to be used when outputting the results
  unsigned int subdivision;

  /// Number of parallel files to generate for output
  unsigned int group_files;

  /// Output file name prefix
  std::string output_name;

  /// Output directory path
  std::string output_path;

  /**
   * @brief Output boundaries flag
   *
   * Controls if the boundaries of the domain are outputted when writing
   * results.
   */
  bool output_boundaries;

  /// Indicator to tell if this is the first assembly of a time step
  bool first_assembly;

  /// The method used to start high-order BDF schemes
  Parameters::SimulationControl::BDFStartupMethods bdf_start_method;

  /**
   * @brief The time scaling used for small time-steps at the startup of the simulation
   *
   * This scaling factor is applied during the startup phase when using
   * high-order BDF methods to gradually transition to the full-order scheme.
   */
  double startup_timestep_scaling;

  /// BDF coefficients used for time-stepping methods
  Vector<double> bdf_coefs;

  /**
   * @brief Update the BDF coefficients
   *
   * It is necessary to update the coefficients when there is a change in the
   * time step values or the time stepping scheme. This method recalculates
   * the BDF coefficients based on the current assembly method and time step
   * history.
   *
   * @see calculate_bdf_coefficients
   */
  void
  update_bdf_coefficients()
  {
    bdf_coefs = calculate_bdf_coefficients(assembly_method, time_step_vector);
  }

public:
  /**
   * @brief Construct a SimulationControl object from parameter structure
   *
   * The simulation control class is constructed from a parameter structure
   * from which it draws its arguments. This structure is not kept internally.
   * This means that required information is copied from the struct to the
   * class.
   *
   * @param[in] param Structure of the parameters for the simulation control
   */
  SimulationControl(const Parameters::SimulationControl &param);

  /**
   * @brief Default destructor.
   */
  virtual ~SimulationControl() = default;


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
   * @brief Pure virtual function to control the progression of the simulation
   *
   * As long as integrate returns true, a simulation should proceed. The
   * criteria used to stop/continue the simulation in integrate is a property of
   * each individual time stepping control. For example, steady-state simulation
   * will proceed until the number of required mesh adaptations has been
   * performed.
   *
   * @return true if the simulation should continue, false otherwise
   */
  virtual bool
  integrate() = 0;

  /**
   * @brief Establishes if a simulation has reached its end
   *
   * The concrete implementation of the class decides what is the stopping
   * criteria (iteration number, time_end reached, residual tolerance, etc.)
   *
   * @return true if the simulation has reached its end condition, false otherwise
   */
  virtual bool
  is_at_end() = 0;

  /**
   * @brief Add a time step and store the previous one in the history
   *
   * This method updates the time step history by adding the new time step
   * to the time_step_vector. The history is maintained for BDF methods that
   * require knowledge of previous time steps.
   *
   * @param[in] p_timestep The new value of the time step for the present
   * iteration
   */
  void
  add_time_step(double p_timestep);


  /**
   * @brief Establish if the iteration is the first iteration or if the simulation has not begun
   *
   * @return true if iteration_number <= 1, false otherwise
   */
  bool
  is_at_start() const
  {
    return iteration_number <= 1;
  }

  /**
   * @brief Establish if the simulation is a steady-state simulation
   *
   * @return true if the method is steady or steady_bdf, false otherwise
   */
  bool
  is_steady() const
  {
    return method ==
             Parameters::SimulationControl::TimeSteppingMethod::steady ||
           method ==
             Parameters::SimulationControl::TimeSteppingMethod::steady_bdf;
  }

  /**
   * @brief Establish if the method is a BDF method
   *
   * @return true if the method is BDF1, BDF2, BDF3, or steady BDF, false otherwise
   */
  bool
  is_bdf() const
  {
    return method == Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
           method ==
             Parameters::SimulationControl::TimeSteppingMethod::steady_bdf ||
           method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 ||
           method == Parameters::SimulationControl::TimeSteppingMethod::bdf3;
  }

  /**
   * @brief Establish if the method is an SDIRK method
   *
   * @return true if the method is sdirk22, sdirk33, or sdirk43, false otherwise
   */
  bool
  is_sdirk() const
  {
    return method ==
             Parameters::SimulationControl::TimeSteppingMethod::sdirk22 ||
           method ==
             Parameters::SimulationControl::TimeSteppingMethod::sdirk33 ||
           method == Parameters::SimulationControl::TimeSteppingMethod::sdirk43;
  }


  /**
   * @brief Calculate the next value of the time step
   *
   * The base function returns the value of the current time step, but derived
   * classes may implement adaptive time stepping based on CFL conditions or
   * other criteria.
   *
   * @return The calculated time step value
   */
  virtual double
  calculate_time_step()
  {
    return time_step;
  }


  /**
   * @brief Print the current progress status of the simulation
   *
   * Pure virtual function that prints simulation progress information to the
   * terminal. The specific format and content depend on the derived class
   * implementation.
   *
   * @param[in] pcout Parallel conditional output stream used to print the
   * information
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
   * @brief Check if the present iteration is an output iteration
   *
   * Determines whether output should be written at the current iteration
   * based on the output control chosen (iteration-based or time-based).
   *
   * @return true if output should be written at this iteration, false otherwise
   */
  virtual bool
  is_output_iteration();

  /**
   * @brief Check if the boundaries of the domain should be outputted when writing results
   *
   * @return true if boundaries should be included in output, false otherwise
   */
  bool
  get_output_boundaries() const
  {
    return output_boundaries;
  }

  /**
   * @brief Check if the present iteration is a verbose iteration
   *
   * Determines if output should be written to the terminal at the current
   * iteration based on the log_frequency setting.
   *
   * @return true if verbose output should be written at this iteration, false otherwise
   */
  virtual bool
  is_verbose_iteration();

  /**
   * @brief Check if this is the first assembly of the present iteration
   *
   * If it indeed is the first assembly, then the first_assembly flag is set to
   * false. This function is used when providing the residual to the simulation
   * control object when using steady-bdf methods. Since the residual must be
   * provided once per time-step (and at the beginning), this function is used
   * to identify that.
   *
   * @return true if this is the first assembly, false otherwise
   *
   * @note This method has a side effect: it sets first_assembly to false after
   * being called when it returns true.
   */
  virtual bool
  is_first_assembly()
  {
    bool return_value = first_assembly;
    first_assembly    = false;
    return return_value;
  }

  /**
   * @brief Define the assembly method to use while integrating
   *
   * Used to start BDF2 and BDF3 schemes with lower-order methods before
   * transitioning to the full-order scheme. This ensures stability during
   * the initial time steps.
   */
  void
  update_assembly_method();

  /**
   * @brief Set the value of the CFL condition
   *
   * The solver calculates and provides the CFL number to the simulation control
   * for use in adaptive time stepping.
   *
   * @param[in] p_CFL Value of the CFL condition calculated by the solver
   */
  void
  set_CFL(const double p_CFL)
  {
    CFL = p_CFL;
  }


  /**
   * @brief Manually force the value of the time step for the present iteration
   *
   * Sets the time step to a specific value. This time step is appended to the
   * time step history.
   *
   * @param[in] new_time_step The new value of the time step
   *
   * @pre new_time_step must be positive (cannot go backward in time)
   */
  void
  set_current_time_step(const double new_time_step)
  {
    Assert(
      time_step > 0,
      ExcMessage(
        "You are trying to set a null or negative time-step in a SimulationControl. This is now allowed, we cannot go backward in time."));
    time_step = new_time_step;
  }

  /**
   * @brief Provide the value of the residual at the beginning of the iteration to the simulation controller
   *
   * Used primarily in adjoint time-stepping to track convergence.
   *
   * @param[in] new_residual Value of the residual at the beginning of an
   * adjoint time-stepping time step
   */
  void
  provide_residual(const double new_residual)
  {
    residual = new_residual;
  }

  /**
   * @brief Get current time step
   *
   * @return Current time step value
   */
  double
  get_time_step() const
  {
    return time_step;
  }

  /**
   * @brief Get current iteration number
   *
   * @return Current iteration number
   */
  unsigned int
  get_iteration_number() const
  {
    return iteration_number;
  }

  /**
   * @brief Get output file name prefix
   *
   * @return Output file name prefix
   */
  std::string
  get_output_name() const
  {
    return output_name;
  }

  /**
   * @brief Get output directory path
   *
   * @return Output directory path
   */
  std::string
  get_output_path() const
  {
    return output_path;
  }

  /**
   * @brief Get number of parallel files to generate for output
   *
   * @return Number of parallel output files
   */
  unsigned int
  get_group_files() const
  {
    return group_files;
  }

  /**
   * @brief Get log precision
   *
   * @return Number of significant digits displayed on the standard outputs
   */
  unsigned int
  get_log_precision() const
  {
    return log_precision;
  }

  /**
   * @brief Get current simulation time
   *
   * @return Current simulation time
   */
  double
  get_current_time() const
  {
    return current_time;
  }

  /**
   * @brief Get previous simulation time
   *
   * @return Previous simulation time
   */
  double
  get_previous_time() const
  {
    return previous_time;
  }

  /**
   * @brief Get time step history vector
   *
   * @return Vector containing previous time steps
   */
  std::vector<double>
  get_time_steps_vector() const
  {
    return time_step_vector;
  }

  /**
   * @brief Get current CFL number
   *
   * @return Current CFL condition value
   */
  double
  get_CFL() const
  {
    return CFL;
  }

  /**
   * @brief Get step number (alias for get_iteration_number)
   *
   * @return Current iteration/step number
   */
  unsigned int
  get_step_number() const
  {
    return iteration_number;
  }

  /**
   * @brief Get number of mesh subdivisions for output
   *
   * @return Number of mesh subdivisions to be used when outputting the results
   */
  unsigned int
  get_number_subdivision() const
  {
    return subdivision;
  }

  /**
   * @brief Get log frequency
   *
   * @return Frequency at which the status of the simulation is written to the terminal
   */
  unsigned int
  get_log_frequency() const
  {
    return log_frequency;
  }

  /**
   * @brief Get current assembly method
   *
   * @return Current time stepping method used for assembly
   */
  Parameters::SimulationControl::TimeSteppingMethod
  get_assembly_method() const
  {
    return assembly_method;
  }

  /**
   * @brief Set the assembly method
   *
   * @param[in] a_method The time stepping method to use for assembly
   */
  void
  set_assembly_method(
    const Parameters::SimulationControl::TimeSteppingMethod a_method)
  {
    assembly_method = a_method;
  }

  /**
   * @brief Get number of previous solutions required in assembly
   *
   * Returns the number of previous solution vectors required for the current
   * time stepping method (e.g., 1 for BDF1, 2 for BDF2, etc.).
   *
   * @return Number of previous solutions needed for the current method
   */
  unsigned int
  get_number_of_previous_solution_in_assembly() const
  {
    return number_of_previous_solutions(method);
  }

  /**
   * @brief Get simulation times for current and previous time steps
   *
   * Creates a vector containing the current time and previous times based on
   * the time step history. The times are calculated by subtracting the
   * accumulated time steps from the current time.
   *
   * @return Vector of simulation times (current and previous)
   */
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

  /**
   * @brief Get BDF coefficients
   *
   * Returns the BDF coefficients used for time integration. These coefficients
   * are updated by update_bdf_coefficients() when the time step or assembly
   * method changes.
   *
   * @return Const reference to the BDF coefficients vector
   */
  Vector<double> const &
  get_bdf_coefficients()
  {
    return bdf_coefs;
  }

  /**
   * @brief Indicate if the simulation uses adaptive time stepping or not
   *
   * Pure virtual function to be implemented by derived classes based on their
   * specific time stepping strategy.
   *
   * @return true if the simulation has adaptive time stepping, false otherwise
   */
  virtual bool
  is_adaptive_time_stepping() const = 0;

  /**
   * @brief Save the simulation control information to a checkpoint file
   *
   * Saves the time step vector, the CFL value, the time and the iteration
   * number to enable restart capabilities.
   *
   * @param[in] prefix The prefix of the checkpoint file for the simulation
   */
  virtual void
  save(const std::string &prefix);

  /**
   * @brief Read the simulation control information from a checkpoint file
   *
   * Reads and updates the time step vector, the CFL value, the time and the
   * iteration number from a checkpoint file to restart the simulation.
   *
   * @param[in] prefix The prefix of the checkpoint file for the simulation
   */
  virtual void
  read(const std::string &prefix);

  /**
   * @brief Read and return simulation control information from checkpoint without updating
   *
   * Reads and returns the simulation control information from the checkpoint
   * file without updating the internal simulation control state. Useful for
   * querying checkpoint information without modifying the current simulation
   * state.
   *
   * @param[in] prefix The prefix of the checkpoint file for the simulation
   *
   * @return A vector containing the last checkpointed iteration number and time step
   */
  std::vector<double>
  get_checkpointed_simulation_control_info(const std::string &prefix);
};


/**
 * @brief Simulation control for transient simulations
 *
 * This class extends SimulationControl to handle transient (time-dependent)
 * simulations. It implements adaptive time stepping based on CFL conditions,
 * time-based output control, and ensures the simulation progresses until the
 * end time is reached.
 *
 * Key features:
 * - Adaptive time stepping with CFL-based control
 * - Time-based or iteration-based output control
 * - Specific output times capability
 * - Maximum time step limiting
 */
class SimulationControlTransient : public SimulationControl
{
protected:
  /// Enable adaptive time stepping
  bool adapt;

  /**
   * @brief Time step scaling factor for adaptive time stepping
   *
   * Used to scale the time step based on the ratio of max_CFL to current CFL
   * when adaptive time stepping is enabled.
   */
  double adaptative_time_step_scaling;

  /// Maximum allowed time step
  double max_dt;

  /// Time of the last output
  double time_last_output;

  /// Output frequency for time-based output control
  double output_time_frequency;

  /**
   * @brief Specific output times for time output control
   *
   * Vector containing specific times at which output should be generated.
   */
  std::vector<double> output_times_vector;

  /// Counter to move between output times given in output_times_vector
  unsigned int output_times_counter;

  /**
   * @brief Flag indicating whether there are more output times to check
   *
   * Set to true when all specific output times in output_times_vector have been
   * reached.
   */
  bool no_more_output_times;

  /**
   * @brief Flag to override the time-step with the set value upon restart
   */
  bool override_time_step_on_restart;

  /**
   * @brief Time interval for output of transient iterations
   *
   * Used with either time output control or iteration control to define
   * intervals for output generation.
   */
  std::vector<double> output_time_interval;

  /// Output control type: iteration-based or time-based
  Parameters::SimulationControl::OutputControl output_control;

  /**
   * @brief Calculate the next value of the time step
   *
   * If adaptation is enabled, the time step is calculated to ensure that the
   * CFL condition is bound by the maximal CFL value. The new time step is
   * scaled by adaptative_time_step_scaling. If the
   * time_step_independent_of_end_time is set to false, the time step is
   * adjusted to meet exactly the end time; the default is to not modify the
   * time step.
   *
   * @return The calculated time step value
   */
  virtual double
  calculate_time_step() override;

public:
  /**
   * @brief Construct a transient simulation control object
   *
   * @param[in] param Structure of the parameters for the simulation control
   */
  SimulationControlTransient(const Parameters::SimulationControl &param);

  /**
   * @brief Print the current progress status of the simulation
   *
   * @param[in] pcout Parallel conditional output stream used to print the
   * information
   */
  virtual void
  print_progression(const ConditionalOStream &pcout) override;

  /**
   * @brief Indicate if the simulation uses adaptive time stepping or not
   *
   * @return true if the simulation has adaptive time stepping, false otherwise
   */
  bool
  is_adaptive_time_stepping() const override
  {
    return adapt;
  }

  /**
   * @brief Proceed with the simulation until the end condition is reached
   *
   * Advances the simulation by incrementing time and iteration number.
   *
   * @return true if the simulation should continue, false otherwise
   */
  virtual bool
  integrate() override;

  /**
   * @brief End the simulation when the end time is reached
   *
   * @return true if the current time has reached or exceeded the end time, false otherwise
   */
  virtual bool
  is_at_end() override;

  /**
   * @brief Check if the current iteration is an output iteration
   *
   * Output iterations are determined based on the value of the output time,
   * the output interval, and the output frequency within the interval.
   *
   * @return true if output should be written at this iteration, false otherwise
   */
  virtual bool
  is_output_iteration() override;

  /**
   * @brief Save the simulation control information to a checkpoint file
   *
   * Saves the time step vector, the CFL value, the time, the iteration number
   * and the index of the output times vector if time control is used.
   *
   * @param[in] prefix The prefix of the checkpoint file for the simulation
   */
  void
  save(const std::string &prefix) override;

  /**
   * @brief Read the simulation control information from a checkpoint file
   *
   * Reads and updates the time step vector, the CFL value, the time, the
   * iteration number and the index of the output times vector if time control
   * is used. Allows changing the time step if adaptive time stepping is
   * disabled.
   *
   * @param[in] prefix The prefix of the checkpoint file for the simulation
   */
  void
  read(const std::string &prefix) override;
};

/**
 * @brief Transient simulation control tailored for the Discrete Element Method (DEM)
 *
 * This class specializes SimulationControlTransient for DEM simulations.
 * It provides DEM-specific progress reporting and output formatting.
 */
class SimulationControlTransientDEM : public SimulationControlTransient
{
public:
  /**
   * @brief Construct a DEM transient simulation control object
   *
   * @param[in] param Structure of the parameters for the simulation control
   */
  SimulationControlTransientDEM(const Parameters::SimulationControl &param);

  /**
   * @brief Print the current progress status of the DEM simulation
   *
   * Provides DEM-specific formatting and information in the progress output.
   *
   * @param[in] pcout Parallel conditional output stream used to print the
   * information
   */
  virtual void
  print_progression(const ConditionalOStream &pcout) override;
};

/**
 * @brief Simulation control for steady-state simulations
 *
 * This class handles steady-state simulations where time integration is not
 * required. The simulation progresses through iterations until the specified
 * number of mesh adaptations is reached or the maximum number of iterations
 * is completed.
 *
 * @note Steady-state simulations do not use adaptive time stepping.
 */
class SimulationControlSteady : public SimulationControl
{
public:
  /**
   * @brief Construct a steady-state simulation control object
   *
   * @param[in] param Structure of the parameters for the simulation control
   */
  SimulationControlSteady(const Parameters::SimulationControl &param);

  /**
   * @brief Print the current progress status of the steady-state simulation
   *
   * @param[in] pcout Parallel conditional output stream used to print the
   * information
   */
  virtual void
  print_progression(const ConditionalOStream &pcout) override;

  /**
   * @brief Indicate if the simulation uses adaptive time stepping
   *
   * @return false (steady-state simulations do not use adaptive time stepping)
   */
  bool
  is_adaptive_time_stepping() const override
  {
    return false;
  }

  /**
   * @brief Proceed with the simulation until the end condition is reached
   *
   * Advances the simulation by incrementing the iteration number.
   *
   * @return true if the simulation should continue, false otherwise
   */
  virtual bool
  integrate() override;

  /**
   * @brief End the simulation when the number of mesh adaptations is reached
   *
   * @return true if the required number of mesh adaptations has been completed, false otherwise
   */
  virtual bool
  is_at_end() override;
};

/**
 * @brief Simulation control for adjoint steady-state simulations
 *
 * This class manages adjoint simulations that use a pseudo-transient approach
 * to converge to a steady state. The simulation terminates when the residual
 * reaches a specified tolerance (stop_tolerance).
 *
 * The class inherits from SimulationControlTransient to leverage time-stepping
 * mechanisms, but applies them in the context of an adjoint problem where time
 * is a pseudo-parameter for convergence.
 */
class SimulationControlAdjointSteady : public SimulationControlTransient
{
public:
  /**
   * @brief Construct an adjoint steady-state simulation control object
   *
   * @param[in] param Structure of the parameters for the simulation control
   */
  SimulationControlAdjointSteady(const Parameters::SimulationControl &param);

  /**
   * @brief Print the current progress status of the adjoint simulation
   *
   * @param[in] pcout Parallel conditional output stream used to print the
   * information
   */
  virtual void
  print_progression(const ConditionalOStream &pcout) override;

  /**
   * @brief End the simulation when the desired residual tolerance is reached
   *
   * @return true if the residual has reached the stop_tolerance, false otherwise
   */
  virtual bool
  is_at_end() override;
};


/**
 * @brief Simulation control for ray tracing simulations
 *
 * This class manages ray tracing simulations where photons are traced through
 * a computational domain. The simulation terminates when all photons have
 * exited the triangulation or been absorbed.
 *
 * Ray tracing simulations do not use time stepping in the traditional sense;
 * instead, iterations represent photon propagation steps until all photons
 * have completed their trajectories.
 */
class SimulationControlRayTracing : public SimulationControl
{
public:
  /**
   * @brief Construct a ray tracing simulation control object
   *
   * @param[in] param Structure of the parameters for the simulation control
   * @param[in] input_photon_handler The particle handler containing the photons
   * to be traced
   */
  SimulationControlRayTracing(
    const Parameters::SimulationControl &param,
    Particles::ParticleHandler<3>       &input_photon_handler);

  /**
   * @brief Proceed with the simulation until the end condition is reached
   *
   * Advances the ray tracing simulation by incrementing the iteration number.
   *
   * @return true if the simulation should continue, false otherwise
   */
  virtual bool
  integrate() override;

  /**
   * @brief End the simulation when there are no longer any photons in the triangulation
   *
   * @return true if all photons have exited or been absorbed, false otherwise
   */
  bool
  is_at_end() override;

  /**
   * @brief Indicate if the simulation uses adaptive time stepping
   *
   * @return false (ray tracing simulations do not use adaptive time stepping)
   */
  bool
  is_adaptive_time_stepping() const override
  {
    return false;
  }

  /**
   * @brief Print the current progress status of the ray tracing simulation
   *
   * @param[in] pcout Parallel conditional output stream used to print the
   * information
   */
  void
  print_progression(const ConditionalOStream &pcout) override;

  /**
   * @brief Reference to the photon handler used during the ray-tracing simulation
   *
   * This handler manages the photons being traced through the domain.
   */
  Particles::ParticleHandler<3> &photon_handler;
};

#endif
