// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_postprocessing_probes_h
#define lethe_postprocessing_probes_h

#include <core/parameters.h>
#include <core/simulation_control.h>
#include <core/utilities.h>
#include <core/vector.h>

#include <ranges>

template <int dim>
class PostprocessingProbes
{
public:
  /**
   * @brief Object that evaluated variables (e.g., velocity) at specified points
   * in the domain.
   *
   * @param[in] simulation_control SimulationControl object which is used to
   * obtain the current time and current iteration number.
   * @param[in] probing_points_parameters Set of parameters describing the
   * specified probing point and variables of interest.
   * @param[in] output_folder Path to the output folder where the .dat files
   * are saved.
   *
   * @remark At the moment, probing for only the velocity, pressure, phase, and
   * temperature variables are implemented. Implementation for other variables
   * will follow.
   */
  PostprocessingProbes(std::shared_ptr<SimulationControl> simulation_control,
                       const Parameters::PostProcessing<dim>::ProbingPoints
                                                   &probing_points_parameters,
                       const double                 output_frequency,
                       const std::string           &output_folder,
                       const Parameters::Verbosity &verbosity);

  /**
   * @brief Prepares probe tables by setting the entry formats and precisions
   * for all specified probes.
   */
  inline void
  prepare_probe_tables()
  {
    const auto &probing_points_per_variable =
      probing_points_parameters.probing_points_per_variable;

    for (auto &current_table : probe_tables)
      {
        current_table.declare_column("time");
        current_table.set_precision("time", table_precision);
        current_table.set_scientific("time", true);
      }

    auto keys_range = std::views::keys(probing_points_per_variable);

    for (auto variable_it = keys_range.begin(); variable_it != keys_range.end();
         ++variable_it)
      {
        const auto &probing_points_of_variable =
          probing_points_parameters.probing_points_per_variable.at(
            *variable_it);

        // Go trough tables and prepare them
        for (unsigned int i = 0; i < probing_points_of_variable.ids.size(); ++i)
          {
            TableHandler     &current_table = probe_tables[i];
            const std::string column_name =
              get_variable_string_with_underscores(*variable_it);

            if (implemented_scalar_variables.contains(*variable_it))
              {
                current_table.declare_column(column_name);
                current_table.set_precision(column_name, table_precision);
                current_table.set_scientific(column_name, true);
              }
            else if (implemented_vector_variables.contains(*variable_it))
              {
                // Only velocity is implemented at the moment
                current_table.declare_column(column_name + "_x");
                current_table.set_precision(column_name + "_x",
                                            table_precision);
                current_table.set_scientific(column_name + "_x", true);
                current_table.declare_column(column_name + "_y");
                current_table.set_precision(column_name + "_y",
                                            table_precision);
                current_table.set_scientific(column_name + "_y", true);
                if constexpr (dim == 3)
                  {
                    current_table.declare_column(column_name + "_z");
                    current_table.set_precision(column_name + "_z",
                                                table_precision);
                    current_table.set_scientific(column_name + "_z", true);
                  }
              }
          }
      }
  }

  /**
   * @brief Evaluates scalar variables (e.g., pressure) at specified points
   * in the domain using evaluate_values_at_points and fills appropriate probe
   * tables.
   *
   * @tparam VectorType Type of vector of the solution vector.
   *
   * @param[in] triangulation Triangulation object.
   * @param[in] mapping Mapping of the domain.
   * @param[in] dof_handler DoF handler associated to the solution field.
   * @param[in] present_solution Vector containing the scalar solution field.
   * @param[in] variable Variable of interest.
   * @param[in] pcout Parallel console output stream.
   * @param[out] evaluated_scalar_values Vector of evaluated scalar values.
   */
  template <typename VectorType>
  void
  postprocess_probes(const Triangulation<dim> &triangulation,
                     const Mapping<dim>       &mapping,
                     const DoFHandler<dim>    &dof_handler,
                     const VectorType         &present_solution,
                     const Variable            variable,
                     const ConditionalOStream &pcout,
                     std::vector<double>      &evaluated_scalar_values);

  /**
   * @brief Evaluates vector variables (e.g., velocity) at specified points
   * in the domain using evaluate_values_at_points and fills appropriate probe
   * tables.
   *
   * @tparam VectorType Type of vector of the solution vector.
   *
   * @param[in] triangulation Triangulation object.
   * @param[in] mapping Mapping of the domain.
   * @param[in] dof_handler DoF handler associated to the solution field.
   * @param[in] present_solution Vector containing the vector solution field.
   * @param[in] variable Variable of interest.
   * @param[in] pcout Parallel console output stream.
   * @param[out] evaluated_vector_values Vector of evaluated vector values.
   *
   * @remark At the moment, only the velocity variable is implemented.
   */
  template <typename VectorType>
  void
  postprocess_probes(
    const Triangulation<dim>            &triangulation,
    const Mapping<dim>                  &mapping,
    const DoFHandler<dim>               &dof_handler,
    const VectorType                    &present_solution,
    const Variable                       variable,
    const ConditionalOStream            &pcout,
    std::vector<Tensor<1, dim, double>> &evaluated_vector_values);

  /**
   * @brief Writes probe table into output .dat file.
   */
  inline void
  write_probing_points_tables()
  {
    if (simulation_control->get_iteration_number() % output_frequency == 0)
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          for (unsigned int i = 0; i < probe_tables.size(); ++i)
            {
              std::string file_path =
                output_folder +
                probing_points_parameters.probing_points_output_names[i] +
                ".dat";
              std::ofstream output(file_path.c_str());

              probe_tables[i].write_text(output);
            }
        }
  }

  /**
   * @brief Clears the content of the set indicating if the current time has
   * been added to a probe table. This is reset at every postprocessing time
   * step.
   */
  inline void
  reset_current_time_set()
  {
    probe_tables_with_current_time.clear();
  };

  /**
   * @brief Returns a vector of references to TableHandler objects that needs to
   * be serialized/deserialized for the postprocessed probes.
   *
   * @return Structure containing a vector of references to TableHandler objects
   * that needs to be serialized/deserialized for the postprocessed probes, and
   * their corresponding filenames.
   */
  std::vector<OutputStructTableHandler>
  gather_tables();

  /**
   * @brief Serializes all probe tables.
   */
  inline void
  write_checkpoint()
  {
    const std::vector<OutputStructTableHandler> &table_output_structs =
      this->gather_tables();
    serialize_tables_vector(table_output_structs, MPI_COMM_WORLD);
  }

  /**
   * @brief Deserializes all probe tables.
   */
  inline void
  read_checkpoint()
  {
    std::vector<OutputStructTableHandler> table_output_structs =
      this->gather_tables();
    deserialize_tables_vector(table_output_structs, MPI_COMM_WORLD);
  }

  /// Set of valid scalar variables
  const std::set<Variable> implemented_scalar_variables =
    {Variable::pressure, Variable::phase, Variable::temperature};

  /// Set of valid vector variables
  const std::set<Variable> implemented_vector_variables = {Variable::velocity};

private:
  /**
   * @brief Checks if the table of a specified probing point (@p probe_id) has
   * already a row entry for the current time step.
   *
   * @param[in] probe_id Identifier of the probing point investigated.
   *
   * @return
   *  - @p true if there already is an entry for the current time step.
   *  - @p false if there is no row for the current time step.
   */
  inline bool
  check_if_table_contains_current_time(const unsigned int &probe_id)
  {
    if (probe_tables_with_current_time.contains(probe_id))
      return true;
    else // No entry yet for the current time step
      {
        // Add the ID of the probing point to the set
        probe_tables_with_current_time.insert(probe_id);
        return false;
      }
  };

  /// Simulation control for current time and current iteration number.
  std::shared_ptr<SimulationControl> simulation_control;

  /** Set of parameters describing the specified probing point and variables of
   * interest */
  const Parameters::PostProcessing<dim>::ProbingPoints
    probing_points_parameters;

  /// Postprocessing output frequency
  const unsigned int output_frequency;

  /// Path to the output folder where all .dat are saved
  const std::string output_folder;

  /// Postprocessing verbosity
  const Parameters::Verbosity verbosity;

  /** Contains IDs of probe for which the current time has been added to the
   * table to avoid duplication of time entries when a probing point contains
   * multiple variables */
  std::set<unsigned int> probe_tables_with_current_time;

  /// Tables of probe measures
  std::vector<TableHandler> probe_tables;

  /// Precision of table entries
  const unsigned int table_precision = 6;
};

#endif
