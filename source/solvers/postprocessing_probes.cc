// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/postprocessing_probes.h>

template <int dim>
PostprocessingProbes<dim>::PostprocessingProbes(
  std::shared_ptr<SimulationControl> simulation_control,
  const Parameters::PostProcessing<dim>::ProbingPoints
                              &probing_points_parameters,
  const double                 output_frequency,
  const std::string           &output_folder,
  const Parameters::Verbosity &verbosity)
  : simulation_control(simulation_control)
  , probing_points_parameters(probing_points_parameters)
  , output_frequency(output_frequency)
  , output_folder(output_folder)
  , verbosity(verbosity)
  , probe_tables(probing_points_parameters.number_of_probing_points)
{
  prepare_probe_tables();
};

template <int dim>
template <typename VectorType>
void
PostprocessingProbes<dim>::postprocess_probes(
  const Triangulation<dim> &triangulation,
  const Mapping<dim>       &mapping,
  const DoFHandler<dim>    &dof_handler,
  const VectorType         &present_solution,
  const Variable            variable,
  const ConditionalOStream &pcout,
  std::vector<double>      &evaluated_scalar_values)
{
  const auto &probing_points_of_the_variable =
    probing_points_parameters.probing_points_per_variable.at(variable);

  // Some scalar field solutions are not contiguous with the solution vector.
  // We identify here the correct component of the solution vector.
  const unsigned int first_selected_component =
    (variable == Variable::pressure)                         ? dim :
    (variable == Variable::chemical_potential_cahn_hilliard) ? 1 :
                                                               0;

  evaluate_values_at_points(triangulation,
                            mapping,
                            dof_handler,
                            present_solution,
                            probing_points_of_the_variable.points,
                            evaluated_scalar_values,
                            first_selected_component);

  // Get variable strings
  const std::string variable_name = get_variable_string(variable);
  const std::string column_name =
    get_variable_string_with_underscores(variable);

  const double current_time = simulation_control->get_current_time();

  if (verbosity == Parameters::Verbosity::verbose)
    announce_string(pcout, variable_name + " probes");

  // For each probing point fill the correct table
  for (unsigned int i = 0; i < probing_points_of_the_variable.ids.size(); ++i)
    {
      const unsigned int &id = probing_points_of_the_variable.ids[i];
      const std::string  &table_name =
        probing_points_parameters.probing_points_output_names[id];
      const double &evaluated_value = evaluated_scalar_values[i];

      // Console output
      if (verbosity == Parameters::Verbosity::verbose)
        {
          pcout << "At probe " << id << " ("
                << Patterns::Tools::Convert<Point<dim>>::to_string(
                     probing_points_of_the_variable.points[i])
                << ") - " << table_name << std::endl;
          pcout << " Evaluated " << variable_name
                << " value: " << evaluated_value << std::endl;
        }

      TableHandler &current_table = probe_tables[i];

      // Add entries to table
      if (!check_if_table_contains_current_time(id))
        current_table.add_value("time", current_time);
      current_table.add_value(column_name, evaluated_value);
    }
}

template void
PostprocessingProbes<2>::postprocess_probes(
  const Triangulation<2>   &triangulation,
  const Mapping<2>         &mapping,
  const DoFHandler<2>      &dof_handler,
  const GlobalVectorType   &present_solution,
  const Variable            variable,
  const ConditionalOStream &pcout,
  std::vector<double>      &evaluated_scalar_values);

template void
PostprocessingProbes<3>::postprocess_probes(
  const Triangulation<3>   &triangulation,
  const Mapping<3>         &mapping,
  const DoFHandler<3>      &dof_handler,
  const GlobalVectorType   &present_solution,
  const Variable            variable,
  const ConditionalOStream &pcout,
  std::vector<double>      &evaluated_scalar_values);

template void
PostprocessingProbes<2>::postprocess_probes(
  const Triangulation<2>      &triangulation,
  const Mapping<2>            &mapping,
  const DoFHandler<2>         &dof_handler,
  const GlobalBlockVectorType &present_solution,
  const Variable               variable,
  const ConditionalOStream    &pcout,
  std::vector<double>         &evaluated_scalar_values);

template void
PostprocessingProbes<3>::postprocess_probes(
  const Triangulation<3>      &triangulation,
  const Mapping<3>            &mapping,
  const DoFHandler<3>         &dof_handler,
  const GlobalBlockVectorType &present_solution,
  const Variable               variable,
  const ConditionalOStream    &pcout,
  std::vector<double>         &evaluated_scalar_values);

#ifndef LETHE_USE_LDV
template void
PostprocessingProbes<2>::postprocess_probes(
  const Triangulation<2>                           &triangulation,
  const Mapping<2>                                 &mapping,
  const DoFHandler<2>                              &dof_handler,
  const LinearAlgebra::distributed::Vector<double> &present_solution,
  const Variable                                    variable,
  const ConditionalOStream                         &pcout,
  std::vector<double>                              &evaluated_scalar_values);

template void
PostprocessingProbes<3>::postprocess_probes(
  const Triangulation<3>                           &triangulation,
  const Mapping<3>                                 &mapping,
  const DoFHandler<3>                              &dof_handler,
  const LinearAlgebra::distributed::Vector<double> &present_solution,
  const Variable                                    variable,
  const ConditionalOStream                         &pcout,
  std::vector<double>                              &evaluated_scalar_values);
#endif

template <int dim>
template <typename VectorType>
void
PostprocessingProbes<dim>::postprocess_probes(
  const Triangulation<dim>            &triangulation,
  const Mapping<dim>                  &mapping,
  const DoFHandler<dim>               &dof_handler,
  const VectorType                    &present_solution,
  const Variable                       variable,
  const ConditionalOStream            &pcout,
  std::vector<Tensor<1, dim, double>> &evaluated_vector_values)
{
  AssertThrow(variable == Variable::velocity,
              ExcMessage(
                "Only implemented for velocity for now. Implementations for "
                "TimeHarmonicMaxwell variables will be done later."));

  const auto &probing_points_of_the_variable =
    probing_points_parameters.probing_points_per_variable.at(variable);

  const unsigned int first_selected_component = 0;

  evaluate_values_at_points<dim>(triangulation,
                                 mapping,
                                 dof_handler,
                                 present_solution,
                                 probing_points_of_the_variable.points,
                                 evaluated_vector_values,
                                 first_selected_component);

  // Get variable strings
  const std::string variable_name = get_variable_string(variable);
  const std::string column_name_prefix =
    get_variable_string_with_underscores(variable);

  const double current_time = simulation_control->get_current_time();

  if (verbosity == Parameters::Verbosity::verbose)
    announce_string(pcout, variable_name + " probes");

  // For each probing point fill the correct table
  for (unsigned int i = 0; i < probing_points_of_the_variable.ids.size(); ++i)
    {
      const unsigned int &id = probing_points_of_the_variable.ids[i];
      const std::string  &table_name =
        probing_points_parameters.probing_points_output_names[id];
      const Tensor<1, dim, double> &evaluated_value =
        evaluated_vector_values[i];

      // Console output
      if (verbosity == Parameters::Verbosity::verbose)
        {
          pcout << "At probe " << id << " ("
                << Patterns::Tools::Convert<Point<dim>>::to_string(
                     probing_points_of_the_variable.points[i])
                << ") - " << table_name << std::endl;
          pcout << " Evaluated " << variable_name << " vector: "
                << Patterns::Tools::Convert<Tensor<1, dim, double>>::to_string(
                     evaluated_value)
                << std::endl;
        }


      TableHandler &current_table = probe_tables[i];

      // Add entries to table
      if (!check_if_table_contains_current_time(id))
        current_table.add_value("time", current_time);
      current_table.add_value(column_name_prefix + "_x", evaluated_value[0]);
      current_table.add_value(column_name_prefix + "_y", evaluated_value[1]);
      if constexpr (dim == 3)
        current_table.add_value(column_name_prefix + "_z", evaluated_value[2]);
    }
}

template void
PostprocessingProbes<2>::postprocess_probes(
  const Triangulation<2>            &triangulation,
  const Mapping<2>                  &mapping,
  const DoFHandler<2>               &dof_handler,
  const GlobalVectorType            &present_solution,
  const Variable                     variable,
  const ConditionalOStream          &pcout,
  std::vector<Tensor<1, 2, double>> &evaluated_vector_values);

template void
PostprocessingProbes<3>::postprocess_probes(
  const Triangulation<3>            &triangulation,
  const Mapping<3>                  &mapping,
  const DoFHandler<3>               &dof_handler,
  const GlobalVectorType            &present_solution,
  const Variable                     variable,
  const ConditionalOStream          &pcout,
  std::vector<Tensor<1, 3, double>> &evaluated_vector_values);

template void
PostprocessingProbes<2>::postprocess_probes(
  const Triangulation<2>            &triangulation,
  const Mapping<2>                  &mapping,
  const DoFHandler<2>               &dof_handler,
  const GlobalBlockVectorType       &present_solution,
  const Variable                     variable,
  const ConditionalOStream          &pcout,
  std::vector<Tensor<1, 2, double>> &evaluated_vector_values);

template void
PostprocessingProbes<3>::postprocess_probes(
  const Triangulation<3>            &triangulation,
  const Mapping<3>                  &mapping,
  const DoFHandler<3>               &dof_handler,
  const GlobalBlockVectorType       &present_solution,
  const Variable                     variable,
  const ConditionalOStream          &pcout,
  std::vector<Tensor<1, 3, double>> &evaluated_vector_values);

#ifndef LETHE_USE_LDV
template void
PostprocessingProbes<2>::postprocess_probes(
  const Triangulation<2>                           &triangulation,
  const Mapping<2>                                 &mapping,
  const DoFHandler<2>                              &dof_handler,
  const LinearAlgebra::distributed::Vector<double> &present_solution,
  const Variable                                    variable,
  const ConditionalOStream                         &pcout,
  std::vector<Tensor<1, 2, double>>                &evaluated_vector_values);

template void
PostprocessingProbes<3>::postprocess_probes(
  const Triangulation<3>                           &triangulation,
  const Mapping<3>                                 &mapping,
  const DoFHandler<3>                              &dof_handler,
  const LinearAlgebra::distributed::Vector<double> &present_solution,
  const Variable                                    variable,
  const ConditionalOStream                         &pcout,
  std::vector<Tensor<1, 3, double>>                &evaluated_vector_values);
#endif

template <int dim>
std::vector<OutputStructTableHandler>
PostprocessingProbes<dim>::gather_tables()
{
  std::vector<OutputStructTableHandler> table_output_structs;

  std::string prefix = output_folder;
  std::string suffix = ".checkpoint";

  for (unsigned int i = 0; i < probe_tables.size(); ++i)
    {
      const auto &filename =
        probing_points_parameters.probing_points_output_names[i];

      table_output_structs.emplace_back(probe_tables[i],
                                        prefix + filename + suffix);
    }

  return table_output_structs;
}

template class PostprocessingProbes<2>;
template class PostprocessingProbes<3>;
