#include <deal.II/base/point.h>

#include <rpt/rpt.h>
#include <rpt/rpt_utilities.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>

template <int dim>
RPT<dim>::RPT(RPTCalculatingParameters &RPTparameters)
  : rpt_parameters(RPTparameters)
{
  // Seed the random number generator
  srand(rpt_parameters.rpt_param.seed);
}

template <int dim>
void
RPT<dim>::setup_and_calculate()
{
  // Reading and assigning positions of detectors and particles
  detectors = assign_detector_positions<dim>(rpt_parameters.detector_param);
  particle_positions = assign_particle_positions<dim>(
    rpt_parameters.rpt_param.particle_positions_file);

  // Calculate count for every particle-detector pair and transfer it to
  // calculated_counts
  calculate_counts();

  // Export results in .csv if enable
  if (rpt_parameters.rpt_param.export_counts)
    export_data();

  // Calculate the cost function for parameters tuning
  if (rpt_parameters.tuning_param.tuning)
    {
      std::vector<double> measured_counts =
        read_counts<dim>(rpt_parameters.tuning_param.experimental_file);
      AssertThrow(
        measured_counts.size() == particle_positions.size(),
        ExcMessage(
          "Prior tuning parameters, the number of particle positions provided"
          " has to be the same number of counts of experimental data. "
          "Note : The experimental counts also have to be at the same positions which can not be verified."))

        double cost_function =
          calculate_cost_function(measured_counts, calculated_counts);
      std::cout << cost_function << std::endl;
    }
}

template <int dim>
void
RPT<dim>::calculate_counts()
{
  // Calculate count for every particle-detector pair
  for (unsigned int i_particle = 0; i_particle < particle_positions.size();
       i_particle++)
    {
      for (unsigned int i_detector = 0; i_detector < detectors.size();
           i_detector++)
        {
          // Create the particle-detector interaction object
          ParticleDetectorInteractions<dim> particle_detector_interactions(
            particle_positions[i_particle],
            detectors[i_detector],

            rpt_parameters.rpt_param);

          double count = particle_detector_interactions.calculate_count();
          calculated_counts.push_back(count);

          if (rpt_parameters.rpt_param.verbosity ==
              Parameters::Verbosity::verbose)
            std::cout << "Count for particle position " << i_particle
                      << " and detector " << i_detector << " : " << count
                      << std::endl;
        }
    }
}


template <int dim>
void
RPT<dim>::export_data()
{
  // Open a file
  std::ofstream myfile;
  std::string   sep;

  std::string filename = rpt_parameters.rpt_param.export_counts_file;

  // check if an extension is specified in the filename, if not add ".csv"
  std::size_t csv_file = filename.find(".csv");
  std::size_t dat_file = filename.find("dat");
  if ((csv_file == std::string::npos) && (dat_file == std::string::npos))
    filename += ".csv";

  myfile.open(filename);
  if (filename.substr(filename.find_last_of(".") + 1) == ".dat")
    {
      myfile
        << "particle_positions_x particle_positions_y particle_positions_z detector_id counts"
        << std::endl;
      sep = " ";
    }
  else // .csv is default
    {
      myfile
        << "particle_positions_x,particle_positions_y,particle_positions_z,detector_id,counts"
        << std::endl;
      sep = ",";
    }

  for (unsigned int i_particle = 0; i_particle < particle_positions.size();
       i_particle++)
    {
      for (unsigned int i_detector = 0; i_detector < detectors.size();
           i_detector++)
        {
          myfile
            << particle_positions[i_particle].get_position()[0] << sep
            << particle_positions[i_particle].get_position()[1] << sep
            << particle_positions[i_particle].get_position()[2] << sep
            << detectors[i_detector].get_id() << sep
            << calculated_counts[i_particle * detectors.size() + i_detector]
            << std::endl;
        }
    }
  myfile.close();
}

template <int dim>
double
RPT<dim>::calculate_cost_function(std::vector<double> &measured_counts,
                                  std::vector<double> &calculated_counts)
{
  double cost_function = 0;

  if (rpt_parameters.tuning_param.cost_function_type ==
      Parameters::RPTTuningParameters::CostFunctionType::larachi)
    {
      for (unsigned int i = 0; i < measured_counts.size(); i++)
        {
          cost_function +=
            std::pow((calculated_counts[i] - measured_counts[i]) /
                       (calculated_counts[i] + measured_counts[i]),
                     2);
        }
    }
  else if (rpt_parameters.tuning_param.cost_function_type ==
           Parameters::RPTTuningParameters::CostFunctionType::l1)
    {
      for (unsigned int i = 0; i < measured_counts.size(); i++)
        cost_function += abs(calculated_counts[i] - measured_counts[i]);

      cost_function /= measured_counts.size();
    }
  else if (rpt_parameters.tuning_param.cost_function_type ==
           Parameters::RPTTuningParameters::CostFunctionType::l2)
    {
      for (unsigned int i = 0; i < measured_counts.size(); i++)
        cost_function +=
          std::pow((calculated_counts[i] - measured_counts[i]), 2);

      cost_function /= measured_counts.size();
    }

  return cost_function;
}

template class RPT<3>;
