/*
 * WriteVTU.cpp
 *
 *  Created on: Oct 1, 2019
 *      Author: shahab
 */

#include <dem/write_vtu.h>

using namespace dealii;

template <int dim> WriteVTU<dim>::WriteVTU() {}

template <int dim>
void WriteVTU<dim>::write_master_files(
    Visualization<dim> &visulization_object,
    const DEMSolverParameters<dim> &dem_parameters) {
  // Defining properties as local variable
  const auto visualization_info = dem_parameters.outputProperties;

  // Getting the general names of simulation results
  std::string solution_file_prefix = visualization_info.general_file_prefix;
  std::vector<std::string> filenames;
  filenames.push_back(solution_file_prefix + ".vtu");
  const std::string pvtu_master_filename = (solution_file_prefix + ".pvtu");

  // Getting the output folder name creating it
  const auto output_folder = "./" + visualization_info.output_folder;
  mkdir(output_folder.c_str(), 0777);
  const auto directory_path = output_folder + "/";

  // Writing .pvtu file
  std::ofstream pvtu_master((directory_path + pvtu_master_filename).c_str());
  visulization_object.write_pvtu_record(pvtu_master, filenames);
}

template <int dim>
void WriteVTU<dim>::writeVTUFiles(
    Visualization<dim> &visulization_object, const int &DEM_step,
    const double &DEM_time, const DEMSolverParameters<dim> &dem_parameters) {
  // const unsigned int n_processes =
  // Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());
  // const unsigned int n_processes = 1;

  // Defining properties as local variable
  const auto visualization_info = dem_parameters.outputProperties;

  // Defining the file name
  const auto output_folder = "./" + visualization_info.output_folder;
  const auto output_full_name =
      output_folder + "/" + visualization_info.result_prefix + "_";
  std::string filename =
      ((output_full_name.c_str() + Utilities::int_to_string(DEM_step, 4)));

  // Writing the .vtu output files
  std::ofstream output((filename + ".vtu"));
  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.cycle = DEM_step;
  vtk_flags.time = DEM_time;
  visulization_object.set_flags(vtk_flags);
  visulization_object.write_vtu(output);
}

template class WriteVTU<2>;
template class WriteVTU<3>;
