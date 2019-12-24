/*
 * WriteVTU.cpp
 *
 *  Created on: Oct 1, 2019
 *      Author: shahab
 */

#include "dem/write_vtu.h"

#include <deal.II/base/data_out_base.h>

#include <fstream>

#include "dem/visualization.h"


using namespace dealii;

WriteVTU::WriteVTU()
{
}

// this function writes the pvtu file:
void
WriteVTU::write_master_files(const Visualization &data_out)
{
    std::string              solution_file_prefix = "Globals";
    std::vector<std::string> filenames;
    filenames.push_back(solution_file_prefix + ".vtu");
  const std::string pvtu_master_filename = (solution_file_prefix + ".pvtu");
  std::ofstream     pvtu_master(("particles/" + pvtu_master_filename).c_str());

  data_out.write_pvtu_record(pvtu_master, filenames);
}

void
WriteVTU::writeVTUFiles(Visualization &visObj, int DEM_step, float DEM_time)
{
  // const unsigned int n_processes =
  // Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());
  // const unsigned int n_processes = 1;
    std::string filename =
      (("particles/Out_" + Utilities::int_to_string(DEM_step, 4)));
    std::ofstream         output((filename + ".vtu"));
    DataOutBase::VtkFlags vtk_flags;
    vtk_flags.cycle = DEM_step;
    vtk_flags.time  = DEM_time;
    visObj.set_flags(vtk_flags);
    visObj.write_vtu(output);

}
