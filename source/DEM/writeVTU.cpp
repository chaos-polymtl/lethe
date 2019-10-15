/*
 * WriteVTU.cpp
 *
 *  Created on: Oct 1, 2019
 *      Author: shahab
 */

#include "writeVTU.h"

#include <deal.II/base/data_out_base.h>
#include <fstream>
#include "visualization.h"


using namespace dealii;

WriteVTU::WriteVTU() {
	// TODO Auto-generated constructor stub

}

//this function writes the pvtu file:
void WriteVTU::write_master_files (const Visualization &data_out, const std::string &solution_file_prefix, const std::vector<std::string> &filenames)
    {

      const std::string pvtu_master_filename = (solution_file_prefix + ".pvtu");
      std::ofstream pvtu_master (( "particles/" + pvtu_master_filename).c_str());

      data_out.write_pvtu_record (pvtu_master, filenames);
    }

void WriteVTU::writeVTUFiles (int DEM_step, float DEM_time, Visualization data_out)
	{


    //const unsigned int n_processes = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());
	//const unsigned int n_processes = 1;

    std::ofstream filename(("particles/" +  Utilities::int_to_string (DEM_step, 4) + ".vtu"));

    // pass time step number and time as metadata into the output file
    DataOutBase::VtkFlags vtk_flags;
    vtk_flags.cycle = DEM_step;
    vtk_flags.time = DEM_time;

    data_out.set_flags (vtk_flags);

    // There is a if definition for variable group_files in particle.cc, here I used 1 temporarily:
    //int group_files = 1;


    data_out.write_vtu(filename);





	}
