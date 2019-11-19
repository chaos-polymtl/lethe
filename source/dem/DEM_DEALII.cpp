

#include <deal.II/particles/particle.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/particles/particle_handler.h>
//#include "tests.h"
#include <deal.II/base/array_view.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/base/data_out_base.h>

#include <deal.II/base/partitioner.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <tuple>
#include <algorithm>
#include "demIterator.h"
#include "particleInsertion.h"
#include "readInputScript.h"
#include "visualization.h"
#include "writeVTU.h"
#include "contactSearch.h"
//#include "integration.h"



using namespace dealii;


int main(int argc, char *argv[]) {

	//considering id, type, diameter, density, x, y, z, vx, vy, vz, ax, ay, az, fx, fy, fz, wx, wy, wz, each particle has 19 properties and 9 fields.
    const unsigned int n_properties = 19;
    const unsigned int n_fileds = 8;
    std::vector<std::tuple<std::string,int>> properties(n_properties);
    properties[0]  = std::make_tuple("id",1);   properties[1]  = std::make_tuple("type",1); properties[2]  = std::make_tuple("Diam",1);	properties[3]  = std::make_tuple("Dens",1);
    properties[4]  = std::make_tuple("Pos",3);  properties[5]  = std::make_tuple("Pos",1);  properties[6]  = std::make_tuple("Pos",1);
    properties[7]  = std::make_tuple("Vel",3);  properties[8]  = std::make_tuple("Vel",1);  properties[9]  = std::make_tuple("Vel",1);
    properties[10] = std::make_tuple("a",3); 	properties[11] = std::make_tuple("a",1);    properties[12] = std::make_tuple("a",1);
    properties[13] = std::make_tuple("f",3);    properties[14] = std::make_tuple("f",1);    properties[15] = std::make_tuple("f",1);
    properties[16] = std::make_tuple("w",3); 	properties[17] = std::make_tuple("w",1); 	properties[18] = std::make_tuple("w",1);

	Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
	//MPILogInitAll all; //question

	//total number of particles in the system
	int nPart = 0;
	//DEM clock
	int DEM_step = 0;
	float DEM_time = 0;

	ReadInputScript readInput;
	ParticleInsertion ins1(readInput);

	parallel::distributed::Triangulation<3,3> tr(MPI_COMM_WORLD);
	Particles::Particle<3> particle;
	GridGenerator::hyper_cube(tr,-0.1,0.1);
	int numRef = 3;
	tr.refine_global(numRef);
	//MappingQ<3,3> mapping(3);
	  	  	  	  	  	  	  	  auto &mappinggg = StaticMappingQ1<3,3>::mapping;

	Particles::ParticleHandler<3,3> particle_handler(tr, mappinggg , n_properties);
	Particles::PropertyPool pool(n_properties);


	int cellNum = pow(pow(2,numRef),3);
	std::pair <std::vector<std::set<Triangulation<3>::active_cell_iterator>>,std::vector<Triangulation<3>::active_cell_iterator>> cellNeighbor;
	ContactSearch cs1;
	cellNeighbor=cs1.findCellNeighbors(cellNum, tr);


	DEM_iterator iter1;
	std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>,Particles::ParticleIterator<3,3>>, std::vector<double>, double, std::vector<double>, double, std::vector<double>, std::vector<double>, double, double >> contactInfo;
	//Integration	 Integ1;

//////////////////////////////add out of sim box -> delete particle
	//Insertion phase:
	while (DEM_step < readInput.tInsertion)
	{

		if (nPart < readInput.nTotal) //number < total number
			{
			if (fmod(DEM_step,readInput.insertFrequncy) == 1)
				{
					ins1.uniformInsertion(particle_handler, tr, readInput, nPart, pool, particle);
				}
			}

		if (fmod(DEM_step,readInput.writeFrequency) == 1)
			{
				Visualization visObj;
				visObj.build_patches(particle_handler, n_fileds, n_properties, properties);
				std::string particle_file_prefix = "Globals";
				std::vector<std::string> filenames;
				filenames.push_back (particle_file_prefix + ".vtu");
				WriteVTU writObj;
				writObj.write_master_files(visObj, particle_file_prefix, filenames);

				std::string filename = (("particles/Out_" +  Utilities::int_to_string (DEM_step, 4)));
				std::ofstream output((filename + ".vtu"));
			    DataOutBase::VtkFlags vtk_flags;
			    vtk_flags.cycle = DEM_step;
			    vtk_flags.time = DEM_time;
			    visObj.set_flags (vtk_flags);
			    visObj.write_vtu(output);
			}

		iter1.engine(nPart, particle_handler, tr, DEM_step, DEM_time, readInput, cellNeighbor, contactInfo);
		std::cout<<DEM_step<<std::endl;

	}

	//Operation phase:std::vector<Point<3>> points;
	while (DEM_step < readInput.tFinal)
	{
		if (fmod(DEM_step,readInput.writeFrequency) == 1)
			{
				Visualization visObj;
				visObj.build_patches(particle_handler, n_fileds, n_properties, properties);
				std::string filename = (("particles/Out_" +  Utilities::int_to_string (DEM_step, 4)));
				std::ofstream output((filename + ".vtu"));
			    DataOutBase::VtkFlags vtk_flags;
			    vtk_flags.cycle = DEM_step;
			    vtk_flags.time = DEM_time;
			    visObj.set_flags (vtk_flags);
			    visObj.write_vtu(output);
			}

		iter1.engine(nPart, particle_handler, tr, DEM_step, DEM_time, readInput, cellNeighbor, contactInfo);
		std::cout<<DEM_step<<std::endl;
	}


//if tf<tins add error
	//yek bar call kardane neighbor cell baraye kole sim kafie, dorostesh kon
	//checke nahayi baraye dorost kar kardane hame chi, aval integration ke az inja call kardanesh hazf shod, badesh contact det





	return 0;
}
