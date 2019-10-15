/*
 * ParticleInsertion.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: meteor
 */

#include "particleInsertion.h"

#include "iostream"
#include <fstream>
#include <string>
			//#include "WriteOutput.h"

#include <deal.II/distributed/tria.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/base/array_view.h>
#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/particles/property_pool.h>
#include "readInputScript.h"

using namespace dealii;

ParticleInsertion::ParticleInsertion(ReadInputScript readInput)
{
	int n_exp;

	n_exp = int((readInput.ins_x_max-readInput.ins_x_min)/(2*readInput.diameter)) * int((readInput.ins_y_max-readInput.ins_y_min)/(2*readInput.diameter)) * int((readInput.ins_z_max-readInput.ins_z_min)/(2*readInput.diameter));
	if (readInput.nInsert > n_exp)
		std::cout<<"The inserted number of particles (" << readInput.nInsert <<  ") is higher than maximum expected number of particles (" << n_exp << ")" << std::endl;
} // add error here


void ParticleInsertion::uniformInsertion(Particles::ParticleHandler<3,3> &particle_handler,
		  const Triangulation<3,3> &tr, ReadInputScript readInput, int &nPart, Particles::PropertyPool &pool, Particles::Particle<3> &particle)
{
	///typename Particles::PropertyPool::Handle handle = pool.allocate_properties_array();
	int nx = int((readInput.ins_x_max-readInput.ins_x_min) / (2*readInput.diameter));
	int ny = int((readInput.ins_y_max-readInput.ins_y_min) / (2*readInput.diameter));
	int nz = int((readInput.ins_z_max-readInput.ins_z_min) / (2*readInput.diameter));
	int nP = 0;


		for (int i=0; i < nx; ++i)
			for (int j=0; j < ny; ++j)
				for (int k=0; k < nz; ++k)
					if(nP < readInput.nInsert)
						{
						Point<3> 		position;
						Point<3>    	reference_position;
						unsigned int    id;

						position[0] = readInput.ins_x_min + (readInput.diameter/2) + (i * 2 * readInput.diameter);
						position[1] = readInput.ins_y_min + (readInput.diameter/2) + (j * 2 * readInput.diameter);
						position[2] = readInput.ins_z_min + (readInput.diameter/2) + (k * 2 * readInput.diameter);
						id = i * ny * nz + j * nz + k + nPart + 1;

						Particles::Particle<3> particle(position, reference_position, id);
						Triangulation<3,3>:: active_cell_iterator cell = GridTools::find_active_cell_around_point(tr, particle.get_location());

						if (id == 95)
						{

							std::cout<< GridTools::find_active_cell_around_point(tr, particle.get_location())->id() << std::endl;
						}

						Particles::ParticleIterator<3,3> pit = particle_handler.insert_particle(particle, cell);

						particle.set_property_pool(pool);

						pit->get_properties()[0] = id;
						pit->get_properties()[1] = 1;
						pit->get_properties()[2] = readInput.diameter;
						pit->get_properties()[3] = readInput.density;
						//Position
						pit->get_properties()[4] = position[0];
						pit->get_properties()[5] = position[1];
						pit->get_properties()[6] = position[2];
						//Velocity
						pit->get_properties()[7] = 0;
						pit->get_properties()[8] = 0;
						pit->get_properties()[9] = 0;
						//Acceleration
						pit->get_properties()[10] = 0 + readInput.g[0];
						pit->get_properties()[11] = 0 + readInput.g[1];
						pit->get_properties()[12] = 0 + readInput.g[2];
						//Force
						pit->get_properties()[13] = 0;
						pit->get_properties()[14] = 0;
						pit->get_properties()[15] = 0;
						//w
						pit->get_properties()[16] = 0;
						pit->get_properties()[17] = 0;
						pit->get_properties()[18] = 0;

						++nP;
						}

	nPart = nPart + readInput.nInsert;
}
