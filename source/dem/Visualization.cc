/*
 * Visualization.cpp
 *
 *  Created on: Oct 1, 2019
 *      Author: shahab
 */

#include "dem/Visualization.h"
#include <vector>
using namespace dealii;

Visualization::Visualization() {


}



void Visualization :: build_patches(const dealii::Particles::ParticleHandler<3,3> &particle_handler)
      {
//link it to the input, so the user defines what parameters to be used
// id, type, size, density, position, velocity, a, f, w
	//int numberofprop = 5;
	int nProp = 5;
	std::vector<std::string> properties(nProp);
	properties[0] = "id"; properties[1] = "type"; properties[2] = "size";
	properties[3] = "Position"; properties[4] = "v"; //properties[5] = "Position";
	//properties[6] = "vel"; properties[7] = "vel"; properties[8] = "vel";
	//properties[7] = "f"; properties[8] = "w";

	std::vector<int> propElements(nProp);
	propElements[0] = 1; propElements[1] = 1; propElements[2] = 1; propElements[3] = 3;
	propElements[4] = 3;// propElements[5] = 3; propElements[6] = 3;
	//propElements[7] = 3; propElements[8] = 3;



	dataset_names.reserve(nProp+1);

    for (int n = 0; n < properties.size(); ++n)
    {
    	dataset_names.push_back(properties[n]);
    }



    //dataset_names.reserve(nProp+1);
   // dataset_names.push_back("id");
    //dataset_names.push_back("pos");


int totalNumofComp = 0;

//numberofprop replace
    for (unsigned int field_index = 0; field_index < properties.size() ; ++field_index)
       {
         const unsigned n_components = propElements[field_index];
         const std::string field_name = properties[field_index];

         // HDF5 only supports 3D vector output, therefore only treat output fields as vector if we
         // have a dimension of 3 and 3 components.
         const bool field_is_vector = (propElements[field_index] < 2)
                                      ?
                                      n_components == propElements[field_index]
                                      :
                                      n_components == 3;


         // If it is a 1D element, or a vector, print just the name, otherwise append the index after an underscore
        if (n_components == 1)
        {
          totalNumofComp = totalNumofComp + n_components;
      //    for (unsigned int component_index=0; component_index<n_components; ++component_index)
          //   dataset_names.push_back(field_name);
        }

         if (n_components == 3)
         {
           totalNumofComp = totalNumofComp + 2 * n_components;
          // for (unsigned int component_index=0; component_index<n_components; ++component_index)
            // dataset_names.push_back(field_name + "_" + Utilities::to_string(component_index));
         }
        }



     // Second store which of these data fields are vectors

     for (unsigned int field_index = 0; field_index < properties.size() ; ++field_index)
       {

         const unsigned n_components = propElements[field_index];

         // If the property has dim components, we treat it as vector
         if (n_components == 3)
           {
        	 //field_position
             const unsigned int field_position = field_index-1;
             const std::string field_name = properties[field_index];


			#if DEAL_II_VERSION_GTE(9,1,0)
							vector_datasets.push_back(std::make_tuple(field_position+1,
																	  field_position+n_components,
																	  field_name,
																	  DataComponentInterpretation::component_is_part_of_vector));
			#else
							vector_datasets.push_back(std::make_tuple(field_position+1,
																	  field_position+n_components,
																	  field_name));
			#endif


           }
       }



        // Third build the actual patch data
        patches.resize(particle_handler.n_locally_owned_particles());

        typename dealii::Particles::ParticleHandler<3,3>::particle_iterator particle = particle_handler.begin();

        for (unsigned int i=0; particle != particle_handler.end(); ++particle, ++i)
          {
            patches[i].vertices[0] = particle->get_location();
            patches[i].patch_index = i;
            patches[i].n_subdivisions = 1;
            //**********************************
            //patches[i].data.reinit(property_information.n_components()+1,1);  //***********************************
            patches[i].data.reinit(totalNumofComp-9,1);

           // patches[i].data(0,0) = particle->get_id();
            if (particle->has_properties())
             {
                const ArrayView<const double> propertiess = particle->get_properties();


                for (unsigned int property_index = 0; property_index < propertiess.size(); ++property_index){
                  patches[i].data(property_index,0) = propertiess[property_index];
                  std::cout<< patches[i].data(property_index,0)<<std::endl;}
              }
          }
      }




   const std::vector<DataOutBase::Patch<0,3> > &
   Visualization::get_patches () const
   {
     return patches;
   }

   std::vector< std::string >
   Visualization::get_dataset_names () const
   {
     return dataset_names;
   }

#if DEAL_II_VERSION_GTE(9,1,0)

   std::vector<
   std::tuple<unsigned int,
       unsigned int, std::string,
       DataComponentInterpretation::DataComponentInterpretation> >
   	   Visualization::get_nonscalar_data_ranges () const
   {
     return vector_datasets;
   }
#else

   std::vector<std::tuple<unsigned int, unsigned int, std::string> >
   Visualization::get_vector_data_ranges () const
   {
     return vector_datasets;
   }
#endif









Visualization::~Visualization() {


}
