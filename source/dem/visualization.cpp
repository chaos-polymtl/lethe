/*
 * Visualization.cpp
 *
 *  Created on: Oct 1, 2019
 *      Author: shahab
 */

#include "visualization.h"

#include <vector>
using namespace dealii;

Visualization::Visualization() {
}


void Visualization :: build_patches(const dealii::Particles::ParticleHandler<3,3> &particle_handler, const unsigned int n_fileds, const unsigned int n_properties, std::vector<std::tuple<std::string,int>> properties)
      {

	dataset_names.reserve(n_fileds);

    for (unsigned int n = 0; n < properties.size(); ++n)
    {
    	dataset_names.push_back(std::get<0>(properties[n]));
    }

    for (unsigned int field_index = 0; field_index < properties.size() ; ++field_index)
       {
         const std::string field_name = std::get<0>(properties[field_index]);
        }

     // Second store which of these data fields are vectors
     for (unsigned int field_index = 0; field_index < properties.size() ; ++field_index)
       {
         const unsigned n_components = std::get<1>(properties[field_index]);
         if (n_components == 3)
           {
             const unsigned int field_position = field_index;
             const std::string field_name = std::get<0>(properties[field_index]);
			 vector_datasets.push_back(std::make_tuple(field_position, field_position+n_components-1, field_name, DataComponentInterpretation::component_is_part_of_vector));
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
            patches[i].data.reinit(n_properties,1);

            if (particle->has_properties())
             {
                const ArrayView<const double> props = particle->get_properties();

                for (unsigned int property_index = 0; property_index < props.size(); ++property_index)
                  patches[i].data(property_index,0) = props[property_index];
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

   std::vector<std::tuple<unsigned int, unsigned int, std::string, DataComponentInterpretation::DataComponentInterpretation>> Visualization::get_nonscalar_data_ranges () const
   {
     return vector_datasets;
   }


Visualization::~Visualization() {
}
