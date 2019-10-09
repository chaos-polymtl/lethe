/*
 * Visualization.h
 *
 *  Created on: Oct 1, 2019
 *      Author: shahab
 */

#include <deal.II/particles/particle_handler.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/particles/particle.h>
#include <tuple>
#include <vector>

using namespace dealii;

#ifndef VISUALIZATION_H_
#define VISUALIZATION_H_


class Visualization : public dealii::DataOutInterface<0,3>
{
public:
	Visualization();

	void build_patches(const Particles::ParticleHandler<3,3> &);

	~Visualization();


private:
          /**
           * Implementation of the corresponding function of the base class.
           */
          virtual const std::vector<DataOutBase::Patch<0,3> > &
          get_patches () const;

          /**
           * Implementation of the corresponding function of the base class.
           */
          virtual std::vector< std::string >
          get_dataset_names () const;

#if DEAL_II_VERSION_GTE(9,1,0)
          virtual
          std::vector<
          std::tuple<unsigned int,
              unsigned int,
              std::string,
              DataComponentInterpretation::DataComponentInterpretation> >
              get_nonscalar_data_ranges () const;
#else
          virtual
          std::vector<std::tuple<unsigned int, unsigned int, std::string> >
          get_vector_data_ranges() const;
#endif

          /**
           * Output information that is filled by build_patches() and
           * written by the write function of the base class.
           */
          std::vector<DataOutBase::Patch<0,3> > patches;

          /**
           * A list of field names for all data components stored in patches.
           */
          std::vector<std::string> dataset_names;

          /**
           * Store which of the data fields are vectors.
           */
#if DEAL_II_VERSION_GTE(9,1,0)
          std::vector<
          std::tuple<unsigned int,
              unsigned int,
              std::string,
              DataComponentInterpretation::DataComponentInterpretation> >
              vector_datasets;
#else
          std::vector<std::tuple<unsigned int, unsigned int, std::string> >
          vector_datasets;
#endif


};
#endif
