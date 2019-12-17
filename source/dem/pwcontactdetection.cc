/*
 * pwcontactdetection.cpp
 *
 *  Created on: Nov 26, 2019
 *      Author: shahab
 */
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/particle_handler.h>


#include "dem/pwcontactdetection.h"
#include <iostream>

using namespace dealii;

pwcontactdetection::pwcontactdetection() {

}

void pwcontactdetection::boundaryCellsAndFaces(const Triangulation<3,3> &tr, std::vector<std::tuple<int, Triangulation<3>::active_cell_iterator, int, Point<3>, Point<3>>> &boundaryCellInfo)
{
	std::vector<std::pair<int, Triangulation<3>::active_cell_iterator>> searchPair;
	const FE_Q<3,3> fe(1);
	QGauss<3>   quadrature_formula(3);
	QGauss<3-1> face_quadrature_formula(3);
	unsigned int n_face_q_points = face_quadrature_formula.size();

	   FEFaceValues<3> fe_face_values (fe, face_quadrature_formula,  update_values | update_quadrature_points |  update_normal_vectors);
		for (Triangulation<3>::active_cell_iterator cell = tr.begin_active(); cell != tr.end(); ++cell)
		{
			for (unsigned int face_id = 0; face_id < GeometryInfo<3>::faces_per_cell; ++face_id)
			{
		          if (cell->face(face_id)->at_boundary() == true)
		          {
		        	  fe_face_values.reinit(cell, face_id);

		  	          for (unsigned int f_q_point = 0; f_q_point <n_face_q_points; ++f_q_point)
		  	                {

		  	                  //const Tensor<1, 3> &Nv =  fe_face_values.normal_vector(f_q_point);
		  	                Point<3> surfNormal; //question? why doesnt it work in one line?!
		  	                surfNormal = fe_face_values.normal_vector(f_q_point);
		  	                surfNormal = -1 * surfNormal;
		  	                Point<3> quadPoint = fe_face_values.quadrature_point(3);
		  	                int boundID = cell->face(face_id)->boundary_id();
		  	                std::tuple cellInfo = std::make_tuple(boundID, cell, face_id, surfNormal, quadPoint);
		  	                std::pair cellInfoPair = std::make_pair(boundID, cell);

							auto it = std::find(searchPair.begin(), searchPair.end(), cellInfoPair);
							if(it == searchPair.end())
							{
		  	                boundaryCellInfo.push_back(cellInfo);
							}

							searchPair.push_back(std::make_pair(boundID, cell));
		  	                }

		          }
			}
		}
}

std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>, int>, Point<3>, Point<3>>> pwcontactdetection::pwcontactlist(std::vector<std::tuple<int, Triangulation<3>::active_cell_iterator, int, Point<3>, Point<3>>> boundaryCellInfo, Particles::ParticleHandler<3,3> &particle_handler)
{
	std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>, int>, Point<3>, Point<3>>> pwContactList;
	std::tuple<std::pair<Particles::ParticleIterator<3,3>, int>, Point<3>, Point<3>> pwInfoTuple;
	std::vector<std::pair<Particles::ParticleIterator<3,3>, int>> searchPair;


	for (unsigned int i = 0; i < boundaryCellInfo.size(); ++i)
	{
		Triangulation<3>::active_cell_iterator workingCell = std::get<1>(boundaryCellInfo[i]);
		const Particles::ParticleHandler<3,3>::particle_iterator_range particle_range = particle_handler.particles_in_cell(workingCell);


		for (typename Particles::ParticleHandler<3,3>::particle_iterator_range::iterator partIter = particle_range.begin(); partIter != particle_range.end();++partIter)
		{
			std::pair pwPair = std::make_pair(partIter, std::get<0>(boundaryCellInfo[i]));

			auto it4 = std::find(searchPair.begin(), searchPair.end(), pwPair);
			if(it4 == searchPair.end())
			{
				pwInfoTuple = std::make_tuple(pwPair, std::get<3>(boundaryCellInfo[i]), std::get<4>(boundaryCellInfo[i]));
				pwContactList.push_back(pwInfoTuple);
				searchPair.push_back(pwPair);
			}
		}
	}

return pwContactList;
}




void pwcontactdetection::pwFineSearch(std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>, int> , Point<3>, Point<3>>> pwContactList , Particles::ParticleHandler<3,3> &particle_handler , std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>, int>, Point<3>, Point<3>, double, double, double, Point<3>, double >> &pwContactInfo, float dt)
{

	std::tuple<std::pair<Particles::ParticleIterator<3,3>, int>, Point<3>, Point<3>, double, double, double, Point<3>, double> pwInfoTuple;
	std::vector<std::pair<Particles::ParticleIterator<3,3>, int>> pwSearchPair;
	if (!pwSearchPair.empty())
	{
		pwSearchPair.clear();
	}
//reread this part:
	for (unsigned int i = 0; i < pwContactInfo.size(); i++)
	{
		Point<3> partCenter = std::get<0>(pwContactInfo[i]).first->get_location();
		Point<3> partVector;
		partVector = partCenter - std::get<2>(pwContactInfo[i]);
		Point<3> projection = findProjection(partVector , std::get<1>(pwContactInfo[i]));

		double distance = ((std::get<0>(pwContactInfo[i]).first->get_properties()[2])/2) - (sqrt(projection.square())) ;
		if (distance > 0)
		{

			Point<3> normVec = std::get<1>(pwContactInfo[i]);
			Point<3> pVel = {std::get<0>(pwContactInfo[i]).first->get_properties()[7], std::get<0>(pwContactInfo[i]).first->get_properties()[8], std::get<0>(pwContactInfo[i]).first->get_properties()[9]};
			Point<3> pOmega = {std::get<0>(pwContactInfo[i]).first->get_properties()[16], std::get<0>(pwContactInfo[i]).first->get_properties()[17], std::get<0>(pwContactInfo[i]).first->get_properties()[18]};
			Point<3> relVel = pVel + crossProduct((((std::get<0>(pwContactInfo[i]).first->get_properties()[2])/2) * pOmega) , normVec);
			double normRelVel = dotProduct(relVel,normVec);
			Point<3> relNormVel = normRelVel * normVec;
			Point<3> relTangVel;
			relTangVel = relVel - relNormVel;
			Point<3> tangVec = {0, 0, 0};
			double relTangVelVal = vecValue(relTangVel);
			if(relTangVelVal != 0)
			{
				tangVec = relTangVel/relTangVelVal;
			}
			double tangRelVel = dotProduct(relVel , tangVec);
			double tangOverlap = std::get<5>(pwContactInfo[i]) + (tangRelVel * dt);

			pwInfoTuple = std::make_tuple(std::get<0>(pwContactInfo[i]), normVec, std::get<2>(pwContactInfo[i]), distance, normRelVel, tangOverlap, tangVec, tangRelVel);
			pwContactInfo[i] = pwInfoTuple;


		}

		else
		{

			pwContactInfo.erase(pwContactInfo.begin()+i);

		}
		// find how to perform vec = std::get<0>(anothervec)
		pwSearchPair.push_back(std::get<0>(pwContactInfo[i]));
	}


	for (unsigned int i = 0; i < pwContactList.size(); i++)
		{

		Point<3> partCenter = std::get<0>(pwContactList[i]).first->get_location();
		Point<3> partVector;
		partVector = partCenter - std::get<2>(pwContactList[i]);
		Point<3> projection = findProjection(partVector , std::get<1>(pwContactList[i]));
		double distance = ((std::get<0>(pwContactList[i]).first->get_properties()[2])/2) - (sqrt(projection.square())) ;

		auto it4 = std::find(pwSearchPair.begin(), pwSearchPair.end(), std::get<0>(pwContactList[i]));
		if(it4 == pwSearchPair.end())
		{
		if(distance > 0)
			{



			Point<3> normVec = std::get<1>(pwContactList[i]);


			Point<3> pVel = {std::get<0>(pwContactList[i]).first->get_properties()[7], std::get<0>(pwContactList[i]).first->get_properties()[8], std::get<0>(pwContactList[i]).first->get_properties()[9]};
			Point<3> pOmega = {std::get<0>(pwContactList[i]).first->get_properties()[16], std::get<0>(pwContactList[i]).first->get_properties()[17], std::get<0>(pwContactList[i]).first->get_properties()[18]};
			Point<3> relVel = pVel + crossProduct((((std::get<0>(pwContactList[i]).first->get_properties()[2])/2) * pOmega) , normVec);
			double normRelVel = dotProduct(relVel,normVec);
			Point<3> relNormVel = normRelVel * normVec;
			Point<3> relTangVel;
			relTangVel = relVel - relNormVel;
			Point<3> tangVec = {0, 0, 0};



			double relTangVelVal = vecValue(relTangVel);
			if(relTangVelVal != 0)
			{
				tangVec = relTangVel/relTangVelVal;
			}


			double tangRelVel = dotProduct(relVel , tangVec);
			double tangOverlap = 0;

																  //   << " " ;

			pwInfoTuple = std::make_tuple(std::get<0>(pwContactList[i]), normVec, std::get<2>(pwContactList[i]), distance, normRelVel, tangOverlap, tangVec, tangRelVel);

			pwContactInfo.push_back(pwInfoTuple);




			}

		}

		}


}





Point<3> pwcontactdetection::findProjection(Point<3> pointA, Point<3> pointB)
{
	Point<3> pointC;
    pointC = ((dotProduct(pointA,pointB)) / (pointB.square())) * pointB;

    return pointC;
	}


double pwcontactdetection::dotProduct(Point<3> A, Point<3> B)
{
	return (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
	}


Point<3> pwcontactdetection::crossProduct(Point<3> A, Point<3> B)
		{
			return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
		}



double pwcontactdetection::vecValue(Point<3> A)
{
	return (sqrt(pow(A[0],2)+pow(A[1],2)+pow(A[2],2)));
	}










