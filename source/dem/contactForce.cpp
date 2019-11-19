/*
 * contactForce.cpp
 *
 *  Created on: Oct 31, 2019
 *      Author: shahab
 */

#include "contactForce.h"
#include "demIterator.h"

ContactForce::ContactForce()
{}


void ContactForce::linearCF(std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>, Particles::ParticleIterator<3, 3>>, std::vector<double>, double, std::vector<double>, double, std::vector<double>, std::vector<double>, double, double>> contactInfo , Particles::ParticleHandler<3,3> &particle_handler, ReadInputScript readInput)
{
	for (unsigned int i = 0; contactInfo.size(); i++)
	{
		{
			std::vector<double> normalForce  = {vecSubtract(numVecProd((-1.0 * readInput.kn * std::get<2>(contactInfo[i])), std::get<3>(contactInfo[i])) , numVecProd((readInput.ethan * std::get<4>(contactInfo[i])) , std::get<3>(contactInfo[i])))};

			std::get<0>(contactInfo[i]).first->get_properties()[13] = normalForce[0];
			std::get<0>(contactInfo[i]).first->get_properties()[14] = normalForce[1];
			std::get<0>(contactInfo[i]).first->get_properties()[15] = normalForce[2];

			std::vector<double> tangForce = { vecSubtract(numVecProd((-1.0 * readInput.kt * std::get<8>(contactInfo[i])) , std::get<6>(contactInfo[i])) , numVecProd((readInput.ethat * std::get<7>(contactInfo[i])) , std::get<6>(contactInfo[i])) ) };
			if(vecValue(tangForce) < (readInput.mu * vecValue(normalForce)))
			{
				std::get<0>(contactInfo[i]).first->get_properties()[13] = std::get<0>(contactInfo[i]).first->get_properties()[13] + tangForce[0];
				std::get<0>(contactInfo[i]).first->get_properties()[14] = std::get<0>(contactInfo[i]).first->get_properties()[14] + tangForce[1];
				std::get<0>(contactInfo[i]).first->get_properties()[15] = std::get<0>(contactInfo[i]).first->get_properties()[15] + tangForce[2];
			}
			else
			{
				std::vector<double> coulumbTangForce = {numVecProd((readInput.mu * vecValue(normalForce) * sgn(std::get<8>(contactInfo[i]))) , std::get<6>(contactInfo[i])) };

				std::get<0>(contactInfo[i]).first->get_properties()[13] = std::get<0>(contactInfo[i]).first->get_properties()[13] + coulumbTangForce[0];
				std::get<0>(contactInfo[i]).first->get_properties()[14] = std::get<0>(contactInfo[i]).first->get_properties()[14] + coulumbTangForce[1];
				std::get<0>(contactInfo[i]).first->get_properties()[15] = std::get<0>(contactInfo[i]).first->get_properties()[15] + coulumbTangForce[2];
			}

			std::get<0>(contactInfo[i]).second->get_properties()[13] = -1.0 * std::get<0>(contactInfo[i]).first->get_properties()[13];
			std::get<0>(contactInfo[i]).second->get_properties()[14] = -1.0 * std::get<0>(contactInfo[i]).first->get_properties()[14];
			std::get<0>(contactInfo[i]).second->get_properties()[15] = -1.0 * std::get<0>(contactInfo[i]).first->get_properties()[15];

			// does it have access to particle handler or I need to loop over particle handler and search for contactInfo particles
	}
}
}



	double ContactForce::dotProduct(std::vector<double> A, std::vector<double> B)
	{
		return (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
		}

	std::vector<double> ContactForce::crossProduct(std::vector<double> A, std::vector<double> B)
			{
				return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
			}

	std::vector<double> ContactForce::vecSubtract(std::vector<double> A, std::vector<double> B)
			{
				return {A[0] - B[0], A[1] - B[1], A[2] - B[2]};
			}

	std::vector<double> ContactForce::vecAdd(std::vector<double> A, std::vector<double> B)
	{
		return {A[0] + B[0], A[1] + B[1], A[2] + B[2]};
		}

	std::vector<double> ContactForce::numVecProd(double A, std::vector<double> B)
			{
		return {A * B[0], A * B[1], A * B[2]};
			}

	double ContactForce::vecValue(std::vector<double> A)
	{
		return (sqrt(pow(A[0],2)+pow(A[1],2)+pow(A[2],2)));
		}

	int ContactForce::sgn(float a)
	{
		int b;
			if (a > 0){
				b = 1;
			}else if (a < 0){
				b = -1;
			}else if (a == 0) {
				b = 0;
			}
	}


