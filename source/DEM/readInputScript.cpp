/*
 * ReadInputScript.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: meteor
 */

#include "readInputScript.h"

#include <iostream>
#include <fstream>

using namespace std;

ReadInputScript::ReadInputScript() {
	int n = 3;
	int m = 4;

	ifstream input("./InputScript.txt");

	if (input.is_open())
	{
		while (!input.eof())
		{
			for (int j = 0; j < n-1; ++j) {
			    input.ignore(1000, '\n');
			  }
			  input >> numOfParams >> diameter >> density >> CoR;
			for (int j = 0; j < m-1; ++j) {
			  	input.ignore(1000, '\n');
			  }
			  input >> ins_x_min >> ins_x_max >> ins_y_min >> ins_y_max >> ins_z_min >> ins_z_max;
			for (int j = 0; j < m-1; ++j) {
				 input.ignore(1000, '\n');
			  }
			  input >> x_min >> x_max >> y_min >> y_max >> z_min >> z_max;
			for (int j = 0; j < m-1; ++j) {
				 input.ignore(1000, '\n');
			  }
				 input >> dt >> tFinal;
			for (int j = 0; j < m-1; ++j) {
				 input.ignore(1000, '\n');
			  }
				 input >> tInsertion >> nTotal >> nInsert >> insertFrequncy;
			for (int j = 0; j < m-1; ++j) {
				 input.ignore(1000, '\n');
			  }
				 input >> g[0] >> g[1] >> g[2];
			for (int j = 0; j < m-1; ++j) {
				 input.ignore(1000, '\n');
			 }
				 input >> writeFrequency;
		}
// error handling

	}
	else std::cout<<"input script is not open"<<std::endl;

}


