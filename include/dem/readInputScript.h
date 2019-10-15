/*
 * ReadInputScript.h
 *
 *  Created on: Sep 24, 2019
 *      Author: meteor
 */
#ifndef READINPUTSCRIPT_H_
#define READINPUTSCRIPT_H_

class ReadInputScript {
public:
	float diameter, density, CoR;
	float ins_x_min, ins_y_min, ins_z_min, ins_x_max, ins_y_max, ins_z_max;
	float x_min, y_min, z_min, x_max, y_max, z_max;
	float dt;
	int tInsertion, tFinal, nTotal, nInsert, insertFrequncy, writeFrequency;
	int	  numOfParams;
	float g[3];

	ReadInputScript();

};

#endif /* READINPUTSCRIPT_H_ */
