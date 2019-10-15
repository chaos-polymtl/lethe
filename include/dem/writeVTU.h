/*
 * WriteVTU.h
 *
 *  Created on: Oct 1, 2019
 *      Author: shahab
 */
#include "visualization.h"

#ifndef CMAKEFILES_WRITEVTU_H_
#define CMAKEFILES_WRITEVTU_H_

class WriteVTU {
public:
	WriteVTU();
	void write_master_files (const Visualization &data_out, const std::string &solution_file_prefix, const std::vector<std::string> &filenames);
	void writeVTUFiles (int, float, Visualization data_out);
};

#endif /* CMAKEFILES_WRITEVTU_H_ */
