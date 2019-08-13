
#ifndef LETHE_PSDHANDLING_H
#define LETHE_PSDHANDLING_H

#include "parameters.h"

using namespace dealii;

class PVDHandler {
public:
  std::vector<std::pair<double, std::string>> times_and_names_;
  void save(std::string filename);
  void read(std::string filename);
  void append(double time, std::string pvtu_filename) {
    times_and_names_.push_back(
        std::pair<double, std::string>(time, pvtu_filename));
  }
  unsigned size() { return times_and_names_.size(); }
};

#endif
