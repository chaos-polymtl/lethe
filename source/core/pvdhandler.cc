
#include "core/pvdhandler.h"

#include <fstream>

using namespace dealii;

void
PVDHandler::save(std::string prefix)
{
  std::string   filename = prefix + ".pvdhandler";
  std::ofstream output(filename.c_str());
  output << times_and_names_.size() << std::endl;
  output << "Time File" << std::endl;
  for (unsigned int i = 0; i < times_and_names_.size(); ++i)
    {
      output << times_and_names_[i].first << " " << times_and_names_[i].second
             << std::endl;
    }
}

void
PVDHandler::read(std::string prefix)
{
  times_and_names_.clear();
  std::string   filename = prefix + ".pvdhandler";
  std::ifstream input(filename.c_str());
  if (!input)
    {
      throw("Unable to open file");
    }
  std::string  buffer;
  unsigned int size;
  input >> size;

  std::getline(input, buffer);
  std::getline(input, buffer);
  for (unsigned int i = 0; i < size; ++i)
    {
      double      time;
      std::string filename;
      input >> time;
      input >> filename;
      append(time, filename);
    }

  if (size != times_and_names_.size())
    throw std::runtime_error("Error when reading pvd restart file ");
}
