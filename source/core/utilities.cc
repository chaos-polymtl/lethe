
#include <core/utilities.h>

#if __GNUC__ > 7
#  include <filesystem>
#endif

template <int dim, typename T>
TableHandler
make_table_scalars_tensors(
  const std::vector<T>              &independent_vector,
  const std::string                 &independent_column_name,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string>    &dependent_column_name,
  const unsigned int                 display_precision)
{
  AssertDimension(independent_vector.size(), dependent_vector.size());
  AssertDimension(dependent_column_name.size(), dim);


  TableHandler table;

  for (unsigned int i = 0; i < dependent_vector.size(); ++i)
    {
      table.add_value(independent_column_name, independent_vector[i]);
      for (unsigned int d = 0; d < dim; ++d)
        table.add_value(dependent_column_name[d], dependent_vector[i][d]);
    }

  table.set_precision(independent_column_name, display_precision);
  for (unsigned int d = 0; d < dim; ++d)
    table.set_precision(dependent_column_name[d], display_precision);

  return table;
}

template <int dim, typename T>
TableHandler
make_table_scalars_tensors(
  const std::vector<T>                           &independent_vector,
  const std::string                              &independent_column_name,
  const std::vector<std::vector<Tensor<1, dim>>> &dependent_vectors,
  const std::vector<std::string>                 &dependent_column_name,
  const unsigned int                              display_precision)
{
  AssertDimension(dependent_column_name.size(), dependent_vectors.size() * dim);

  TableHandler table;
  unsigned int vect_index = 0;

  for (unsigned int i = 0; i < dependent_vectors[0].size(); ++i)
    table.add_value(independent_column_name, independent_vector[i]);

  for (auto &vect : dependent_vectors)
    {
      AssertDimension(independent_vector.size(), vect.size());
      for (unsigned int i = 0; i < vect.size(); ++i)
        {
          for (unsigned int d = 0; d < dim; ++d)
            {
              table.add_value(dependent_column_name[d + vect_index],
                              vect[i][d]);
              table.set_precision(dependent_column_name[d + vect_index],
                                  display_precision);
            }
        }
      vect_index += dim;
    }

  table.set_precision(independent_column_name, display_precision);

  return table;
}


template <int dim>
TableHandler
make_table_tensors_tensors(
  const std::vector<Tensor<1, dim>> &independent_vector,
  const std::vector<std::string>    &independent_column_name,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string>    &dependent_column_name,
  const unsigned int                 display_precision)
{
  AssertDimension(independent_vector.size(), dependent_vector.size());
  AssertDimension(independent_column_name.size(), dim);
  AssertDimension(dependent_column_name.size(), dim);

  TableHandler table;

  for (unsigned int i = 0; i < dependent_vector.size(); ++i)
    {
      for (unsigned int d = 0; d < dim; ++d)
        {
          table.add_value(independent_column_name[d], independent_vector[i][d]);
          table.add_value(dependent_column_name[d], dependent_vector[i][d]);
        }
    }

  for (unsigned int d = 0; d < dim; ++d)
    {
      table.set_precision(independent_column_name[d], display_precision);
      table.set_precision(dependent_column_name[d], display_precision);
    }

  return table;
}

template <int dim>
TableHandler
make_table_tensors_scalars(
  const std::vector<Tensor<1, dim>> &independent_vector,
  const std::vector<std::string>    &independent_column_name,
  const std::vector<double>         &dependent_vector,
  const std::string                 &dependent_column_name,
  const unsigned int                 display_precision)
{
  AssertDimension(independent_vector.size(), dependent_vector.size());
  AssertDimension(independent_column_name.size(), dim);


  TableHandler table;

  for (unsigned int i = 0; i < dependent_vector.size(); ++i)
    {
      table.add_value(dependent_column_name, dependent_vector[i]);
      for (unsigned int d = 0; d < dim; ++d)
        table.add_value(independent_column_name[d], independent_vector[i][d]);
    }

  table.set_precision(dependent_column_name, display_precision);
  for (unsigned int d = 0; d < dim; ++d)
    table.set_precision(independent_column_name[d], display_precision);

  return table;
}

void
fill_table_from_file(TableHandler     &table,
                     const std::string file_name,
                     const std::string delimiter)
{
  table.clear();
  std::ifstream myfile(file_name);
  // open the file
  if (myfile.is_open())
    {
      std::string              line;
      std::vector<std::string> vector_of_column_names;
      std::vector<double>      line_of_data;
      unsigned int             line_count = 0;

      while (std::getline(myfile, line))
        {
          // read the line and clean the resulting vector
          std::vector<std::string> list_of_words_base =
            Utilities::split_string_list(line, delimiter);
          std::vector<std::string> list_of_words_clean;
          for (unsigned int i = 0; i < list_of_words_base.size(); ++i)
            {
              if (list_of_words_base[i] != "")
                {
                  list_of_words_clean.push_back(list_of_words_base[i]);
                }
            }
          //  If it's the first line, we only initialize the variable names.
          if (line_count != 0)
            {
              line_of_data = Utilities::string_to_double(list_of_words_clean);
              for (unsigned int i = 0; i < line_of_data.size(); ++i)
                {
                  table.add_value(vector_of_column_names[i], line_of_data[i]);
                }
            }
          else
            {
              // the line contains words we assume these are the column
              vector_of_column_names = list_of_words_clean;
            }
          ++line_count;
        }
      myfile.close();
    }
  else
    std::cout << "Unable to open file";
}

void
fill_vectors_from_file(std::map<std::string, std::vector<double>> &map,
                       std::string                                 file,
                       const std::string                           delimiter)
{
  // fill a pair, first being a vector of vector name and the second being the
  // vector of vector associated with the vector name.


  std::ifstream myfile(file);
  // open the file.
  if (myfile.is_open())
    {
      std::string              line;
      std::vector<std::string> column_names;
      std::vector<double>      line_of_data;
      unsigned int             line_count = 0;

      while (std::getline(myfile, line))
        {
          // read the line and clean the resulting vector.
          std::vector<std::string> list_of_words_base =
            Utilities::split_string_list(line, delimiter);
          std::vector<std::string> list_of_words_clean;
          for (unsigned int i = 0; i < list_of_words_base.size(); ++i)
            {
              if (list_of_words_base[i] != "")
                {
                  list_of_words_clean.push_back(list_of_words_base[i]);
                }
            }
          // check if the line is contained words or numbers.
          if (line_count != 0)
            {
              line_of_data = Utilities::string_to_double(list_of_words_clean);
              for (unsigned int i = 0; i < line_of_data.size(); ++i)
                {
                  map[column_names[i]].push_back(line_of_data[i]);
                }
            }
          else
            {
              // the line contains words, we assume these are the columns names.
              column_names = list_of_words_clean;
              for (unsigned int i = 0; i < list_of_words_clean.size(); ++i)
                {
                  std::vector<double> base_vector;
                  map[column_names[i]] = base_vector;
                }
            }
          ++line_count;
        }
      myfile.close();
    }
  else
    std::cout << "Unable to open file";
}

void
fill_string_vectors_from_file(std::map<std::string, std::vector<std::string>> &map,
                       std::string                                 file,
                       const std::string                           delimiter)
{
  // fill a pair, first being a vector of vector name and the second being the
  // vector of vector associated with the vector name.


  std::ifstream myfile(file);
  // open the file.
  if (myfile.is_open())
    {
      std::string              line;
      std::vector<std::string> column_names;
      std::vector<std::string>      line_of_data;
      unsigned int             line_count = 0;

      while (std::getline(myfile, line))
        {
          // remove spaces
          //line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
          // read the line and clean the resulting vector.
          std::vector<std::string> list_of_words_base =
            Utilities::split_string_list(line, delimiter);
          std::vector<std::string> list_of_words_clean;
          for (unsigned int i = 0; i < list_of_words_base.size(); ++i)
            {
              if (list_of_words_base[i] != "")
                {
                  list_of_words_clean.push_back(list_of_words_base[i]);
                }
            }
          // check if the line is contained words or numbers.
          if (line_count != 0)
            {
              line_of_data = list_of_words_clean;
              for (unsigned int i = 0; i < line_of_data.size(); ++i)
                {
                  map[column_names[i]].push_back(line_of_data[i]);
                }
            }
          else
            {
              // the line contains words, we assume these are the columns names.
              column_names = list_of_words_clean;
              for (unsigned int i = 0; i < list_of_words_clean.size(); ++i)
                {
                  std::vector<std::string> base_vector;
                  map[column_names[i]] = base_vector;
                }
            }
          ++line_count;
        }
      myfile.close();
    }
  else
    std::cout << "Unable to open file";
}

void
create_output_folder(const std::string &dirname)
{
#if __GNUC__ > 7
  std::filesystem::create_directory(dirname);
#else
  // mkdir(dirname.c_str(), 0755);
#endif
}

template TableHandler
make_table_scalars_tensors(
  const std::vector<double>       &independent_values,
  const std::string               &independent_column_name,
  const std::vector<Tensor<1, 2>> &dependent_vector,
  const std::vector<std::string>  &dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<double>       &independent_values,
  const std::string               &independent_column_name,
  const std::vector<Tensor<1, 3>> &dependent_vector,
  const std::vector<std::string>  &dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<double>                    &independent_values,
  const std::string                            &independent_column_name,
  const std::vector<std::vector<Tensor<1, 2>>> &dependent_vector,
  const std::vector<std::string>               &dependent_column_name,
  const unsigned int                            display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<double>                    &independent_values,
  const std::string                            &independent_column_name,
  const std::vector<std::vector<Tensor<1, 3>>> &dependent_vector,
  const std::vector<std::string>               &dependent_column_name,
  const unsigned int                            display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<unsigned int> &independent_values,
  const std::string               &independent_column_name,
  const std::vector<Tensor<1, 2>> &dependent_vector,
  const std::vector<std::string>  &dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<unsigned int> &independent_values,
  const std::string               &independent_column_name,
  const std::vector<Tensor<1, 3>> &dependent_vector,
  const std::vector<std::string>  &dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<unsigned int>              &independent_values,
  const std::string                            &independent_column_name,
  const std::vector<std::vector<Tensor<1, 2>>> &dependent_vector,
  const std::vector<std::string>               &dependent_column_name,
  const unsigned int                            display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<unsigned int>              &independent_values,
  const std::string                            &independent_column_name,
  const std::vector<std::vector<Tensor<1, 3>>> &dependent_vector,
  const std::vector<std::string>               &dependent_column_name,
  const unsigned int                            display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<int>          &independent_values,
  const std::string               &independent_column_name,
  const std::vector<Tensor<1, 2>> &dependent_vector,
  const std::vector<std::string>  &dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<int>          &independent_values,
  const std::string               &independent_column_name,
  const std::vector<Tensor<1, 3>> &dependent_vector,
  const std::vector<std::string>  &dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<int>                       &independent_values,
  const std::string                            &independent_column_name,
  const std::vector<std::vector<Tensor<1, 2>>> &dependent_vector,
  const std::vector<std::string>               &dependent_column_name,
  const unsigned int                            display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<int>                       &independent_values,
  const std::string                            &independent_column_name,
  const std::vector<std::vector<Tensor<1, 3>>> &dependent_vector,
  const std::vector<std::string>               &dependent_column_name,
  const unsigned int                            display_precision);

template TableHandler
make_table_tensors_tensors(
  const std::vector<Tensor<1, 2>> &independent_values,
  const std::vector<std::string>  &independent_column_name,
  const std::vector<Tensor<1, 2>> &dependent_vector,
  const std::vector<std::string>  &dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_tensors_tensors(
  const std::vector<Tensor<1, 3>> &independent_values,
  const std::vector<std::string>  &independent_column_name,
  const std::vector<Tensor<1, 3>> &dependent_vector,
  const std::vector<std::string>  &dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_tensors_scalars(
  const std::vector<Tensor<1, 2>> &independent_vector,
  const std::vector<std::string>  &independent_column_name,
  const std::vector<double>       &dependent_values,
  const std::string               &dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_tensors_scalars(
  const std::vector<Tensor<1, 3>> &independent_vector,
  const std::vector<std::string>  &independent_column_name,
  const std::vector<double>       &dependent_values,
  const std::string               &dependent_column_name,
  const unsigned int               display_precision);


std::string
get_last_value_of_parameter(const std::string &file_name,
                            const std::string &parameter_name)
{
  std::string return_value;

  std::ifstream x_file(file_name);
  AssertThrow(x_file.fail() == false, ExcIO());

  while (x_file)
    {
      // Get one line and then match a regex to it that matches the parameter
      // we are looking for. Before we do that, strip spaces from the front
      // and back of the line:
      std::string line;
      std::getline(x_file, line);

      while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
        line.erase(0, 1);
      while ((line.size() > 0) &&
             (line[line.size() - 1] == ' ' || line[line.size() - 1] == '\t'))
        line.erase(line.size() - 1, std::string::npos);

      std::match_results<std::string::const_iterator> matches;
      const std::string                               regex =
        "set[ \t]+" + parameter_name + "[ \t]*=[ \t]*(.*)";
      if (std::regex_match(line, matches, std::regex(regex)))
        {
          // Since the line as a whole matched, the 'matches' variable needs to
          // contain two entries: [0] denotes the whole string, and [1] the
          // one that was matched by the '(.*)' expression.
          Assert(matches.size() == 2, dealii::ExcInternalError());
          return_value = std::string(matches[1].first, matches[1].second);
        }
    }

  return return_value;
}


int
get_max_value_of_parameter(const std::string &file_name,
                           const std::string &parameter_name)
{
  std::string return_string;
  int         return_value = -100000;

  std::ifstream x_file(file_name);
  AssertThrow(x_file.fail() == false, ExcIO());

  while (x_file)
    {
      // Get one line and then match a regex to it that matches the parameter
      // we are looking for. Before we do that, strip spaces from the front
      // and back of the line:
      std::string line;
      std::getline(x_file, line);

      while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
        line.erase(0, 1);
      while ((line.size() > 0) &&
             (line[line.size() - 1] == ' ' || line[line.size() - 1] == '\t'))
        line.erase(line.size() - 1, std::string::npos);

      std::match_results<std::string::const_iterator> matches;
      const std::string                               regex =
        "set[ \t]+" + parameter_name + "[ \t]*=[ \t]*(.*)";
      if (std::regex_match(line, matches, std::regex(regex)))
        {
          // Since the line as a whole matched, the 'matches' variable needs to
          // contain two entries: [0] denotes the whole string, and [1] the
          // one that was matched by the '(.*)' expression.
          Assert(matches.size() == 2, dealii::ExcInternalError());
          return_string = std::string(matches[1].first, matches[1].second);
          return_value =
            std::max(return_value, Utilities::string_to_int(return_string));
        }
    }

  return return_value;
}



unsigned int
get_dimension(const std::string &file_name)
{
  const std::string dimension =
    get_last_value_of_parameter(file_name, "dimension");

  if (dimension.size() > 0)
    {
      // Extracted from ASPECT
      // A common problem is that people have .prm files that were generated
      // on Windows, but then run this on Linux where the line endings are
      // different. This is pernicious because it means that the conversion
      // of a string such as "2\r" to an integer fails, but if we print
      // this string, it comes out completely garbled because it contains
      // a carriage-return without a newline -- so the error message looks
      // like this:
      //
      //    >.  While reading the dimension from the input file, ASPECT found a
      //    string that can not be converted to an integer: <2
      //
      // Note how the end of the error message overwrites the beginning
      // of the line.
      //
      // To avoid this kind of error, specifically test up front that the
      // text in question does not contain '\r' characters. If we are on
      // linux, then this kind of character would means that the line endings
      // are wrong. On the other hand, if we are on windows, then the
      // getline command we have used in finding 'dimension' would have
      // filtered it out. So its presence points to a problem.

      AssertThrow(
        dimension.find('\r') == std::string::npos,
        dealii::ExcMessage(
          "It appears that your input file uses Windows-style "
          "line endings ('\\r\\n') but you are running on a system where "
          "the C++ run time environment expects input files to have "
          "Unix-style line endings ('\\n'). You need to convert your "
          "input file to use the correct line endings before running "
          "ASPECT with it."));
      try
        {
          return dealii::Utilities::string_to_int(dimension);
        }
      catch (...)
        {
          AssertThrow(false,
                      dealii::ExcMessage(
                        "While reading the dimension from the input file, "
                        "Lethe found a string that can not be converted to "
                        "an integer: <" +
                        dimension + ">."));
          return 0; // we should never get here.
        }
    }
  else
    {
      AssertThrow(
        false,
        dealii::ExcMessage(
          "While reading the dimension from the input file, "
          "Lethe found a value that is neither 2 or 3. Since August 2023, "
          "Lethe requires that the user explicitly specify the dimension of the problem within the parameter file. This can be achieved by adding set dimension = 2 or set dimension = 3 within the parameter file"));
    }
}


int
get_max_number_of_boundary_conditions(const std::string &file_name)
{
  int max_number_of_boundary_conditions =
    get_max_value_of_parameter(file_name, "number");

  AssertThrow(
    max_number_of_boundary_conditions >= 0,
    dealii::ExcMessage(
      "Your parameter file does not contain any indication for the number of boundary conditions for any physics supported by Lethe. Since November 2023, Lethe requires that a \"boundary conditions\" subsection is present with at least \"number=0\" "));

  return std::max(max_number_of_boundary_conditions, 0);
}
