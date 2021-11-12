
#include <core/utilities.h>


template <int dim, typename T>
TableHandler
make_table_scalars_tensors(
  const std::vector<T> &             independent_vector,
  const std::string &                independent_column_name,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string> &   dependent_column_name,
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


template <int dim>
TableHandler
make_table_tensors_tensors(
  const std::vector<Tensor<1, dim>> &independent_vector,
  const std::vector<std::string> &   independent_column_name,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string> &   dependent_column_name,
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
  const std::vector<std::string> &   independent_column_name,
  const std::vector<double> &        dependent_vector,
  const std::string &                dependent_column_name,
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
fill_table_from_file(TableHandler &    table,
                     std::string       file,
                     const std::string delimiter)
{
  table.clear();
  std::ifstream myfile(file);
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

template TableHandler
make_table_scalars_tensors(
  const std::vector<double> &      independent_values,
  const std::string &              independent_column_name,
  const std::vector<Tensor<1, 2>> &dependent_vector,
  const std::vector<std::string> & dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<double> &      independent_values,
  const std::string &              independent_column_name,
  const std::vector<Tensor<1, 3>> &dependent_vector,
  const std::vector<std::string> & dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<unsigned int> &independent_values,
  const std::string &              independent_column_name,
  const std::vector<Tensor<1, 2>> &dependent_vector,
  const std::vector<std::string> & dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<unsigned int> &independent_values,
  const std::string &              independent_column_name,
  const std::vector<Tensor<1, 3>> &dependent_vector,
  const std::vector<std::string> & dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<int> &         independent_values,
  const std::string &              independent_column_name,
  const std::vector<Tensor<1, 2>> &dependent_vector,
  const std::vector<std::string> & dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_scalars_tensors(
  const std::vector<int> &         independent_values,
  const std::string &              independent_column_name,
  const std::vector<Tensor<1, 3>> &dependent_vector,
  const std::vector<std::string> & dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_tensors_tensors(
  const std::vector<Tensor<1, 2>> &independent_values,
  const std::vector<std::string> & independent_column_name,
  const std::vector<Tensor<1, 2>> &dependent_vector,
  const std::vector<std::string> & dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_tensors_tensors(
  const std::vector<Tensor<1, 3>> &independent_values,
  const std::vector<std::string> & independent_column_name,
  const std::vector<Tensor<1, 3>> &dependent_vector,
  const std::vector<std::string> & dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_tensors_scalars(
  const std::vector<Tensor<1, 2>> &independent_vector,
  const std::vector<std::string> & independent_column_name,
  const std::vector<double> &      dependent_values,
  const std::string &              dependent_column_name,
  const unsigned int               display_precision);

template TableHandler
make_table_tensors_scalars(
  const std::vector<Tensor<1, 3>> &independent_vector,
  const std::vector<std::string> & independent_column_name,
  const std::vector<double> &      dependent_values,
  const std::string &              dependent_column_name,
  const unsigned int               display_precision);
