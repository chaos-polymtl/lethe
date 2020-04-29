
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
