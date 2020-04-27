
#include <core/utilities.h>


template <int dim>
TableHandler
  make_table_from_vector_of_tensors(std::vector<Tensor<1, dim>> vector,
                                    std::vector<std::string>    column_names,
                                    const unsigned int display_precision)
{
  TableHandler table;
  for (unsigned int i = 0; i < vector.size(); ++i)
    {
      table.add_value(column_names[0], i);

      table.add_value(column_names[1], vector[i][0]);
      table.set_precision(column_names[1], display_precision);

      table.add_value(column_names[2], vector[i][1]);
      table.set_precision(column_names[2], display_precision);
      if (dim == 3)
        {
          table.add_value(column_names[3], vector[i][2]);
          table.set_precision(column_names[3], display_precision);
        }
    }
  return table;
}



template TableHandler
  make_table_from_vector_of_tensors(std::vector<Tensor<1, 2>> vector,
                                    std::vector<std::string>  column_names,
                                    const unsigned int display_precision);

template TableHandler
  make_table_from_vector_of_tensors(std::vector<Tensor<1, 3>> vector,
                                    std::vector<std::string>  column_names,
                                    const unsigned int display_precision);
