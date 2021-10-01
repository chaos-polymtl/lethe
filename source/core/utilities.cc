
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
fill_table_from_file(TableHandler & table, std::string file,const std::string delimiter ){
  std::string   line;
  table.clear();
  std::ifstream myfile(file);
  if (myfile.is_open())
    {
      std::vector<std::string> vector_of_column_names;
      std::vector<double> line_of_data;

      while (std::getline(myfile, line))
        {
          std::vector<std::string> list_of_words_base=
            Utilities::split_string_list(line, delimiter);
          std::vector<std::string> list_of_words_clean;
          for(unsigned int i=0 ; i<list_of_words_base.size();++i)
            {
              if(list_of_words_base[i]!=""){
                  list_of_words_clean.push_back(list_of_words_base[i]);
                }
            }
          //check if the line is contained words or numbers
          std::cout<< "list of word clean"<<list_of_words_clean.size()<<std::endl;
          try{
              line_of_data=Utilities::string_to_double(list_of_words_clean);
              for(unsigned int i=0 ; i<line_of_data.size();++i){
                  table.add_value( vector_of_column_names[i], line_of_data[i]);

                }
            }
          catch(...){
              //the line contains words we assume these are the column
              vector_of_column_names=list_of_words_clean;
            }


        }
      myfile.close();
    }
  else
    std::cout << "Unable to open file";

}

void
fill_vectors_from_file(std::pair<std::vector<std::string>,std::vector<std::vector<double>>>& vectors, std::string file,const std::string delimiter ){
  //fill a pair, first being a vector of vector name and the second being the vector of vector associated with the vector name.

  std::string   line;

  std::ifstream myfile(file);
  std::cout <<"file is open"<<std::endl;
  if (myfile.is_open())
    {
      std::vector<std::string> vector_of_column_names;
      std::vector<double> line_of_data;

      while (std::getline(myfile, line))
        {

          std::vector<std::string> list_of_words_base=
            Utilities::split_string_list(line, delimiter);
          std::vector<std::string> list_of_words_clean;
          for(unsigned int i=0 ; i<list_of_words_base.size();++i)
            {
              if(list_of_words_base[i]!=""){
                  list_of_words_clean.push_back(list_of_words_base[i]);
                }
            }
          //check if the line a word or number
          try{
              line_of_data=Utilities::string_to_double(list_of_words_clean);
              for(unsigned int i=0 ; i<line_of_data.size();++i){
                  vectors.second[i].push_back(line_of_data[i]);
                }
            }
          catch(...){
              //the line contains words, we assume these are the columns names
              vectors.first=list_of_words_clean;
              vectors.second.resize(vectors.first.size());

            }
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
