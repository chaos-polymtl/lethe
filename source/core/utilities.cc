
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


template <int dim>
typename DoFHandler<dim>::active_cell_iterator
find_cell_around_point_with_tree(const DoFHandler<dim> &dof_handler,
                                 Point<dim>             point)
{
  // find cell around point using tree properties instead of looping on all cell

  // define stuff for the search
  MappingQ1<dim> map;
  const auto &   cell_iterator = dof_handler.cell_iterators_on_level(0);
  unsigned int   max_childs    = GeometryInfo<dim>::max_children_per_cell;
  typename DoFHandler<dim>::cell_iterator best_cell_iter;

  bool cell_0_found = false;

  // loop on the cell on lvl 0 of the mesh
  for (const auto &cell : cell_iterator)
    {
      try
        {
          const Point<dim, double> p_cell =
            map.transform_real_to_unit_cell(cell, point);

          const double dist = GeometryInfo<dim>::distance_to_unit_cell(p_cell);

          if (dist == 0)
            {
              // cell on lvl 0 found
              cell_0_found   = true;
              best_cell_iter = cell;
            }
        }
      catch (const typename MappingQGeneric<dim>::ExcTransformationFailed &)
        {}
    }

  if (cell_0_found)
    {
      // the cell on lvl 0 contain the point so now loop on the childs of
      // this cell when we found the child of the cell that containt it we
      // stop and loop if the cell is active if the cell is not active loop
      // on the child of the cell. Repeat
      unsigned int lvl = 0;
      while (best_cell_iter->is_active() == false)
        {
          bool         cell_found = false;
          double       best_dist  = DBL_MAX;
          unsigned int best_index = 0;
          for (unsigned int i = 0; i < max_childs; ++i)
            {
              try
                {
                  const Point<dim, double> p_cell =
                    map.transform_real_to_unit_cell(best_cell_iter->child(i),
                                                    point);
                  const double dist =
                    GeometryInfo<dim>::distance_to_unit_cell(p_cell);
                  bool inside = true;


                  if (dist <= best_dist and inside)
                    {
                      best_dist  = dist;
                      best_index = i;
                      cell_found = true;
                      if (dist == 0)
                        break;
                    }
                }
              catch (
                const typename MappingQGeneric<dim>::ExcTransformationFailed &)
                {}
            }

          best_cell_iter = best_cell_iter->child(best_index);

          if (cell_found == false)
            {
              std::cout << "cell not found" << point << std::endl;
              break;
            }

          lvl += 1;
        }

      //      break;
    }


  return best_cell_iter;
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

template typename DoFHandler<3>::active_cell_iterator
find_cell_around_point_with_tree(const DoFHandler<3> &dof_handler,
                                 Point<3>             point);
template typename DoFHandler<2>::active_cell_iterator
find_cell_around_point_with_tree(const DoFHandler<2> &dof_handler,
                                 Point<2>             point);
