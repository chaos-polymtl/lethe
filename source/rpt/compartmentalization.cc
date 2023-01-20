#include <core/lethe_grid_tools.h>
#include <core/solutions_output.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/block_vector.h>

#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/data_out.h>

#include <bits/stdc++.h>
#include <rpt/compartmentalization.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

template <int dim>
Compartmentalization<dim>::Compartmentalization(
  CPCalculatingParameters &CPparameters)
  : computing_timer(std::cout, TimerOutput::summary, TimerOutput::wall_times)
  , mpi_communicator(MPI_COMM_WORLD)
  , triangulation(mpi_communicator)
  , cp_parameters(CPparameters)
{
  std::cout << cp_parameters.cp_param.subdivisions << std::endl;
  std::cout << cp_parameters.cp_param.CFD_input_velocity << std::endl;
}

template <int dim>
void
Compartmentalization<dim>::generate_cylindrical_grid()
{
  GridOut grid_out;
  GridGenerator::subdivided_cylinder(triangulation, 1, 0.03, 0.06);
  triangulation.refine_global(1);
  std::ofstream output("triangulation.vtk");
  grid_out.write_vtk(triangulation, output);
  std::vector<double>     cell_index_output;
  std::vector<Point<dim>> points;
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // std::vector<double> center;
      Point<3, double> center;
      if (cell->is_locally_owned())
        {
          center = cell->center();
          cells_and_indices.insert({cell->index(), cell});
          points.push_back(center);
          center.clear();
        }
    }

  // output the center of cell to feed the Pyvista and extract the physical
  // cfd or other simulation
  std::string   filename = "cell_centers";
  std::ofstream myfile;
  std::string   sep = " ";
  myfile.open(filename);
  for (auto &point : points)
    {
      for (int j = 0; j < 3; j++)
        {
          myfile << point[j] << sep;
        }
      myfile << "\n";
    }
  myfile.close();
}

template <int dim>
auto
Compartmentalization<dim>::read_electric_field()
  -> std::vector<std::vector<double>>
{
  // Read the electric field and write it in a vector
  std::vector<double> electric_field;
  std::ifstream       fin("input_values_emw");
  double              element;
  while (fin >> element)
    {
      electric_field.push_back(element);
    }

  // Fill the container with the electromagnetic field
  for (const auto &[cell_index, cell] : cells_and_indices)
    {
      std::vector<double> fill_matrix;
      fill_matrix.push_back(cell_index);
      fill_matrix.push_back(electric_field[cell_index]);
      primer_matrix_index_emw.push_back(fill_matrix);
      fill_matrix.clear();
    }
  return primer_matrix_index_emw;
}

template <int dim>
auto
Compartmentalization<dim>::read_velocity_magnitude()
  -> std::map<typename Triangulation<dim>::active_cell_iterator, double>
{
  // Read the velocity magnitude and write it in a vector
  std::vector<double> velocity_magnitude;
  std::ifstream       fin("input_values_velocity_magnitude");
  double              element;
  while (fin >> element)
    {
      velocity_magnitude.push_back(element);
    }

  // this for loop fill the container with the properties
  //  (velocity or electromagnetic field)
  for (const auto &[cell_index, cell] : cells_and_indices)
    {
      primer_matrix_index_velocity.insert(
        {cell, velocity_magnitude[cell_index]});
    }
  return primer_matrix_index_velocity;
}

template <int dim>
auto
Compartmentalization<dim>::sort_agglomeration_deagglomeration_emw()
  -> std::map<int,
              std::vector<typename Triangulation<dim>::active_cell_iterator>>
{
  // sort_agglomeration_deagglomeration_emw, first sort cells based on theis
  // value,
  //  then agglomerate cells in a certain threshold and then deagglomerate
  //  clusters to consider neighbour only neighbour cells in a group
  std::map<double, typename Triangulation<dim>::active_cell_iterator>
    primer_map_cell_index;

  // vector_cell is a vector, later it will store the sorted cell
  //  (cell itself, not the cell_index or value)
  std::vector<typename Triangulation<dim>::active_cell_iterator> vector_cell;

  for (const auto &[cell_index, cell] : cells_and_indices)
    {
      primer_map_cell_index.insert({cell->index(), cell});
    }

  // sorting, will sort cells based on their value and the
  // index associated with the cell will move along with value
  std::sort(primer_matrix_index_emw.begin(),
            primer_matrix_index_emw.end(),
            [](const std::vector<double> &a, const std::vector<double> &b) {
              return a[1] > b[1];
            });

  // now, all the cells are sorted based on their values
  //  (a[0]=cell_index and a[1]=value of physical property)
  // this for loop will find the cell associated with the sorted cell_index and
  // will rearrange it based on the sorted value,
  //  this means after sorting the matrix
  // we need to update and sort the map which
  //  contains only the cell info and not the value
  for (auto &i : primer_matrix_index_emw)
    {
      double     index  = i[0];
      const auto search = primer_map_cell_index.find(index);

      // vector_cell stores only the sorted cells
      vector_cell.push_back(search->second);
    }

  // Since we have at least one compartment the iterator starts from 1
  int number_of_compartments = 1;

  // this "k" iterator will avoid left over the last cell in case
  //  that it must be clustered in one separate group
  unsigned int k = 1;

  // first_set, is a map containing each compartment and their cell
  // before the deagglomeration step
  // since it is possible that some cells locate in same compartment but,
  // they are not connected
  // directly or indirectly and the need deagglomeration
  std::map<int, std::vector<typename Triangulation<dim>::active_cell_iterator>>
    first_set;

  // vector will temporarily contain the cell of each compartment
  std::vector<typename Triangulation<dim>::active_cell_iterator> vector;

  // tol is the threshold of the clustering for difference between the
  // physical property values
  double Tol = 1;
  vector.push_back(vector_cell[0]);

  // max is the highest physical property
  double max = primer_matrix_index_emw[0][1];
  for (unsigned int i = 0; i < primer_matrix_index_emw.size() - 1; i++)
    {
      while (max - primer_matrix_index_emw[i + 1][1] < Tol)
        {
          vector.push_back(vector_cell[i + 1]);
          k = k + 1;
          if ((i + 1) < primer_matrix_index_emw.size() - 1)
            {
              i = i + 1;
            }
          else
            {
              break;
            }
        }
      first_set.insert({number_of_compartments, vector});
      number_of_compartments++;
      vector.clear();
      vector.push_back(vector_cell[i + 1]);
      k   = k + 1;
      max = primer_matrix_index_emw[i + 1][1];
    }
  if (primer_matrix_index_emw.size() == k)
    {
      number_of_compartments++;
      vector.clear();
      vector.push_back(vector_cell[primer_matrix_index_emw.size() - 1]);
      first_set.insert({number_of_compartments, vector});
    }
  double         compartment_id = 1;
  Vector<double> compartments;
  compartments.reinit(triangulation.n_active_cells());
  for (const auto &pair : first_set)
    {
      for (typename Triangulation<dim>::active_cell_iterator d : pair.second)
        {
          compartments[d->active_cell_index()] = compartment_id;
        }
      compartment_id = compartment_id + 1;
    }

  //  Start the deagglomeration
  int key_final_set = 1;

  // This is a vector that contains those cells that are neighbor in first
  //  cluster, will make first cluster out of the first one
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    filling_vector_declustered_set;

  // In this set I put all the cells from the first cluster to deagglomerate
  std::set<typename Triangulation<dim>::active_cell_iterator> set_to_search;

  // In this set I put the neighbors of each cell that I loop over
  std::set<typename Triangulation<dim>::active_cell_iterator> temporary_set;

  for (const auto &pair : first_set)
    {
      // Fill the set with the cells into the first cluster
      for (const typename Triangulation<dim>::active_cell_iterator &d :
           pair.second)
        {
          set_to_search.insert(d);
        }
      while (set_to_search.size() != 0)
        {
          // I consider the first cell into set_to_search as starting point
          //  since it is part of unique cluster
          filling_vector_declustered_set.push_back(*set_to_search.begin());

          // now we should delete it from the set_to_search since it got its
          //  forever home
          set_to_search.erase(*set_to_search.begin());

          // at first, search for the neighbour of this first cell that we
          //  have added manually
          for (unsigned int i = 0; i < 6; i++)
            {
              typename Triangulation<dim>::active_cell_iterator
                neighbour_first_element =
                  filling_vector_declustered_set[0]->neighbor(i);
              if (set_to_search.count(neighbour_first_element) == 1)
                {
                  // add those neighbours which they are in one cluster to the
                  // vector
                  filling_vector_declustered_set.push_back(
                    neighbour_first_element);

                  // add them to temporary set to verify their neighbors
                  temporary_set.insert(neighbour_first_element);
                  set_to_search.erase(neighbour_first_element);
                }
            }
          while (temporary_set.size() != 0)
            {
              for (unsigned int j = 0; j < temporary_set.size(); j++)
                {
                  for (unsigned int i = 0; i < 6; i++)
                    {
                      typename Triangulation<dim>::active_cell_iterator
                           temporary_set_element = *temporary_set.begin();
                      auto neighbour = temporary_set_element->neighbor(i);
                      if (set_to_search.count(neighbour) == 1)
                        {
                          temporary_set.insert(neighbour);
                          filling_vector_declustered_set.push_back(neighbour);
                          set_to_search.erase(neighbour);
                        }
                    }
                  temporary_set.erase(temporary_set.begin());
                }
            }
          final_set.insert({key_final_set, filling_vector_declustered_set});
          filling_vector_declustered_set.clear();
          key_final_set = key_final_set + 1;
        }
    }
  double         compartment_final_id = 1;
  Vector<double> compartments_emw;
  compartments_emw.reinit(triangulation.n_active_cells());
  for (const auto &pair : final_set)
    {
      for (typename Triangulation<dim>::active_cell_iterator d : pair.second)
        {
          compartments_emw[d->active_cell_index()] = compartment_final_id;
        }
      compartment_final_id = compartment_final_id + 1;
    }
  write_file_compartments_first_field(compartments_emw, 1.0, 1);
  return final_set;
}

template <int dim>
void
Compartmentalization<dim>::overlaid_map()
{
  // loop over compartments
  // loop over cells
  // sort their velocity value
  // agglomerate based on velocity
  // deagglomerate if they are not connected
  // add to the overlaid_map
  std::map<double, typename Triangulation<dim>::active_cell_iterator>
    primer_map_cell_index;

  std::vector<typename Triangulation<dim>::active_cell_iterator> vector_cell;
  std::vector<std::vector<double>>                               inner_matrix;

  // loop over compartments
  int key_final_set = 1;
  for (const auto &pair : final_set)
    {
      for (typename Triangulation<dim>::active_cell_iterator d : pair.second)
        {
          primer_map_cell_index.insert({d->index(), d});
          std::vector<double> fill_matrix;
          fill_matrix.push_back(d->index());
          fill_matrix.push_back(primer_matrix_index_velocity.find(d)->second);
          inner_matrix.push_back(fill_matrix);
          fill_matrix.clear();
        }
      std::sort(inner_matrix.begin(),
                inner_matrix.end(),
                [](const std::vector<double> &a, const std::vector<double> &b) {
                  return a[1] > b[1];
                });
      for (auto &i : inner_matrix)
        {
          double index = i[0];

          // const auto search = primer_map_cell_index.find(index);
          // gere vector_cell stores only the sorted cells
          vector_cell.push_back(primer_map_cell_index.find(index)->second);
        }
      int          number_of_compartments = 1;
      unsigned int k                      = 1;
      std::map<int,
               std::vector<typename Triangulation<dim>::active_cell_iterator>>
        first_set;

      // vector will temporarily contain the cell of each compartment
      std::vector<typename Triangulation<dim>::active_cell_iterator> vector;

      // tol is the threshold of the clustering for difference between the
      // physical property values
      double Tol = 0.7;
      vector.push_back(vector_cell[0]);

      // max is the highest physical property
      double max = inner_matrix[0][1];
      for (unsigned int i = 0; i < inner_matrix.size() - 1; i++)
        {
          while (max - inner_matrix[i + 1][1] < Tol)
            {
              vector.push_back(vector_cell[i + 1]);
              k = k + 1;
              if ((i + 1) < inner_matrix.size() - 1)
                {
                  i = i + 1;
                }
              else
                {
                  break;
                }
            }
          first_set.insert({number_of_compartments, vector});
          number_of_compartments++;
          vector.clear();
          vector.push_back(vector_cell[i + 1]);
          k   = k + 1;
          max = inner_matrix[i + 1][1];
        }
      if (inner_matrix.size() == k)
        {
          number_of_compartments++;
          vector.clear();
          vector.push_back(vector_cell[inner_matrix.size() - 1]);
          first_set.insert({number_of_compartments, vector});
        }

      // This is a vector that contains those cells that are neighbor in first
      //  cluster, will make first cluster out of the first one
      std::vector<typename Triangulation<dim>::active_cell_iterator>
        filling_vector_declustered_set;

      // In this set I put all the cells from the first cluster to deagglomerate
      std::set<typename Triangulation<dim>::active_cell_iterator> set_to_search;

      // In this set I put the neighbors of each cell that I loop over
      std::set<typename Triangulation<dim>::active_cell_iterator> temporary_set;

      for (const auto &pair : first_set)
        {
          // Fill the set with the cells into the first cluster
          for (const typename Triangulation<dim>::active_cell_iterator &d :
               pair.second)
            {
              set_to_search.insert(d);
            }
          while (set_to_search.size() != 0)
            {
              // I consider the first cell into set_to_search as starting point
              //  since it is part of unique cluster
              filling_vector_declustered_set.push_back(*set_to_search.begin());

              // now we should delete it from the set_to_search since it got its
              //  forever home
              set_to_search.erase(*set_to_search.begin());

              // at first, search for the neighbour of this first cell that we
              //  have added manually
              for (unsigned int i = 0; i < 6; i++)
                {
                  typename Triangulation<dim>::active_cell_iterator
                    neighbour_first_element =
                      filling_vector_declustered_set[0]->neighbor(i);
                  if (set_to_search.count(neighbour_first_element) == 1)
                    {
                      // add those neighbours which they are in one cluster to
                      // the vector
                      filling_vector_declustered_set.push_back(
                        neighbour_first_element);

                      // add them to temporary set to verify their neighbors
                      temporary_set.insert(neighbour_first_element);
                      set_to_search.erase(neighbour_first_element);
                    }
                }
              while (temporary_set.size() != 0)
                {
                  for (unsigned int j = 0; j < temporary_set.size(); j++)
                    {
                      for (unsigned int i = 0; i < 6; i++)
                        {
                          typename Triangulation<dim>::active_cell_iterator
                               temporary_set_element = *temporary_set.begin();
                          auto neighbour = temporary_set_element->neighbor(i);
                          if (set_to_search.count(neighbour) == 1)
                            {
                              temporary_set.insert(neighbour);
                              filling_vector_declustered_set.push_back(
                                neighbour);
                              set_to_search.erase(neighbour);
                            }
                        }
                      temporary_set.erase(temporary_set.begin());
                    }
                }
              overlaid_set.insert(
                {key_final_set, filling_vector_declustered_set});
              filling_vector_declustered_set.clear();
              key_final_set = key_final_set + 1;
            }
        }
      for (auto &loop : overlaid_set)
        {
          ultimate_map.insert({loop.first, loop.second});
        }
      first_set.clear();
      vector.clear();
      vector_cell.clear();
      filling_vector_declustered_set.clear();
      set_to_search.clear();
      temporary_set.clear();
      inner_matrix.clear();
      primer_map_cell_index.clear();
    }
  double         compartment_final_id = 1;
  Vector<double> compartments_ultimate;
  compartments_ultimate.reinit(triangulation.n_active_cells());
  for (const auto &pair : ultimate_map)
    {
      for (typename Triangulation<dim>::active_cell_iterator d : pair.second)
        {
          compartments_ultimate[d->active_cell_index()] = compartment_final_id;
        }
      compartment_final_id = compartment_final_id + 1;
    }
  write_ultimate_file(compartments_ultimate, 1.0, 1);
}



template <int dim>
void
Compartmentalization<dim>::write_file_compartments_first_field(
  Vector<double> &    compartments_final,
  const double &      time,
  const unsigned int &step_number)
{
  // This function visualize the clusters
  const std::string folder = "./";

  DataOut<dim> data_out;
  data_out.attach_triangulation(triangulation);
  std::string average_solution_names("Compartments_emw");

  data_out.add_data_vector(compartments_final,
                           average_solution_names,
                           DataOut<dim>::type_cell_data);

  data_out.build_patches();


  write_vtu_and_pvd<dim>(grid_pvdhandler,
                         data_out,
                         folder,
                         "emw",
                         time,
                         step_number,
                         1,
                         mpi_communicator);
}

template <int dim>
void
Compartmentalization<dim>::write_ultimate_file(
  Vector<double> &    compartments_ultimate,
  const double &      time,
  const unsigned int &step_number)
{
  const std::string folder = "./";

  DataOut<dim> data_out;
  data_out.attach_triangulation(triangulation);
  std::string average_solution_names("overlaid_map");

  data_out.add_data_vector(compartments_ultimate,
                           average_solution_names,
                           DataOut<dim>::type_cell_data);

  data_out.build_patches();


  write_vtu_and_pvd<dim>(grid_pvdhandler,
                         data_out,
                         folder,
                         "overlaid_map",
                         time,
                         step_number,
                         1,
                         mpi_communicator);
}

template <int dim>
void
Compartmentalization<dim>::test()
{
  generate_cylindrical_grid();

  // This reads the electric field
  read_electric_field();

  // This agglomerate cells based on emw
  sort_agglomeration_deagglomeration_emw();

  // Here it should read the velocity magnitude data
  // and behave with each compartment as control volume and compartmentalize!
  read_velocity_magnitude();
  overlaid_map();
}

template class Compartmentalization<3>;
