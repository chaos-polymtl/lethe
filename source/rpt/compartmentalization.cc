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
{}

template <int dim>
void
Compartmentalization<dim>::generate_cylindrical_grid()
{
  GridOut grid_out;
  GridGenerator::subdivided_cylinder(
    triangulation,
    cp_parameters.cp_param.subdivisions,
    cp_parameters.cp_param.cylinder_radius,
    cp_parameters.cp_param.cylinder_half_length);
  triangulation.refine_global(cp_parameters.cp_param.initial_refinement);
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

  // output the center of cell to feed the Pyvista and extract the electric
  // field and velocity at the cell centers
  // from COMSOL or CFD simulation
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
Compartmentalization<dim>::electric_field_stored_in_vector()
  -> std::vector<std::vector<double>>
{
  std::vector<double> electric_field;
  std::ifstream       fin("input_values_emw");
  double              element;
  while (fin >> element)
    {
      electric_field.push_back(element);
    }


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
Compartmentalization<dim>::electric_field_stored_in_map()
  -> std::map<typename Triangulation<dim>::active_cell_iterator, double>
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
      primer_map_cell_emw.insert({cell, electric_field[cell_index]});
    }
  return primer_map_cell_emw;
}

template <int dim>
auto
Compartmentalization<dim>::read_velocity_vector()
  -> std::map<typename Triangulation<dim>::active_cell_iterator,
              std::vector<double>>
{
  std::ifstream       inputFile{"velocity_vector"};
  std::vector<double> velocity_components{
    std::istream_iterator<double>{inputFile}, {}};
  for (unsigned int i = 0; i < cells_and_indices.size(); i++)
    {
      std::vector<double> velocity_vector;
      for (unsigned int j = 0; j < 3; j++)
        {
          velocity_vector.push_back(velocity_components[0]);
          velocity_components.erase(velocity_components.begin());
        }
      // Push {cell,[vector of velocity]}
      // in cells_and_indices[i], "i" is the key which is cell index
      // cells_and_indices[i] in the world of C++ map, means the element which
      // its key is "i"
      map_of_all_cells_and_velocity_vectors.insert(
        {cells_and_indices[i], velocity_vector});
      velocity_vector.clear();
    }
  return map_of_all_cells_and_velocity_vectors;
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

  //  this for loop fill the container with the properties
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
  // Store the sorted cells (active_cell_iterator)
  std::vector<typename Triangulation<dim>::active_cell_iterator> vector_cell;

  // Sort cells based on their value of the electric fields and the
  // cell index associated with the cell will move along with value
  std::sort(primer_matrix_index_emw.begin(),
            primer_matrix_index_emw.end(),
            [](const std::vector<double> &a, const std::vector<double> &b) {
              return a[1] > b[1];
            });

  // All the cells are sorted based on their values
  // (a[0]=cell_index and a[1]=value of the electric field)
  // The following "for loop" rearrange the cells (active_cell_iterator) to have
  // the sorted version it uses the sorted cell_index to rearrange and sorts the
  // cells
  for (auto &i : primer_matrix_index_emw)
    {
      double     index  = i[0];
      const auto search = cells_and_indices.find(index);
      // vector_cell stores only the sorted cells
      vector_cell.push_back(search->second);
    }

  // Since we have at least one compartment from the beginning, the iterator
  // starts from 1
  int number_of_compartments = 1;

  // this "k" iterator will avoid left over the last cell in case
  // that it must be clustered in one separate group
  unsigned int k = 1;

  // first_set, is a map containing the compartment_ids and cells inside each
  // compartment before the deagglomeration step since it is possible that some
  // cells locate in same compartment but, they are not connected directly or
  // indirectly, the agglomeration is essential
  std::map<int, std::vector<typename Triangulation<dim>::active_cell_iterator>>
    first_set;

  // vector will temporarily contain the cell of each compartment
  std::vector<typename Triangulation<dim>::active_cell_iterator> vector;

  // tol is the threshold of the clustering for difference between the
  // electrical field
  double Tol = cp_parameters.cp_param.electric_field_tolerance;

  // the first cell in the sorted cells goes to first compartment
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
      k++;
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

  //  Start the deagglomeration (compartment_id)
  int key_final_set = 1;

  // This is a vector that contains those cells that are connected in the first
  // cluster
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
          // since it is part of unique cluster
          filling_vector_declustered_set.push_back(*set_to_search.begin());

          // now we should delete it from the set_to_search since it got its
          // forever home
          set_to_search.erase(*set_to_search.begin());

          // at first, search for the neighbour of this first cell that we
          // have added manually
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
          key_final_set++;
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
      compartment_final_id++;
    }
  write_file_compartments_first_field(compartments_emw, 1.0, 1);

  return final_set;
}

template <int dim>
void
Compartmentalization<dim>::overlaid_map()
{
  // loop over each compartment
  // loop over cells inside each compartment
  // sort the cells based on their velocity value
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
          vector_cell.push_back(primer_map_cell_index.find(index)->second);
        }
      int          number_of_compartments = 1;
      unsigned int k                      = 1;
      std::map<int,
               std::vector<typename Triangulation<dim>::active_cell_iterator>>
        first_set;

      // vector will temporarily contain the cell of each compartment
      std::vector<typename Triangulation<dim>::active_cell_iterator> vector;
      double Tol = cp_parameters.cp_param.velocity_tolerance;
      vector.push_back(vector_cell[0]);
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
      first_set.clear();
      vector.clear();
      vector_cell.clear();
      filling_vector_declustered_set.clear();
      set_to_search.clear();
      temporary_set.clear();
      inner_matrix.clear();
      primer_map_cell_index.clear();
    }

  std::vector<double> electric_field_average;
  for (const auto &pair : overlaid_set)
    {
      // loop over cells
      double sum = 0;
      for (typename Triangulation<dim>::active_cell_iterator d : pair.second)
        {
          sum = sum + primer_map_cell_emw[d];
        }
      double average = sum / pair.second.size();
      electric_field_average.push_back(average);
    }

  std::string   filename = "Average Electric Field";
  std::ofstream myfile;
  std::string   sep = " ";
  myfile.open(filename);
  for (auto &EF : electric_field_average)
    {
      myfile << EF << "\n";
    }
  myfile.close();
  std::vector<double> volume;
  // loop over compartments
  for (const auto &pair : overlaid_set)
    {
      // loop over cells
      double sum = 0;
      for (typename Triangulation<dim>::active_cell_iterator d : pair.second)
        {
          sum = sum + d->measure();
        }
      volume.push_back(sum);
    }

  std::string   file_name = "volume";
  std::ofstream my_file;
  std::string   sep_ = " ";
  myfile.open(file_name);
  for (auto &EF : volume)
    {
      myfile << EF << "\n";
    }
  my_file.close();
  double         compartment_final_id = 1;
  Vector<double> compartments_ultimate;
  compartments_ultimate.reinit(triangulation.n_active_cells());
  for (const auto &pair : overlaid_set)
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
Compartmentalization<dim>::flux_calculation_between_compartments()
{
  std::vector<double> neighbour_velocity;
  std::vector<double> main_velocity;
  std::vector<double> average_velocity_vector;
  double              average_velocity;

  // the shared surface between two cells
  double area;

  // the volumetric flux between two cells
  double volumetric_flux;

  // to clarify the direction (+ or -)
  double      This;
  std::string in;
  std::string out;

  //"map_of_all_sets_and_cells" helps to search and find the
  // compartment_id associated to each cell
  //"map_of_all_sets_and_cells" is a map that has the cell
  //(active_cell_iterator) as the key and the associated compartment_id
  std::map<typename Triangulation<dim>::active_cell_iterator, int>
    map_of_all_sets_and_cells;

  // filling the map
  for (const auto &pair : overlaid_set)
    {
      for (const typename Triangulation<dim>::active_cell_iterator &k :
           pair.second)
        {
          map_of_all_sets_and_cells.insert({k, pair.first});
        }
    }
  // loop over compartments
  for (const auto &pair : overlaid_set)
    {
      // define a set containing cell iterator as key element and the associated
      // physical value
      std::set<typename Triangulation<dim>::active_cell_iterator> set_to_search;
      for (const typename Triangulation<dim>::active_cell_iterator &k :
           pair.second)
        {
          set_to_search.insert(k);
        }
      // loop over cells inside each compartment
      for (typename Triangulation<dim>::active_cell_iterator d : pair.second)
        {
          // main cell is the one that we are looping over it now
          Point<3> center_of_main_cell;
          Point<3> center_of_neighbor;
          center_of_main_cell = d->center();
          for (unsigned int i = 0; i < 6; i++)
            {
              std::vector<double> norm_of_cell;
              typename Triangulation<dim>::active_cell_iterator
                neighbours_of_each_cell = d->neighbor(i);

              // this means if the neighbour of the cell is not in the same
              // cluster
              if (set_to_search.count(neighbours_of_each_cell) != 1)
                {
                  // here it find the cluster that this neighbour belongs to
                  for (auto &it : map_of_all_sets_and_cells)
                    {
                      if (it.first == neighbours_of_each_cell)
                        {
                          center_of_neighbor =
                            neighbours_of_each_cell->center();
                          area = neighbours_of_each_cell->face(i)->measure();
                          std::string compartment_id_second =
                            std::to_string(it.second);
                          std::string compartment_id_first =
                            std::to_string(pair.first);
                          out =
                            compartment_id_first + "-" + compartment_id_second;
                          in =
                            compartment_id_second + "-" + compartment_id_first;
                          for (auto &s : map_of_all_cells_and_velocity_vectors)
                            {
                              if (s.first == neighbours_of_each_cell)
                                {
                                  neighbour_velocity = s.second;
                                }
                            }
                          for (auto &ss : map_of_all_cells_and_velocity_vectors)
                            {
                              if (ss.first == d)
                                {
                                  main_velocity = ss.second;
                                }
                            }
                          // calculate normal vector of the cell
                          norm_of_cell.push_back(
                            (center_of_main_cell[0] - center_of_neighbor[0]) /
                            ((center_of_main_cell - center_of_neighbor)
                               .norm()));
                          norm_of_cell.push_back(
                            (center_of_main_cell[1] - center_of_neighbor[1]) /
                            ((center_of_main_cell - center_of_neighbor)
                               .norm()));
                          norm_of_cell.push_back(
                            (center_of_main_cell[2] - center_of_neighbor[2]) /
                            ((center_of_main_cell - center_of_neighbor)
                               .norm()));

                          // calculate average velocity at the common face of
                          // two cells
                          average_velocity_vector.push_back(
                            (neighbour_velocity[0] + main_velocity[0]) / 2);
                          average_velocity_vector.push_back(
                            (neighbour_velocity[1] + main_velocity[1]) / 2);
                          average_velocity_vector.push_back(
                            (neighbour_velocity[2] + main_velocity[2]) / 2);

                          // calculate dot product to see the direction of
                          // velocity vector and its value
                          This = norm_of_cell[0] * average_velocity_vector[0] +
                                 norm_of_cell[1] * average_velocity_vector[1] +
                                 norm_of_cell[2] * average_velocity_vector[2];
                          average_velocity = abs(This);
                          volumetric_flux  = average_velocity * area;

                          // the volumetric_flux/2 (the divided by 2 is because
                          // it counts each two neighbour twice)
                          if (flux.size() == 0)
                            {
                              if (This > 0)
                                {
                                  flux.insert({in, volumetric_flux / 2});
                                }
                              else
                                flux.insert({out, volumetric_flux / 2});
                            }
                          else if (This > 0)
                            {
                              if (flux.find(in) != flux.end())
                                {
                                  flux[in] = flux[in] + volumetric_flux / 2;
                                }
                              else
                                flux.insert({in, volumetric_flux / 2});
                            }
                          else if (flux.find(out) != flux.end())
                            {
                              flux[out] = flux[out] + volumetric_flux / 2;
                            }
                          else
                            flux.insert({out, volumetric_flux / 2});

                          center_of_neighbor.clear();
                          norm_of_cell.clear();
                          average_velocity_vector.clear();
                          neighbour_velocity.clear();
                        }
                    }
                }
            }
          center_of_main_cell.clear();
          main_velocity.clear();
        }
    }
}

template <int dim>
void
Compartmentalization<dim>::inlet_outlet_flux()
{
  std::vector<double> velocity_vector_outlet;
  std::vector<double> velocity_vector_inlet;
  std::vector<double> velocity_vector_at_boundary;
  // loop over compartments
  // loop over cell
  // loop over face

  // loop over compartments
  for (const auto &pair : overlaid_set)
    {
      std::string compartment_at_boundary_id = std::to_string(pair.first);
      std::string in_flow_name = "Inlet flow to " + compartment_at_boundary_id;
      std::string out_flow_name =
        "outlet flow from " + compartment_at_boundary_id;
      double in_flow        = 0;
      double out_flow       = 0;
      double CFD_inlet_flow = cp_parameters.cp_param.CFD_input_velocity;

      // loop over cells inside each compartment
      for (typename Triangulation<dim>::active_cell_iterator d : pair.second)
        {
          // loop over face
          for (unsigned int i = 0; i < 6; i++)
            {
              if (d->face(i)->at_boundary())
                {
                  if (d->face(i)->boundary_id() == 1)
                    {
                      for (unsigned int j = 0; j < 6; j++)
                        {
                          if (d->face(j)->boundary_id() == 0)
                            {
                              in_flow =
                                in_flow + 0; // very close to the wall, noslip;
                            }
                          else
                            in_flow =
                              in_flow +
                              (d->face(i)->measure() *
                               CFD_inlet_flow); // far from wall at the boundary
                          break;
                        }
                    }
                  if (d->face(i)->boundary_id() == 2)
                    {
                      for (auto &s : map_of_all_cells_and_velocity_vectors)
                        {
                          if (s.first == d)
                            {
                              velocity_vector_outlet = s.second;
                            }
                        }
                      // source of error can be here since we are considering
                      // the velocity at the cell center and not the face
                      out_flow =
                        out_flow + (d->face(i)->measure() *
                                    (pow((pow(velocity_vector_outlet[0], 2)) +
                                           (pow(velocity_vector_outlet[1], 2)) +
                                           (pow(velocity_vector_outlet[2], 2)),
                                         0.5)));
                    }
                }
            }
        }
      inlet_flux.insert({in_flow_name, in_flow});
      outlet_flux.insert({out_flow_name, out_flow});
    }
}


template <int dim>
void
Compartmentalization<dim>::output_flux()
{
  // This is a matrix including all fluxes between compartments
  std::vector<std::vector<double>> matrix_of_flux;
  // To fill the matrix
  std::vector<double> filling_matrix;
  // combination includes to connected compartments in form of string
  // fill the filling_matrix with input fluxes

  filling_matrix.push_back(0);
  for (auto &in : inlet_flux)
    {
      filling_matrix.push_back(in.second);
    }
  filling_matrix.push_back(0);
  // push it to matrix_of_flux
  matrix_of_flux.push_back(filling_matrix);
  filling_matrix.clear();

  std::string combination;
  std::string outlet_combination;
  // loop over flux
  for (unsigned int i = 1; i < overlaid_set.size() + 1; i++)
    {
      for (unsigned int j = 1; j < overlaid_set.size() + 1; j++)
        {
          std::string first  = std::to_string(i);
          std::string second = std::to_string(j);
          combination        = first + "-" + second;
          if (flux.find(combination) == flux.end())
            {
              filling_matrix.push_back(0);
            }
          else
            filling_matrix.push_back(flux.find(combination)->second);
        }
      std::string compartment_id = std::to_string(i);
      outlet_combination         = "outlet flow from " + compartment_id;
      filling_matrix.push_back(outlet_flux.find(outlet_combination)->second);
      matrix_of_flux.push_back(filling_matrix);
      filling_matrix.clear();
    }


  std::string   filename = "Flux_Matrix";
  std::ofstream myfile;
  std::string   sep = " ";
  myfile.open(filename);
  for (auto &vector : matrix_of_flux)
    {
      for (unsigned int j = 0; j < overlaid_set.size() + 1; j++)
        {
          myfile << vector[j] << sep;
        }
      myfile << "\n";
    }
  myfile.close();
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
  electric_field_stored_in_map();
  electric_field_stored_in_vector();
  sort_agglomeration_deagglomeration_emw();
  read_velocity_magnitude();
  read_velocity_vector();
  overlaid_map();
  flux_calculation_between_compartments();
  inlet_outlet_flux();
  output_flux();
}

template class Compartmentalization<3>;
