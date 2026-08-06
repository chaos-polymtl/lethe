// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grid_cube_merged.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/fe/fe_q.h>

template <int dim, int spacedim>
GridCubeMerged<dim, spacedim>::GridCubeMerged(const std::string &grid_arguments)
{
  if constexpr (!(dim == 3 && spacedim == 3))
    {
      AssertThrow(
        false,
        ExcMessage(
          "Custom cube merged mesh is only supported in 3D space with 3D elements."));
      return;
    }

  this->grid_arguments = grid_arguments;

  // Separate arguments of the string
  std::vector<std::string> arguments;
  std::stringstream        s_stream(grid_arguments);
  while (s_stream.good())
    {
      std::string substr;
      getline(s_stream, substr, ':');
      arguments.push_back(substr);
    }

  // Arguments declaration
  if (arguments.size() != 4)
    {
      AssertThrow(
        false,
        ExcMessage(
          "Mandatory cylinder parameters are (plate file name : cylinder file name : cylindrical BID plate : cylindrical BID cylinder)"));
    }
  else
    {
      this->plate_file_name          = arguments[0];
      this->cylinder_file_name       = arguments[1];
      this->cylindrical_bid_plate    = Utilities::string_to_int(arguments[2]);
      this->cylindrical_bid_cylinder = Utilities::string_to_int(arguments[3]);
    }
}


template <>
void
GridCubeMerged<3, 3>::make_grid(Triangulation<3, 3> &triangulation)
{
  // Create a temporary triangulation for the extruded plate with hole
  Triangulation<3, 3> tria_plate;
  GridIn<3, 3>        grid_plate;
  grid_plate.attach_triangulation(tria_plate);
  std::ifstream plate_input_file(this->plate_file_name);
  grid_plate.read_msh(plate_input_file);

  // std::cout << "Read plate mesh" << std::endl;
  // std::cout << "Plate boundaries: " << tria_plate.get_boundary_ids().size() << std::endl;

  // Create a temporary triangulation for the inner cylinder
  Triangulation<3, 3> tria_cylinder;
  GridIn<3, 3>        grid_cylinder;
  grid_cylinder.attach_triangulation(tria_cylinder);
  std::ifstream cylinder_input_file(this->cylinder_file_name);
  grid_cylinder.read_msh(cylinder_input_file);
 
  // std::cout << "Read cylinder mesh" << std::endl;
  // std::cout << "Cylinder boundaries: " << tria_cylinder.get_boundary_ids().size() << std::endl;

  // Set cylindrical manifolds
  tria_plate.set_manifold(this->cylindrical_bid_plate,
    CylindricalManifold<3>(Tensor<1, 3>({0, 0, 1}),
    Point<3>({0, 0, 0})));

  for (const auto &cell : tria_plate.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary() && face->boundary_id() == this->cylindrical_bid_plate)
        {
          face->set_all_manifold_ids(this->cylindrical_bid_plate);
          cell->set_manifold_id(this->cylindrical_bid_plate);
        }
    
  tria_cylinder.set_manifold(this->cylindrical_bid_cylinder,
    CylindricalManifold<3>(Tensor<1, 3>({0, 0, 1}),
    Point<3>({0, 0, 0})));
  
  for (const auto &cell : tria_cylinder.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary() && face->boundary_id() == this->cylindrical_bid_cylinder)
        {
          face->set_all_manifold_ids(this->cylindrical_bid_cylinder);
          cell->set_manifold_id(this->cylindrical_bid_cylinder);
        }

  // Merge triangulations
  GridGenerator::merge_triangulations(
    tria_plate, tria_cylinder, triangulation, 1e-8, true, false);

  triangulation.set_manifold(this->cylindrical_bid_plate,
    CylindricalManifold<3>(Tensor<1, 3>({0, 0, 1}),
    Point<3>({0, 0, 0})));
  
  triangulation.set_manifold(this->cylindrical_bid_cylinder,
    CylindricalManifold<3>(Tensor<1, 3>({0, 0, 1}),
    Point<3>({0, 0, 0})));

  // std::cout << "Boundary ids: " << triangulation.get_boundary_ids().size() << std::endl;
  // for (const auto boundary_id : triangulation.get_boundary_ids())
  //   std::cout << boundary_id << ", ";
  
  // std::cout << std::endl;

  // std::cout << "Manifold ids: " << triangulation.get_manifold_ids().size() << std::endl;
  // for (const auto manifold_id : triangulation.get_manifold_ids())
  //   std::cout << manifold_id << ", ";
  
  // std::cout << std::endl;

  // Generate vtu file
  // DataOut<3>       data_out;
  // int mapping_degree = 1;
  // MappingQ<3, 3> mapping(mapping_degree);

  // DataOutBase::VtkFlags flags;
  // flags.write_higher_order_cells = true;
  // data_out.set_flags(flags);
  // data_out.attach_triangulation(triangulation);

  // Vector<double> ranks(triangulation.n_active_cells());
  // ranks = Utilities::MPI::this_mpi_process(triangulation.get_mpi_communicator());
  // data_out.add_data_vector(ranks, "ranks");
  // data_out.build_patches(mapping,
  //                        mapping_degree + 1,
  //                        DataOut<3>::CurvedCellRegion::curved_inner_cells);
  // data_out.write_vtu_in_parallel("out.vtu", triangulation.get_mpi_communicator());

  // Generate vtu file for boundaries
  // DataPostprocessors::BoundaryIds<3> boundary_ids;
  // DataOutFaces<3>                    data_out_faces;
  // FE_Q<3>                            dummy_fe(1);

  // DoFHandler<3> dummy_dof_handler(triangulation);
  // dummy_dof_handler.distribute_dofs(dummy_fe);

  // Vector<double> dummy_solution(dummy_dof_handler.n_dofs());

  // data_out_faces.attach_dof_handler(dummy_dof_handler);
  // data_out_faces.add_data_vector(dummy_solution, boundary_ids);
  // data_out_faces.build_patches();

  // std::ofstream out("boundary_ids.vtu");
  // data_out_faces.write_vtu(out);
}

// Fallback make_grid definition for unsupported template parameters. This
// provides a linker-visible symbol and a clear runtime error when the
// class is instantiated for dim/spacedim combinations that are not
// specialized above.
template <int dim, int spacedim>
void
GridCubeMerged<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> & /*triangulation*/)
{
  AssertThrow(
    false,
    ExcMessage(
      "GridCubeMerged is only implemented for dim = 3 and spacedim = 3."));
}

// Explicit template instantiations
template class GridCubeMerged<2, 2>;
template class GridCubeMerged<2, 3>;
template class GridCubeMerged<3, 3>;
