/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */


#include <core/serial_solid.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/data_out.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <memory.h>

#include <fstream>



template <int dim, int spacedim>
SerialSolid<dim, spacedim>::SerialSolid(
  std::shared_ptr<Parameters::SolidObject<spacedim>> &param,
  unsigned int                                        id)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , param(param)
  , velocity(&param->solid_velocity)
  , id(id)
{
  if (param->solid_mesh.simplex)
    {
      // for simplex meshes
      fe            = std::make_shared<FE_SimplexP<dim, spacedim>>(1);
      solid_mapping = std::make_shared<MappingFE<dim, spacedim>>(*fe);
      solid_tria    = std::make_shared<Triangulation<dim, spacedim>>();

      displacement_fe =
        std::make_shared<FESystem<dim, spacedim>>(FE_SimplexP<dim, spacedim>(1),
                                                  spacedim);
    }
  else
    {
      // Usual case, for quad/hex meshes
      fe            = std::make_shared<FE_Q<dim, spacedim>>(1);
      solid_mapping = std::make_shared<MappingQGeneric<dim, spacedim>>(1);
      solid_tria    = std::make_shared<Triangulation<dim, spacedim>>(
        typename Triangulation<dim, spacedim>::MeshSmoothing(
          Triangulation<dim, spacedim>::smoothing_on_refinement |
          Triangulation<dim, spacedim>::smoothing_on_coarsening));

      displacement_fe =
        std::make_shared<FESystem<dim, spacedim>>(FE_Q<dim, spacedim>(1),
                                                  spacedim);
    }

  // Dof handler associated with particles
  solid_dh.clear();
  solid_dh.reinit(*this->solid_tria);

  // Dof handler associated with mesh displacement
  displacement_dh.clear();
  displacement_dh.reinit(*this->solid_tria);
}

template <int dim, int spacedim>
void
SerialSolid<dim, spacedim>::initial_setup()
{
  // initial_setup is called if the simulation is not a restarted one
  // Set-up of Nitsche triangulation, then particles (order important)
  setup_triangulation(false);
  setup_displacement();
}



template <int dim, int spacedim>
void
SerialSolid<dim, spacedim>::setup_triangulation(const bool restart)
{
  if (param->solid_mesh.type == Parameters::Mesh::Type::gmsh)
    {
      if (param->solid_mesh.simplex)
        {
          auto        comm      = solid_tria->get_communicator();
          std::string file_name = param->solid_mesh.file_name;

          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation_in_groups<dim, spacedim>(
              [file_name](dealii::Triangulation<dim, spacedim> &basetria) {
                GridIn<dim, spacedim> grid_in;
                grid_in.attach_triangulation(basetria);
                std::ifstream input_file(file_name);

                grid_in.read_msh(input_file);
              },
              [](dealii::Triangulation<dim, spacedim> &basetria,
                 const MPI_Comm                        comm,
                 const unsigned int /*group_size*/) {
                GridTools::partition_triangulation(
                  Utilities::MPI::n_mpi_processes(comm), basetria);
              },
              comm,
              Utilities::MPI::n_mpi_processes(comm) /* group size */,
              dealii::Triangulation<dim, spacedim>::none);

          solid_tria->create_triangulation(construction_data);
        }
      else
        { // Grid creation
          GridIn<dim, spacedim> grid_in;
          // Attach triangulation to the grid
          grid_in.attach_triangulation(*solid_tria);
          // Read input gmsh file
          std::ifstream input_file(param->solid_mesh.file_name);
          grid_in.read_msh(input_file);
        }
    }
  else if (param->solid_mesh.type == Parameters::Mesh::Type::dealii)
    {
      if (param->solid_mesh.simplex)
        { // TODO Using dealii generated meshes with simplices in Nitsche solver
          // generates and error. "boundary_id !=
          // numbers::internal_face_boundary_id"
          Triangulation<dim, spacedim> temporary_quad_triangulation;
          GridGenerator::generate_from_name_and_arguments(
            temporary_quad_triangulation,
            param->solid_mesh.grid_type,
            param->solid_mesh.grid_arguments);

          // initial refinement
          const int initial_refinement = param->solid_mesh.initial_refinement;
          temporary_quad_triangulation.refine_global(initial_refinement);
          // flatten the triangulation
          Triangulation<dim, spacedim> flat_temp_quad_triangulation;
          GridGenerator::flatten_triangulation(temporary_quad_triangulation,
                                               flat_temp_quad_triangulation);

          Triangulation<dim, spacedim> temporary_tri_triangulation(
            Triangulation<dim, spacedim>::limit_level_difference_at_vertices);
          GridGenerator::convert_hypercube_to_simplex_mesh(
            flat_temp_quad_triangulation, temporary_tri_triangulation);

          GridTools::partition_triangulation_zorder(
            Utilities::MPI::n_mpi_processes(solid_tria->get_communicator()),
            temporary_tri_triangulation);
          GridTools::partition_multigrid_levels(temporary_tri_triangulation);

          // extract relevant information from distributed triangulation
          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation(
              temporary_tri_triangulation,
              solid_tria->get_communicator(),
              TriangulationDescription::Settings::
                construct_multigrid_hierarchy);
          solid_tria->create_triangulation(construction_data);
        }
      else
        { // deal.ii creates grid and attaches the solid triangulation
          GridGenerator::generate_from_name_and_arguments(
            *solid_tria,
            param->solid_mesh.grid_type,
            param->solid_mesh.grid_arguments);
        }
    }
  else
    throw std::runtime_error(

      "Unsupported mesh type - solid mesh will not be created");

  // Translate the triangulation
  if (param->solid_mesh.translate)
    GridTools::shift(Point<spacedim>(param->solid_mesh.delta_x,
                                     param->solid_mesh.delta_y,
                                     param->solid_mesh.delta_z),
                     *solid_tria);
  if (param->solid_mesh.rotate)
    {
      rotate_grid(param->solid_mesh.angle, param->solid_mesh.axis);
    }

  // Refine the solid triangulation to its initial size
  // NB: solid_tria should not be refined if loaded from a restart file
  // afterwards
  if (!restart)
    {
      if (param->solid_mesh.simplex)
        {
          // Simplex triangulation refinement isn't possible yet
        }
      else
        {
          solid_tria->refine_global(param->solid_mesh.initial_refinement);
        }
      solid_dh.distribute_dofs(*fe);
    }
}

template <int dim, int spacedim>
void
SerialSolid<dim, spacedim>::load_triangulation(const std::string filename_tria)
{
  // Load solid triangulation from given file
  // TODO not functional for now, as not all information as passed with the load
  // function (see dealii documentation for classTriangulation) => change the
  // way the load works, or change the way the solid triangulation is handled
  std::ifstream in_folder(filename_tria.c_str());
  if (!in_folder)
    AssertThrow(false,
                ExcMessage(
                  std::string(
                    "You are trying to restart a previous computation, "
                    "but the restart file <") +
                  filename_tria + "> does not appear to exist!"));

  try
    {
      if (auto solid_tria =
            dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
              get_solid_triangulation().get()))
        {
          solid_tria->load(filename_tria.c_str());
        }
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("Cannot open snapshot mesh file or read the "
                             "solid triangulation stored there."));
    }

  // Initialize dof handler for solid
  solid_dh.distribute_dofs(*fe);
}

template <int dim, int spacedim>
DoFHandler<dim, spacedim> &
SerialSolid<dim, spacedim>::get_solid_dof_handler()
{
  return solid_dh;
}

template <int dim, int spacedim>
DoFHandler<dim, spacedim> &
SerialSolid<dim, spacedim>::get_displacement_dof_handler()
{
  return displacement_dh;
}


template <int dim, int spacedim>
Vector<double> &
SerialSolid<dim, spacedim>::get_displacement_vector()
{
  return displacement;
}

template <int dim, int spacedim>
std::shared_ptr<Triangulation<dim, spacedim>>
SerialSolid<dim, spacedim>::get_solid_triangulation()
{
  return solid_tria;
}

template <int dim, int spacedim>
Function<spacedim> *
SerialSolid<dim, spacedim>::get_solid_velocity()
{
  return velocity;
}


template <>
void
SerialSolid<1, 2>::rotate_grid(double /*angle*/, int /*axis*/)
{
  // Not implemented right now
  throw(std::runtime_error("This is currently not implemented"));
}

template <>
void
SerialSolid<2, 2>::rotate_grid(double /*angle*/, int /*axis*/)
{
  // Not implemented right now
  throw(std::runtime_error("This is currently not implemented"));
}
template <>
void
SerialSolid<3, 3>::rotate_grid(double /*angle*/, int /*axis*/)
{
  // Not implemented right now
  throw(std::runtime_error("This is currently not implemented"));
}

template <>
void
SerialSolid<2, 3>::rotate_grid(double angle, int axis)
{
  GridTools::rotate(angle, axis, *solid_tria);
}


template <int dim, int spacedim>
unsigned int
SerialSolid<dim, spacedim>::get_solid_id()
{
  return id;
}

template <int dim, int spacedim>
void
SerialSolid<dim, spacedim>::move_solid_triangulation(const double time_step,
                                                     const double initial_time)
{
  // First we calculate the displacement that the solid triangulation
  // will feel during this iteration.
  // Then, we apply this displacement to the triangulation
  const unsigned int      dofs_per_cell = displacement_fe->dofs_per_cell;
  const std::vector<bool> locally_owned_vertices =
    GridTools::get_locally_owned_vertices(*this->solid_tria);
  std::vector<bool> vertex_moved(this->solid_tria->n_vertices(), false);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


  for (const auto &cell : displacement_dh.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int vertex = 0;
               vertex < cell->reference_cell().n_vertices();
               ++vertex)
            {
              // Calculate the next position of the particle using the RK4
              // time integration scheme
              // x(t+dt) = x(t) + 1/6 (k1+2*k2+2*k3+k4)
              // k1 = dt * v(t,x(t))
              // k2 = dt * v(t+0.5*dt,x+k1/2)
              // k3 = dt * v(t+0.5*dt,x+k2/2)
              // k4 = dt * vt(t+dt,x+k3)
              // The four stages (1 to 4) are built successively
              // For each stage, the time of the function must be "reset"
              // to ensure adequate evaluation of the RK4 scheme.
              // See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
              // for more details

              const auto       dof_index = cell->vertex_dof_index(vertex, 0);
              Point<spacedim> &vertex_position  = cell->vertex(vertex);
              const unsigned   global_vertex_no = cell->vertex_index(vertex);
              if (vertex_moved[global_vertex_no] ||
                  !locally_owned_vertices[global_vertex_no])
                continue;

              Tensor<1, spacedim> k1;
              velocity->set_time(initial_time);
              for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                k1[comp_i] = velocity->value(vertex_position, comp_i);

              Point<spacedim>     p1 = vertex_position + time_step / 2 * k1;
              Tensor<1, spacedim> k2;
              velocity->set_time(initial_time + time_step / 2);
              for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                k2[comp_i] = velocity->value(p1, comp_i);

              Point<spacedim>     p2 = vertex_position + time_step / 2 * k2;
              Tensor<1, spacedim> k3;
              for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                k3[comp_i] = velocity->value(p2, comp_i);

              Point<spacedim>     p3 = vertex_position + time_step * k3;
              Tensor<1, spacedim> k4;
              velocity->set_time(initial_time + time_step);
              for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                k4[comp_i] = velocity->value(p3, comp_i);

              auto vertex_displacement =
                time_step / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

              vertex_position += vertex_displacement;

              for (unsigned d = 0; d < spacedim; ++d)
                displacement[dof_index + d] =
                  displacement[dof_index + d] + vertex_displacement[d];

              vertex_moved[global_vertex_no] = true;
            }
        }
    }
}


template <int dim, int spacedim>
void
SerialSolid<dim, spacedim>::move_solid_triangulation_with_displacement()
{
  const unsigned int     dofs_per_cell = displacement_fe->dofs_per_cell;
  std::set<unsigned int> dof_vertex_displaced;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : displacement_dh.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int vertex = 0;
               vertex < cell->reference_cell().n_vertices();
               ++vertex)
            {
              if (!dof_vertex_displaced.count(cell->vertex_index(vertex)))
                {
                  const auto dof_index = cell->vertex_dof_index(vertex, 0);
                  Point<spacedim> &vertex_position = cell->vertex(vertex);

                  for (unsigned d = 0; d < spacedim; ++d)
                    {
                      vertex_position[d] =
                        vertex_position[d] + displacement[dof_index + d];
                    }

                  dof_vertex_displaced.insert(cell->vertex_index(vertex));
                }
            }
        }
    }
}


template <int dim, int spacedim>
void SerialSolid<dim, spacedim>::write_checkpoint(std::string /*prefix*/)
{
  // SolutionTransfer<dim, Vector<double>, spacedim> system_trans_vectors(
  //   this->displacement_dh);

  // std::vector<const Vector<double> *> sol_set_transfer;
  // sol_set_transfer.push_back(&displacement);

  // system_trans_vectors.prepare_for_serialization(sol_set_transfer);

  // if (auto tria = dynamic_cast<parallel::distributed::Triangulation<dim> *>(
  //       this->solid_tria.get()))
  //   {
  //     std::string triangulationName = prefix + ".triangulation";
  //     tria->save(prefix + ".triangulation");
  //   }
}

template <int dim, int spacedim>
void SerialSolid<dim, spacedim>::read_checkpoint(std::string /*prefix*/)
{
  // Setup an un-refined triangulation before loading
  // setup_triangulation(true);

  // Read the triangulation from the checkpoint
  // const std::string filename = prefix + ".triangulation";
  // std::ifstream     in(filename.c_str());
  // if (!in)
  //  AssertThrow(false,
  //              ExcMessage(
  //                std::string("You are trying to read a solid triangulation, "
  //                            "but the restart file <") +
  //                filename + "> does not appear to exist!"));

  // try
  //   {
  //     if (auto tria = dynamic_cast<parallel::distributed::Triangulation<dim>
  //     *>(
  //           this->solid_tria.get()))
  //       tria->load(filename.c_str());
  //   }
  // catch (...)
  //   {
  //     AssertThrow(false,
  //                 ExcMessage("Cannot open snapshot mesh file or read the "
  //                            "triangulation stored there."));
  //   }

  // Setup dof-handler for solid and displacement
  // solid_dh.distribute_dofs(*fe);
  // setup_displacement();


  //// Read displacement vector
  // std::vector<Vector<double> *> x_system(1);
  // x_system[0] = &(displacement);

  // SolutionTransfer<dim, Vector<double>, spacedim> system_trans_vectors(
  //   this->displacement_dh);

  // system_trans_vectors.deserialize(x_system);
  // displacement;

  //// Reset triangulation position using displacement vector
  // move_solid_triangulation_with_displacement();


  // We did not checkpoint particles, we re-create them from scratch
}


template <int dim, int spacedim>
void
SerialSolid<dim, spacedim>::setup_displacement()
{
  displacement_dh.distribute_dofs(*this->displacement_fe);

  displacement.reinit(displacement_dh.n_locally_owned_dofs());

  displacement = 0;
}


// Pre-compile the 2D, 3D and the 2D in 3D versions with the types that can
// occur
template class SerialSolid<1, 2>;
template class SerialSolid<2, 3>;
template class SerialSolid<2>;
template class SerialSolid<3>;
