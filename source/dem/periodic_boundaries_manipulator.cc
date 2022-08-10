#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/periodic_boundaries_manipulator.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>

using namespace dealii;

template <int dim>
PeriodicBoundariesManipulator<dim>::PeriodicBoundariesManipulator()
{}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::get_periodic_boundaries_info(
  typename Triangulation<dim>::cell_iterator  cell,
  const unsigned int                          face_id,
  periodic_boundaries_cells_info_struct<dim> &boundaries_information)
{
  // Initialize a simple quadrature for the system. This will be used to
  // obtain a single sample point on the boundaries faces.
  const FE_Q<dim>   fe(1);
  QGauss<dim - 1>   face_quadrature_formula(1);
  FEFaceValues<dim> fe_face_values(fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors);
  Tensor<1, dim>    normal_vector;
  Point<dim>        quad_point;

  // Store information of the cell
  boundaries_information.cell           = cell;
  boundaries_information.boundary_id    = cell->face(face_id)->boundary_id();
  boundaries_information.global_face_id = cell->face_index(face_id);

  // Store information of the periodic cell
  boundaries_information.periodic_cell = cell->periodic_neighbor(face_id);
  boundaries_information.periodic_boundary_id =
    boundaries_information.periodic_cell->face(face_id)->boundary_id();
  boundaries_information.global_periodic_face_id =
    boundaries_information.periodic_cell->face_index(face_id);

  // Find the normal vector of the boundary face and point
  fe_face_values.reinit(cell, face_id);
  boundaries_information.normal_vector = fe_face_values.normal_vector(0);
  boundaries_information.point_on_face = fe_face_values.quadrature_point(0);

  // Find the normal vector of the periodic boundary face and point with the
  // reinit FE face values
  fe_face_values.reinit(boundaries_information.periodic_cell,
                        cell->periodic_neighbor_face_no(face_id));
  normal_vector = fe_face_values.normal_vector(0);
  quad_point    = fe_face_values.quadrature_point(0);
  boundaries_information.periodic_normal_vector = normal_vector;
  boundaries_information.point_on_periodic_face = quad_point;
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::map_periodic_cells(
  const parallel::distributed::Triangulation<dim> &triangulation)
{
  periodic_boundaries_cells_information.clear();

  // Iterating over the active cells in the trangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Iterating over cell faces
          for (const auto &face : cell->face_iterators())
            {
              // Check if face is on the periodic boundary flaged as outlet.
              // Pairs of periodic cells are stored once.
              for (unsigned int &outlet_boundary_id : outlet_boundary_ids)
                {
                  if (face->boundary_id() == outlet_boundary_id)
                    {
                      // Save boundaries information related to the cell on
                      // the outlet boundary of periodic walls.
                      // Information about both boundaries are stored in
                      // periodic_boundary_cells_info_struct
                      periodic_boundaries_cells_info_struct<dim>
                                   boundaries_information;
                      unsigned int face_id = cell->face_iterator_to_index(face);
                      get_periodic_boundaries_info(cell,
                                                   face_id,
                                                   boundaries_information);

                      // Store boundaries information in map with cell id at
                      // outlet as key
                      periodic_boundaries_cells_information.insert(
                        {boundaries_information.cell
                           ->global_active_cell_index(),
                         boundaries_information});
                    }
                }
            }
        }
    }
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::check_and_move_particles(
  const periodic_boundaries_cells_info_struct<dim> &boundaries_cells_content,
  typename Particles::ParticleHandler<dim>::particle_iterator_range
    &        particles_in_cell,
  const bool particles_in_outlet_cell)
{
  for (auto particle = particles_in_cell.begin();
       particle != particles_in_cell.end();
       ++particle)
    {
      // Get the current particle location
      Point<dim> particle_position = particle->get_location();

      // Initialize points and normal vector related of the cell that contains
      // the particle
      Point<dim>     point_on_face, point_on_periodic_face;
      Tensor<1, dim> normal_vector;
      if (particles_in_outlet_cell)
        {
          point_on_face = boundaries_cells_content.point_on_face;
          point_on_periodic_face =
            boundaries_cells_content.point_on_periodic_face;
          normal_vector = boundaries_cells_content.normal_vector;
        }
      else
        {
          point_on_face = boundaries_cells_content.point_on_periodic_face;
          point_on_periodic_face = boundaries_cells_content.point_on_face;
          normal_vector = boundaries_cells_content.periodic_normal_vector;
        }

      // Calculate distance between particle position and the face of
      // boundary cell d = n•(pt_particle - pt_face)
      double distance_with_face =
        scalar_product(particle_position - point_on_face, normal_vector);


      // If distance >= 0, particle is outside of cell (or on face).
      // If so, particle location is modified to get moved into the periodic
      // cell.
      if (distance_with_face >= 0.0)
        {
          Point<dim, double> distance_between_faces;
          distance_between_faces = point_on_periodic_face - point_on_face;

          // Move particle outside the current cell to the periodic cell.
          particle_position += distance_between_faces;
          particle->set_location(particle_position);
        }
    }
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::execute_particles_displacement(
  const Particles::ParticleHandler<dim> &particle_handler)
{
  if (!periodic_boundaries_cells_information.empty())
    {
      for (auto boundaries_cells_information_iterator =
             periodic_boundaries_cells_information.begin();
           boundaries_cells_information_iterator !=
           periodic_boundaries_cells_information.end();
           ++boundaries_cells_information_iterator)
        {
          // Gets the cell and periodic cell from map
          auto boundaries_cells_content =
            boundaries_cells_information_iterator->second;
          auto cell          = boundaries_cells_content.cell;
          auto periodic_cell = boundaries_cells_content.periodic_cell;

          // Checks and executes displacement of particles crossing a periodic
          // wall.
          if (cell->is_locally_owned())
            {
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                         particles_in_cell = particle_handler.particles_in_cell(cell);
              const bool particles_exist_in_main_cell =
                !particles_in_cell.empty();

              if (particles_exist_in_main_cell)
                {
                  bool particles_in_outlet_cell = true;
                  check_and_move_particles(boundaries_cells_content,
                                           particles_in_cell,
                                           particles_in_outlet_cell);
                }
            }

          if (periodic_cell->is_locally_owned())
            {
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_periodic_cell =
                  particle_handler.particles_in_cell(periodic_cell);
              const bool particles_exist_in_periodic_cell =
                !particles_in_periodic_cell.empty();

              if (particles_exist_in_periodic_cell)
                {
                  bool particles_in_outlet_cell = false;
                  check_and_move_particles(boundaries_cells_content,
                                           particles_in_periodic_cell,
                                           particles_in_outlet_cell);
                }
            }
        }
    }
}

template class PeriodicBoundariesManipulator<2>;
template class PeriodicBoundariesManipulator<3>;