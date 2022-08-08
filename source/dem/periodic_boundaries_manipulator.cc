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
PeriodicBoundariesManipulator<dim>::get_boundary_info(
  typename Triangulation<dim>::cell_iterator cell,
  unsigned int                               face_id,
  boundary_cells_info_struct<dim> &          boundary_information)
{
  // Initialize a simple quadrature for on the system. This will be used to
  // obtain a single sample point on the boundary faces.
  const FE_Q<dim>   fe(1);
  QGauss<dim - 1>   face_quadrature_formula(1);
  FEFaceValues<dim> fe_face_values(fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors);

  // Store information
  boundary_information.cell           = cell;
  boundary_information.boundary_id    = cell->face(face_id)->boundary_id();
  boundary_information.global_face_id = cell->face_index(face_id);

  // Finding the normal vector of the boundary face and point
  fe_face_values.reinit(cell, face_id);
  Tensor<1, dim> normal_vector = fe_face_values.normal_vector(0);
  Point<dim>     quad_point    = fe_face_values.quadrature_point(0);

  boundary_information.normal_vector = normal_vector;
  boundary_information.point_on_face = quad_point;
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::map_periodic_cells(
  const parallel::distributed::Triangulation<dim> &triangulation)
{
  periodic_boundary_cells_information.clear();

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
              if (face->boundary_id() == outlet_boundary_id)
                {
                  unsigned int face_id = cell->face_iterator_to_index(face);

                  // Save boundaries information related to the cell on
                  // outlet boundary.
                  boundary_cells_info_struct<dim> boundary_information;
                  get_boundary_info(cell, face_id, boundary_information);

                  // Save boundaries information related to the cell on
                  // Periodic boundary.
                  boundary_cells_info_struct<dim> periodic_boundary_information;
                  typename Triangulation<dim>::active_cell_iterator
                    periodic_cell = cell->periodic_neighbor(face_id);
                  get_boundary_info(periodic_cell,
                                    cell->periodic_neighbor_face_no(face_id),
                                    periodic_boundary_information);

                  // Store both cell information in map with cell id at
                  // outlet as key
                  periodic_boundary_cells_information.insert(
                    {boundary_information.cell->global_active_cell_index(),
                     std::make_pair(boundary_information,
                                    periodic_boundary_information)});
                }
            }
        }
    }
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::check_and_move_particles(
  boundary_cells_info_struct<dim> &cell_1,
  boundary_cells_info_struct<dim> &cell_2,
  typename Particles::ParticleHandler<dim>::particle_iterator_range
    &particles_in_cell)
{
  for (auto particle = particles_in_cell.begin();
       particle != particles_in_cell.end();
       ++particle)
    {
      Point<dim, double> particle_position = particle->get_location();

      // Calculate distance between particle position and the face of
      // boundary cell d = n•(pt_particle - pt_face)
      double distance_with_face =
        scalar_product(particle_position - cell_1.point_on_face,
                       cell_1.normal_vector);


      // If distance >= 0, particle is outside of cell (or on face).
      // If so, particle location is modified to get moved in the periodic cell.
      if (distance_with_face >= 0.0)
        {
          double distance_between_faces =
            cell_2.point_on_face[direction] - cell_1.point_on_face[direction];

          // Move particle outside the current cell to the periodic cell.
          particle_position[direction] += distance_between_faces;
          particle->set_location(particle_position);
        }
    }
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::execute_particles_displacement(
  const Particles::ParticleHandler<dim> &particle_handler)
{
  if (!periodic_boundary_cells_information.empty())
    {
      for (auto boundary_cells_information_iterator =
             periodic_boundary_cells_information.begin();
           boundary_cells_information_iterator !=
           periodic_boundary_cells_information.end();
           ++boundary_cells_information_iterator)
        {
          // Get the cell and periodic cell content from map
          auto boundary_cells_content =
            boundary_cells_information_iterator->second.first;
          auto cell = boundary_cells_content.cell;

          auto periodic_boundary_cells_content =
            boundary_cells_information_iterator->second.second;
          auto periodic_cell = periodic_boundary_cells_content.cell;

          // Since cell pairs are mapped once <cell at outlet, cell at PB>,
          // both cells in pair are checked individually, so we know the
          // cell the particle is coming from and will execute displacement
          // with the locally owned cell.
          if (cell->is_locally_owned())
            {
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                         particles_in_cell = particle_handler.particles_in_cell(cell);
              const bool particles_exist_in_main_cell =
                !particles_in_cell.empty();

              if (particles_exist_in_main_cell)
                {
                  check_and_move_particles(boundary_cells_content,
                                           periodic_boundary_cells_content,
                                           particles_in_cell);
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
                  check_and_move_particles(periodic_boundary_cells_content,
                                           boundary_cells_content,
                                           particles_in_periodic_cell);
                }
            }
        }
    }
}

template class PeriodicBoundariesManipulator<2>;
template class PeriodicBoundariesManipulator<3>;