#include <dem/adaptive_sparse_contacts.h>
#include <dem/load_balancing.h>

template <int dim>
LoadBalancing<dim>::LoadBalancing()
{}

template <int dim>
unsigned int
LoadBalancing<dim>::calculate_total_cell_weight(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const CellStatus                       status,
  const Particles::ParticleHandler<dim> &particle_handler) const
{
  // Assign no weight to cells we do not own.
  if (!cell->is_locally_owned())
    return 0;

  switch (status)
    {
      case CellStatus::cell_will_persist:

      case CellStatus::cell_will_be_refined:
        // If CELL_PERSIST, do as CELL_REFINE
        {
          const unsigned int n_particles_in_cell =
            particle_handler.n_particles_in_cell(cell);
          return n_particles_in_cell * particle_weight;
          break;
        }
      case CellStatus::cell_invalid:
        break;

      case CellStatus::children_will_be_coarsened:
        {
          unsigned int n_particles_in_cell = 0;

          for (unsigned int child_index = 0;
               child_index < GeometryInfo<dim>::max_children_per_cell;
               ++child_index)
            n_particles_in_cell +=
              particle_handler.n_particles_in_cell(cell->child(child_index));

          return n_particles_in_cell * particle_weight;
          break;
        }

      default:
        Assert(false, ExcInternalError());
        break;
    }

  return 0;
}

template <int dim>
unsigned int
LoadBalancing<dim>::calculate_total_cell_weight_with_mobility_status(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const CellStatus                       status,
  const Particles::ParticleHandler<dim> &particle_handler,
  const typename DEM::dem_data_structures<dim>::cell_index_int_map
    mobility_status) const
{
  // Assign no weight to cells we do not own.
  if (!cell->is_locally_owned())
    return 0;

  // Get mobility status of the cell
  const unsigned int cell_mobility_status = mobility_status.at(cell->active_cell_index());

  // Applied a factor on the particle weight regards the mobility status
  // Factor of 1 when mobile cell
  double alpha = 1.0;
  if (cell_mobility_status == AdaptiveSparseContacts<dim>::static_active ||
      cell_mobility_status == AdaptiveSparseContacts<dim>::advected_active)
    {
      alpha = active_load_balancing_factor;
    }
  else if (cell_mobility_status == AdaptiveSparseContacts<dim>::inactive ||
           cell_mobility_status == AdaptiveSparseContacts<dim>::advected)
    {
      alpha = inactive_load_balancing_factor;
    }

  switch (status)
    {
      case dealii::CellStatus::cell_will_persist:

      case dealii::CellStatus::cell_will_be_refined:
        {
          const unsigned int n_particles_in_cell =
            particle_handler.n_particles_in_cell(cell);
          return alpha * n_particles_in_cell * particle_weight;
          break;
        }

      case dealii::CellStatus::cell_invalid:
        break;

      case dealii::CellStatus::children_will_be_coarsened:
        {
          unsigned int n_particles_in_cell = 0;

          for (unsigned int child_index = 0;
               child_index < GeometryInfo<dim>::max_children_per_cell;
               ++child_index)
            n_particles_in_cell +=
              particle_handler.n_particles_in_cell(cell->child(child_index));

          return alpha * n_particles_in_cell * particle_weight;
          break;
        }

      default:
        Assert(false, ExcInternalError());
        break;
    }

  return 0;
}

template class LoadBalancing<2>;
template class LoadBalancing<3>;