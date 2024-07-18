#include <dem/adaptive_sparse_contacts.h>
#include <dem/load_balancing.h>

using namespace dealii;

template <int dim>
LoadBalancing<dim>::LoadBalancing()
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
{}

template <int dim>
inline bool
LoadBalancing<dim>::check_load_balance_once()
{
  if (simulation_control->get_step_number() == load_balance_step)
    {
      return true;
    }

  return false;
}

template <int dim>
inline bool
LoadBalancing<dim>::check_load_balance_frequent()
{
  if ((simulation_control->get_step_number() % load_balance_frequency) == 0)
    {
      return true;
    }

  return false;
}

template <int dim>
inline bool
LoadBalancing<dim>::check_load_balance_dynamic()
{
  if (simulation_control->get_step_number() % dynamic_check_frequency == 0)
    {
      unsigned int maximum_particle_number_on_proc = 0;
      unsigned int minimum_particle_number_on_proc = 0;

      maximum_particle_number_on_proc =
        Utilities::MPI::max(particle_handler->n_locally_owned_particles(),
                            mpi_communicator);
      minimum_particle_number_on_proc =
        Utilities::MPI::min(particle_handler->n_locally_owned_particles(),
                            mpi_communicator);

      // Execute load balancing if difference of load between processors is
      // larger than threshold of the load per processor
      if ((maximum_particle_number_on_proc - minimum_particle_number_on_proc) >
          load_threshold *
            (particle_handler->n_global_particles() / n_mpi_processes))
        {
          return true;
        }
    }

  return false;
}

template <int dim>
inline bool
LoadBalancing<dim>::check_load_balance_with_sparse_contacts()
{
  if (simulation_control->get_step_number() % dynamic_check_frequency == 0)
    {
      // Process to accumulate the load of each process regards the number
      // of cells and particles with their selected weight and with a factor
      // related to the mobility status of the cells
      std::vector<double> process_to_load_weight(n_mpi_processes, 0.0);


      for (const auto &cell : triangulation->active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              // Apply a weight of 1000 to the cell (default value)
              process_to_load_weight[this_mpi_process] += 1000;

              // Get the mobility status of the cell & the number of particles
              const unsigned int cell_mobility_status =
                adaptive_sparse_contacts->check_cell_mobility(cell);
              const unsigned int n_particles_in_cell =
                particle_handler->n_particles_in_cell(cell);

              // Apply a factor on the particle weight regards the
              // mobility status. alpha = 1 by default for mobile cell, but
              // is modified if cell is active or inactive
              double alpha = 1.0;
              if (cell_mobility_status ==
                    AdaptiveSparseContacts<dim>::static_active ||
                  cell_mobility_status ==
                    AdaptiveSparseContacts<dim>::advected_active)
                {
                  alpha = active_status_factor;
                }
              else if (cell_mobility_status ==
                         AdaptiveSparseContacts<dim>::inactive ||
                       cell_mobility_status ==
                         AdaptiveSparseContacts<dim>::advected)
                {
                  alpha = inactive_status_factor;
                }

              // Add the particle weight time the number of particles in the
              // cell to the processor load
              process_to_load_weight[this_mpi_process] +=
                alpha * n_particles_in_cell * particle_weight;
            }
        }

      // Exchange information
      double maximum_load_on_proc = 0.0;
      double minimum_load_on_proc = 0.0;
      double total_load           = 0.0;

      maximum_load_on_proc =
        Utilities::MPI::max(*std::max_element(process_to_load_weight.begin(),
                                              process_to_load_weight.end()),
                            mpi_communicator);

      // Find the minimum load on a process
      // First it finds the minimum load on a process, but since values in
      // the vector that are not on this process are 0.0, it looks for
      // values > 1e-8. After that, it finds the minimum load of all the
      // processors
      minimum_load_on_proc =
        Utilities::MPI::min(*std::min_element(process_to_load_weight.begin(),
                                              process_to_load_weight.end(),
                                              [](double a, double b) {
                                                return (a > 1e-8) ?
                                                         (b > 1e-8 ? a < b :
                                                                     true) :
                                                         false;
                                              }),
                            mpi_communicator);

      // Get the total load
      total_load =
        Utilities::MPI::sum(std::accumulate(process_to_load_weight.begin(),
                                            process_to_load_weight.end(),
                                            0.0),
                            mpi_communicator);

      if ((maximum_load_on_proc - minimum_load_on_proc) >
          load_threshold * (total_load / n_mpi_processes))
        {
          return true;
        }
    }

  // Clear and connect a new cell weight function
  connect_mobility_status_weight_signals();

  return false;
}

template <int dim>
unsigned int
LoadBalancing<dim>::calculate_total_cell_weight(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const CellStatus status) const
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
            particle_handler->n_particles_in_cell(cell);
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
              particle_handler->n_particles_in_cell(cell->child(child_index));

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
  const CellStatus status) const
{
  // Assign no weight to cells we do not own.
  if (!cell->is_locally_owned())
    return 0;

  // Get mobility status of the cell
  const unsigned int cell_mobility_status =
    adaptive_sparse_contacts->check_cell_mobility(cell);

  // Applied a factor on the particle weight regards the mobility status
  // Factor of 1 when mobile cell
  double alpha = 1.0;
  if (cell_mobility_status == AdaptiveSparseContacts<dim>::static_active ||
      cell_mobility_status == AdaptiveSparseContacts<dim>::advected_active)
    {
      alpha = active_status_factor;
    }
  else if (cell_mobility_status == AdaptiveSparseContacts<dim>::inactive ||
           cell_mobility_status == AdaptiveSparseContacts<dim>::advected)
    {
      alpha = inactive_status_factor;
    }

  switch (status)
    {
      case dealii::CellStatus::cell_will_persist:
      case dealii::CellStatus::cell_will_be_refined:
        {
          const unsigned int n_particles_in_cell =
            particle_handler->n_particles_in_cell(cell);
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
              particle_handler->n_particles_in_cell(cell->child(child_index));

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
