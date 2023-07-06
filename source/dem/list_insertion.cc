#include <dem/list_insertion.h>

using namespace DEM;

DeclException2(ListSizeCoherence,
               int,
               int,
               << "Incoherent particle position lists n=" << arg1
               << ", m=" << arg2);

// The constructor of list insertion class does not accomplish anything other
// than check if the position list are of the coherent size and to
// create the insertion_points member
template <int dim>
ListInsertion<dim>::ListInsertion(
  const DEMSolverParameters<dim> &dem_parameters)
  : remaining_particles_of_each_type(
      dem_parameters.lagrangian_physical_properties.number.at(0))
{
  // Inializing current inserting particle type
  current_inserting_particle_type = 0;

  const auto &list_x = dem_parameters.insertion_info.list_x;
  const auto &list_y = dem_parameters.insertion_info.list_y;
  const auto &list_z = dem_parameters.insertion_info.list_z;


  Assert(list_x.size() == list_y.size(),
         ListSizeCoherence(list_x.size(), list_y.size()));
  if (dim == 3)
    {
      Assert(list_y.size() == list_z.size(),
             ListSizeCoherence(list_y.size(), list_z.size()));
    }

  // Generate vector of insertion position
  for (unsigned int i = 0; i < list_x.size(); ++i)
    {
      if (dim == 2)
        insertion_points.emplace_back(Point<dim>({list_x[i], list_y[i]}));
      else
        insertion_points.emplace_back(
          Point<dim>({list_x[i], list_y[i], list_z[i]}));
    }
}

// The main insertion function. Insert_global_function is used to insert the
// particles
template <int dim>
void
ListInsertion<dim>::insert(
  Particles::ParticleHandler<dim> &                particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim> &                 dem_parameters)
{
  // TODO refactor into a function call
  if (remaining_particles_of_each_type == 0 &&
      current_inserting_particle_type !=
        dem_parameters.lagrangian_physical_properties.particle_type_number - 1)
    {
      remaining_particles_of_each_type =
        dem_parameters.lagrangian_physical_properties.number.at(
          ++current_inserting_particle_type);
    }

  if (remaining_particles_of_each_type > 0)
    {
      unsigned int n_total_particles_to_insert = insertion_points.size();
      n_total_particles_to_insert =
        std::min(remaining_particles_of_each_type, n_total_particles_to_insert);

      // All processors except 0 will not insert particles
      MPI_Comm communicator = triangulation.get_communicator();
      auto this_mpi_process = Utilities::MPI::this_mpi_process(communicator);
      const unsigned int n_particles_to_insert_this_proc =
        this_mpi_process == 0 ? n_total_particles_to_insert : 0;

      std::vector<Point<dim>> insertion_points_on_proc_this_step;
      insertion_points_on_proc_this_step.reserve(
        n_particles_to_insert_this_proc);

      // Because the list insertion is made to insert only a few particles
      // only processor 0 manages the insertion of these particles
      if (this_mpi_process == 0)
        {
          for (unsigned int p = 0; p < n_particles_to_insert_this_proc; ++p)
            {
              insertion_points_on_proc_this_step.emplace_back(
                insertion_points[p]);
            }
        }

      // Obtain global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);


      // Assign inserted particles properties
      this->assign_particle_properties(dem_parameters,
                                       n_particles_to_insert_this_proc,
                                       current_inserting_particle_type,
                                       this->particle_properties);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(
        insertion_points_on_proc_this_step,
        global_bounding_boxes,
        this->particle_properties);

      // Update number of particles remaining to be inserted
      remaining_particles_of_each_type -= n_total_particles_to_insert;


      ConditionalOStream pcout(
        std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);
      this->print_insertion_info(n_total_particles_to_insert,
                                 remaining_particles_of_each_type,
                                 current_inserting_particle_type,
                                 pcout);
    }
}

template class ListInsertion<2>;
template class ListInsertion<3>;
