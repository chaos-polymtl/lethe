#include <dem/write_checkpoint.h>

using namespace dealii;

template <int dim>
void
write_checkpoint(
  TimerOutput                                             &computing_timer,
  const DEMSolverParameters<dim>                          &parameters,
  std::shared_ptr<SimulationControl>                      &simulation_control,
  PVDHandler                                              &particles_pvdhandler,
  PVDHandler                                              &grid_pvdhandler,
  parallel::distributed::Triangulation<dim>               &triangulation,
  Particles::ParticleHandler<dim>                         &particle_handler,
  std::shared_ptr<Insertion<dim>>                         &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solid_objects,
  const ConditionalOStream                                &pcout,
  MPI_Comm                                                &mpi_communicator)
{
  TimerOutput::Scope timer(computing_timer, "write_checkpoint");

  pcout << "Writing restart file" << std::endl;

  std::string prefix = parameters.restart.filename;
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      simulation_control->save(prefix);
      particles_pvdhandler.save(prefix);

      if (parameters.post_processing.Lagrangian_post_processing)
        {
          grid_pvdhandler.save(prefix + "_postprocess_data");
        }
    }

  // Prepare the particle handler for checkpointing
  particle_handler.prepare_for_serialization();

  std::ostringstream            oss;
  boost::archive::text_oarchive oa(oss, boost::archive::no_header);
  oa << particle_handler;
  std::string triangulation_name = prefix + ".triangulation";
  triangulation.save(prefix + ".triangulation");

  // Write additional particle information for deserialization
  std::string   particle_filename = prefix + ".particles";
  std::ofstream output(particle_filename.c_str());
  output << oss.str() << std::endl;

  // Prepare the insertion object for checkpointing
  // insertion_object.prepare_for_serialization();

  std::string   insertion_object_filename = prefix + ".insertion_object";
  std::ofstream oss_insertion_obj(insertion_object_filename);
  boost::archive::text_oarchive oa_insertion_obj(oss_insertion_obj,
                                                 boost::archive::no_header);
  insertion_object->serialize(oa_insertion_obj, 0);

  // Checkpoint the serial solid objects one by one
  for (unsigned int i = 0; i < solid_objects.size(); ++i)
    {
      solid_objects[i]->write_checkpoint(prefix);
    }
}

template void
write_checkpoint(TimerOutput                             &computing_timer,
                 const DEMSolverParameters<2>            &parameters,
                 std::shared_ptr<SimulationControl>      &simulation_control,
                 PVDHandler                              &particles_pvdhandler,
                 PVDHandler                              &grid_pvdhandler,
                 parallel::distributed::Triangulation<2> &triangulation,
                 Particles::ParticleHandler<2>           &particle_handler,
                 std::shared_ptr<Insertion<2>>           &insertion_object,
                 std::vector<std::shared_ptr<SerialSolid<1, 2>>> &solid_objects,
                 const ConditionalOStream                        &pcout,
                 MPI_Comm &mpi_communicator);

template void
write_checkpoint(TimerOutput                             &computing_timer,
                 const DEMSolverParameters<3>            &parameters,
                 std::shared_ptr<SimulationControl>      &simulation_control,
                 PVDHandler                              &particles_pvdhandler,
                 PVDHandler                              &grid_pvdhandler,
                 parallel::distributed::Triangulation<3> &triangulation,
                 Particles::ParticleHandler<3>           &particle_handler,
                 std::shared_ptr<Insertion<3>>           &insertion_object,
                 std::vector<std::shared_ptr<SerialSolid<2, 3>>> &solid_objects,
                 const ConditionalOStream                        &pcout,
                 MPI_Comm &mpi_communicator);
