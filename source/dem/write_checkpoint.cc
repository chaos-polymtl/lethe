#include <dem/write_checkpoint.h>

using namespace dealii;

template <int dim>
void
write_checkpoint(TimerOutput &                       computing_timer,
                 const DEMSolverParameters<dim> &    parameters,
                 std::shared_ptr<SimulationControl> &simulation_control,
                 PVDHandler &                        particles_pvdhandler,
                 parallel::distributed::Triangulation<dim> &triangulation,
                 Particles::ParticleHandler<dim> &          particle_handler,
                 const ConditionalOStream &                 pcout,
                 MPI_Comm &                                 mpi_communicator)
{
  TimerOutput::Scope timer(computing_timer, "write_checkpoint");

  pcout << "Writing restart file" << std::endl;

  std::string prefix = parameters.restart.filename;
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      simulation_control->save(prefix);
      particles_pvdhandler.save(prefix);
    }

  triangulation.signals.pre_distributed_save.connect(std::bind(
    &Particles::ParticleHandler<dim>::register_store_callback_function,
    &particle_handler));

  std::ostringstream            oss;
  boost::archive::text_oarchive oa(oss, boost::archive::no_header);
  oa << particle_handler;
  std::string triangulation_name = prefix + ".triangulation";
  triangulation.save(prefix + ".triangulation");

  // Write additional particle information for deserialization
  std::string   particle_filename = prefix + ".particles";
  std::ofstream output(particle_filename.c_str());
  output << oss.str() << std::endl;
}

template void
write_checkpoint(TimerOutput &                            computing_timer,
                 const DEMSolverParameters<2> &           parameters,
                 std::shared_ptr<SimulationControl> &     simulation_control,
                 PVDHandler &                             particles_pvdhandler,
                 parallel::distributed::Triangulation<2> &triangulation,
                 Particles::ParticleHandler<2> &          particle_handler,
                 const ConditionalOStream &               pcout,
                 MPI_Comm &                               mpi_communicator);

template void
write_checkpoint(TimerOutput &                            computing_timer,
                 const DEMSolverParameters<3> &           parameters,
                 std::shared_ptr<SimulationControl> &     simulation_control,
                 PVDHandler &                             particles_pvdhandler,
                 parallel::distributed::Triangulation<3> &triangulation,
                 Particles::ParticleHandler<3> &          particle_handler,
                 const ConditionalOStream &               pcout,
                 MPI_Comm &                               mpi_communicator);
