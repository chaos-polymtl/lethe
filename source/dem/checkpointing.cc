#include <dem/checkpointing.h>

using namespace dealii;

template <int dim>
Checkpointing<dim>::Checkpointing()
{}

template <int dim>
void
Checkpointing<dim>::read_checkpoint(
  TimerOutput &                              computing_timer,
  const DEMSolverParameters<dim> &           parameters,
  std::shared_ptr<SimulationControl> &       simulation_control,
  PVDHandler &                               particles_pvdhandler,
  parallel::distributed::Triangulation<dim> &triangulation,
  Particles::ParticleHandler<dim> &          particle_handler)
{
  TimerOutput::Scope timer(computing_timer, "read_checkpoint");
  std::string        prefix = parameters.restart.filename;
  simulation_control->read(prefix);
  particles_pvdhandler.read(prefix);

  triangulation.signals.post_distributed_load.connect(
    std::bind(&Particles::ParticleHandler<dim>::register_load_callback_function,
              &particle_handler,
              true));

  // Gather particle serialization information
  std::string   particle_filename = prefix + ".particles";
  std::ifstream input(particle_filename.c_str());
  AssertThrow(input, ExcFileNotOpen(particle_filename));

  std::string buffer;
  std::getline(input, buffer);
  std::istringstream            iss(buffer);
  boost::archive::text_iarchive ia(iss, boost::archive::no_header);

  ia >> particle_handler;

  const std::string filename = prefix + ".triangulation";
  std::ifstream     in(filename.c_str());
  if (!in)
    AssertThrow(false,
                ExcMessage(
                  std::string(
                    "You are trying to restart a previous computation, "
                    "but the restart file <") +
                  filename + "> does not appear to exist!"));

  try
    {
      triangulation.load(filename.c_str());
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("Cannot open snapshot mesh file or read the "
                             "triangulation stored there."));
    }
}

template <int dim>
void
Checkpointing<dim>::write_checkpoint(
  TimerOutput &                              computing_timer,
  const DEMSolverParameters<dim> &           parameters,
  std::shared_ptr<SimulationControl> &       simulation_control,
  PVDHandler &                               particles_pvdhandler,
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

template class Checkpointing<2>;
template class Checkpointing<3>;
