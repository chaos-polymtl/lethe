#include <dem/read_checkpoint.h>

using namespace dealii;

template <int dim>
void
read_checkpoint(
  TimerOutput                                             &computing_timer,
  const DEMSolverParameters<dim>                          &parameters,
  std::shared_ptr<SimulationControl>                      &simulation_control,
  PVDHandler                                              &particles_pvdhandler,
  PVDHandler                                              &grid_pvdhandler,
  parallel::distributed::Triangulation<dim>               &triangulation,
  Particles::ParticleHandler<dim>                         &particle_handler,
  std::shared_ptr<Insertion<dim>>                         &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solid_surfaces)
{
  TimerOutput::Scope timer(computing_timer, "read_checkpoint");
  std::string        prefix = parameters.restart.filename;
  simulation_control->read(prefix);
  particles_pvdhandler.read(prefix);

  if (parameters.post_processing.Lagrangian_post_processing)
    {
      grid_pvdhandler.read(prefix + "_postprocess_data");
    }

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


  // Unpack the information in the particle handler
  particle_handler.deserialize();


  // Load insertion object
  std::string   insertion_object_filename = prefix + ".insertion_object";
  std::ifstream iss_insertion_obj(insertion_object_filename);
  boost::archive::text_iarchive ia_insertion_obj(iss_insertion_obj,
                                                 boost::archive::no_header);
  insertion_object->deserialize(ia_insertion_obj, 0);

  // Load solid surfaces
  for (unsigned int i = 0; i < solid_surfaces.size(); ++i)
    {
      solid_surfaces[i]->read_checkpoint(prefix);
    }
}

template void
read_checkpoint(
  TimerOutput                                     &computing_timer,
  const DEMSolverParameters<2>                    &parameters,
  std::shared_ptr<SimulationControl>              &simulation_control,
  PVDHandler                                      &particles_pvdhandler,
  PVDHandler                                      &grid_pvdhandler,
  parallel::distributed::Triangulation<2>         &triangulation,
  Particles::ParticleHandler<2>                   &particle_handler,
  std::shared_ptr<Insertion<2>>                   &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<1, 2>>> &solid_surfaces);

template void
read_checkpoint(
  TimerOutput                                     &computing_timer,
  const DEMSolverParameters<3>                    &parameters,
  std::shared_ptr<SimulationControl>              &simulation_control,
  PVDHandler                                      &particles_pvdhandler,
  PVDHandler                                      &grid_pvdhandler,
  parallel::distributed::Triangulation<3>         &triangulation,
  Particles::ParticleHandler<3>                   &particle_handler,
  std::shared_ptr<Insertion<3>>                   &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<2, 3>>> &solid_surfaces);
