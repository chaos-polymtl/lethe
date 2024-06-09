/**
 * @brief Displacing volume solid object (dim=3,spacedim=3).
 */

// Deal.II

#include <deal.II/base/parameter_acceptor.h>

#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <core/serial_solid.h>
#include <core/solid_objects_parameters.h>


// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim, int spacedim>
void
test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  // RigidSolidObject
  auto param = std::shared_ptr<Parameters::RigidSolidObject<spacedim>>();
  param->output_bool                   = false;
  param->solid_mesh.type               = Parameters::Mesh::Type::dealii;
  param->solid_mesh.grid_type          = "hyper_cube";
  param->solid_mesh.grid_arguments     = "-0.5 : 0.5 : false";
  param->solid_mesh.initial_refinement = 0;
  param->solid_mesh.simplex            = false;
  param->solid_mesh.translation        = Tensor<1, 3>({0.5, 0.5, 0.5});
  param->solid_mesh.rotation_axis      = Tensor<1, 3>({1., 0., 0.});
  param->solid_mesh.rotation_angle     = -0.39269908169; //  0.125 * pi
  param->center_of_rotation            = Point<3>({0., 0., 0.});

  // Functions
  param->translational_velocity        = std::make_shared<Function<dim>>();
  param->angular_velocity              = std::make_shared<Function<dim>>();


  // Output
  for (auto particle_wall_contact_list_iterator =
         particle_wall_contact_list.begin();
       particle_wall_contact_list_iterator != particle_wall_contact_list.end();
       ++particle_wall_contact_list_iterator)
    {
      auto particle_wall_candidate_content =
        &particle_wall_contact_list_iterator->second;
      for (auto particle_wall_candidate_content_iterator =
             particle_wall_candidate_content->begin();
           particle_wall_candidate_content_iterator !=
           particle_wall_candidate_content->end();
           ++particle_wall_candidate_content_iterator)
        {
          auto contact_information =
            particle_wall_candidate_content_iterator->second;
          auto particle_information = std::get<0>(contact_information);
          deallog << "Particle " << particle_information->get_id()
                  << " is located in a boundary cell" << std::endl;
        }
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test<3, 3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
