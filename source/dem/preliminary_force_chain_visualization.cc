
#include <core/manifolds.h>
#include <core/solutions_output.h>

#include <dem/data_containers.h>
#include <dem/dem.h>
#include <dem/distributions.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/find_contact_detection_step.h>
#include <dem/gear3_integrator.h>
#include <dem/input_parameter_inspection.h>
#include <dem/insertion_file.h>
#include <dem/insertion_list.h>
#include <dem/insertion_plane.h>
#include <dem/insertion_volume.h>
#include <dem/particle_wall_nonlinear_force.h>
#include <dem/post_processing.h>
#include <dem/read_checkpoint.h>
#include <dem/read_mesh.h>
#include <dem/set_particle_particle_contact_force_model.h>
#include <dem/set_particle_wall_contact_force_model.h>
#include <dem/velocity_verlet_integrator.h>
#include <dem/write_checkpoint.h>
#include <dem/visualization.h>
#include <deal.II/base/table_handler.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_out.h>

#include <sys/stat.h>

#include <sstream>



template <int dim>
void
DEMForceChain<dim>::write_output_results()
{
  TimerOutput::Scope t(this->computing_timer, "Output VTU");

  const std::string folder = "./output/";
  const std::string particles_solution_name =
    "Force_chain_test";
  const unsigned int iter        = 1;
  const double       time        = 0;
  const unsigned int group_files = 1;

  // Write particles
  Visualization<dim> particle_data_out;
  particle_data_out.build_patches(particle_handler,
                                  properties_class.get_properties_name());

  write_vtu_and_pvd<0, dim>(particles_pvdhandler,
                            particle_data_out,
                            folder,
                            particles_solution_name,
                            time,
                            iter,
                            group_files,
                            mpi_communicator);

  if (simulation_control->get_output_boundaries())
    {
      DataOutFaces<dim> data_out_faces;

      // Setup background dofs to initiate right sized boundary_id vector
      Vector<float> boundary_id(background_dh.n_dofs());

      // Attach the boundary_id to data_out_faces object
      BoundaryPostprocessor<dim> boundary;
      data_out_faces.attach_dof_handler(background_dh);
      data_out_faces.add_data_vector(boundary_id, boundary);
      data_out_faces.build_patches();

      write_boundaries_vtu<dim>(
        data_out_faces, folder, time, iter, this->mpi_communicator);
    }

  // Write all solid objects
  for (const auto &solid_object : solids)
    solid_object->write_output_results(simulation_control);
}


template <int dim>
void
DEMForceChain<dim>::solve()
{
  
    write_output_results();
        
}

template class DEMForceChain<2>;
template class DEMForceChain<3>;