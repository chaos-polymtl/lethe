//
// Created by lucka on 2021-10-11.
//
#include <fem-dem/ib_particles_dem.h>


template <int dim>
void
IBParticlesDEM<dim>::initialize(
  SimulationParameters<dim> p_nsparam, MPI_Comm&     mpi_communicator_input,std::vector<IBParticle<dim>> particles){
  parameters = p_nsparam;
  mpi_communicator= mpi_communicator_input;
  dem_particles=particles;
  boundary_cells.resize(dem_particles.size());
}
template <int dim>
void
IBParticlesDEM<dim>::update_particles(
  std::vector<IBParticle<dim>> particles,double time){
  dem_particles=particles;
  cfd_time=time;

}


template <int dim>
void
IBParticlesDEM<dim>::calculate_pp_contact_force(
  const double                &dt_dem,
  std::vector<Tensor<1, dim>> &contact_force,
  std::vector<Tensor<1, 3>>   &contact_torque)
{
  for (auto &particle_one : dem_particles)
    {
      for (auto &particle_two : dem_particles)
        {
          if (particle_one.particle_id != particle_two.particle_id and
              particle_one.particle_id < particle_two.particle_id)
            {
              const Point<dim> particle_one_location = particle_one.position;
              const Point<dim> particle_two_location = particle_two.position;
              ContactTangentialHistory contact_history;
              ContactTangentialHistory contact_info;
              try
                {
                  contact_info = pp_contact_map[particle_one.particle_id]
                                               [particle_two.particle_id];
                }
              catch (...)
                {
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_history.tangential_overlap[d]           = 0;
                      contact_history.tangential_relative_velocity[d] = 0;
                    }
                  pp_contact_map[particle_one.particle_id]
                                [particle_two.particle_id] = contact_history;
                }

              double effective_mass = (particle_one.mass * particle_two.mass) /
                                      (particle_one.mass + particle_two.mass);
              double effective_radius =
                (particle_one.radius * particle_two.radius) /
                (particle_one.radius + particle_two.radius);


              // Calculation of normal overlap
              double normal_overlap =
                (particle_one.radius + particle_two.radius) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0)
                // This means that the adjacent particles are in contact
                {
                  // Calculation of the contact vector (vector from particle one
                  // to particle two
                  auto contact_vector =
                    particle_two_location - particle_one_location;

                  // Using contact_vector, the contact normal vector is obtained
                  auto normal_unit_vector =
                    contact_vector / contact_vector.norm();
                  Tensor<1, 3> normal;
                  if (dim == 2)
                    {
                      normal[0] = normal_unit_vector[0];
                      normal[1] = normal_unit_vector[1];
                    }
                  if (dim == 3)
                    {
                      normal[0] = normal_unit_vector[0];
                      normal[1] = normal_unit_vector[1];
                      normal[2] = normal_unit_vector[2];
                    }

                  // Defining relative contact velocity
                  Tensor<1, dim> contact_relative_velocity;

                  if (dim == 3)
                    {
                      // Calculation of contact relative velocity
                      Tensor<1, 3> rotational_velocity = cross_product_3d(
                        (particle_one.radius * particle_one.omega +
                         particle_two.radius * particle_two.omega),
                        normal);

                      contact_relative_velocity[0] =
                        (particle_one.velocity - particle_two.velocity)[0] +
                        rotational_velocity[0];
                      contact_relative_velocity[1] =
                        (particle_one.velocity - particle_two.velocity)[1] +
                        rotational_velocity[1];
                      contact_relative_velocity[2] =
                        (particle_one.velocity - particle_two.velocity)[2] +
                        rotational_velocity[2];
                    }
                  else
                    {
                      // TODO: Correct this after correcting dem_2d
                      contact_relative_velocity =
                        particle_one.velocity - particle_two.velocity;
                    }

                  auto normal_relative_velocity_value =
                    contact_relative_velocity * normal_unit_vector;
                  Tensor<1, dim> normal_relative_velocity =
                    normal_relative_velocity_value * normal_unit_vector;

                  // Calculation of tangential relative velocity
                  Tensor<1, dim> tangential_relative_velocity =
                    contact_relative_velocity - normal_relative_velocity;

                  Tensor<1, dim> modified_tangential_overlap =
                    contact_info.tangential_overlap +
                    contact_info.tangential_relative_velocity * dt_dem;

                  // Updating the contact_info container based on the new
                  // calculated values
                  contact_history.tangential_overlap =
                    modified_tangential_overlap;
                  contact_history.tangential_relative_velocity =
                    tangential_relative_velocity;
                  pp_contact_map[particle_one.particle_id]
                                [particle_two.particle_id] = contact_history;


                  const double effective_mass =
                    (particle_one.mass * particle_two.mass) /
                    (particle_one.mass + particle_two.mass);
                  const double effective_radius =
                    (particle_one.radius * particle_two.radius) /
                    (particle_one.radius + particle_two.radius);
                  const double effective_youngs_modulus =
                    (particle_one.youngs_modulus *
                     particle_two.youngs_modulus) /
                    ((particle_two.youngs_modulus *
                      (1 - particle_one.poisson_ratio *
                             particle_one.poisson_ratio)) +
                     (particle_one.youngs_modulus *
                      (1 - particle_two.poisson_ratio *
                             particle_two.poisson_ratio)) +
                     DBL_MIN);

                  const double effective_shear_modulus =
                    (particle_one.youngs_modulus *
                     particle_two.youngs_modulus) /
                    (2 * ((particle_two.youngs_modulus *
                           (2 - particle_one.poisson_ratio) *
                           (1 + particle_one.poisson_ratio)) +
                          (particle_one.youngs_modulus *
                           (2 - particle_two.poisson_ratio) *
                           (1 + particle_two.poisson_ratio))) +
                     DBL_MIN);


                  const double effective_coefficient_of_restitution =
                    2 * particle_one.restitution_coefficient *
                    particle_two.restitution_coefficient /
                    (particle_one.restitution_coefficient +
                     particle_two.restitution_coefficient + DBL_MIN);

                  const double effective_coefficient_of_friction =
                    2 * particle_one.friction_coefficient *
                    particle_two.friction_coefficient /
                    (particle_one.friction_coefficient +
                     particle_two.friction_coefficient + DBL_MIN);


                  const double effective_coefficient_of_rolling_friction =
                    2 * particle_one.rolling_friction_coefficient *
                    particle_two.rolling_friction_coefficient /
                    (particle_one.rolling_friction_coefficient +
                     particle_two.rolling_friction_coefficient + DBL_MIN);

                  const double restitution_coefficient_particle_log =
                    std::log(effective_coefficient_of_restitution);

                  const double model_parameter_beta =
                    restitution_coefficient_particle_log /
                    sqrt(restitution_coefficient_particle_log *
                           restitution_coefficient_particle_log +
                         9.8696);

                  const double radius_times_overlap_sqrt =
                    sqrt(effective_radius * normal_overlap);

                  const double model_parameter_sn =
                    2 * effective_youngs_modulus * radius_times_overlap_sqrt;

                  const double model_parameter_st =
                    8 * effective_youngs_modulus * radius_times_overlap_sqrt;

                  // Calculation of normal and tangential spring and dashpot
                  // constants using particle properties
                  const double normal_spring_constant =
                    0.66665 * model_parameter_sn;
                  const double normal_damping_constant =
                    -1.8257 * model_parameter_beta *
                    sqrt(model_parameter_sn * effective_mass);
                  const double tangential_spring_constant =
                    8 * effective_shear_modulus * radius_times_overlap_sqrt +
                    DBL_MIN;
                  const double tangential_damping_constant =
                    normal_damping_constant *
                    sqrt(model_parameter_st / model_parameter_sn);

                  // Calculation of normal force using spring and dashpot normal
                  // forces
                  Tensor<1, dim> normal_force =
                    ((normal_spring_constant * normal_overlap) *
                     normal_unit_vector) +
                    ((normal_damping_constant *
                      normal_relative_velocity_value) *
                     normal_unit_vector);


                  // Calculation of tangential force using spring and dashpot
                  // tangential forces. Since we need dashpot tangential force
                  // in the gross sliding again, we define it as a separate
                  // variable
                  Tensor<1, dim> damping_tangential_force =
                    tangential_damping_constant *
                    contact_info.tangential_relative_velocity;
                  Tensor<1, dim> tangential_force =
                    (tangential_spring_constant *
                     contact_info.tangential_overlap) +
                    damping_tangential_force;

                  const double coulomb_threshold =
                    effective_coefficient_of_friction * normal_force.norm();

                  // Check for gross sliding
                  if (tangential_force.norm() > coulomb_threshold)
                    {
                      // Gross sliding occurs and the tangential overlap and
                      // tangnetial force are limited to Coulumb's criterion
                      contact_info.tangential_overlap =
                        (coulomb_threshold *
                           (tangential_force /
                            (tangential_force.norm() + DBL_MIN)) -
                         damping_tangential_force) /
                        (tangential_spring_constant + DBL_MIN);

                      tangential_force = (tangential_spring_constant *
                                          contact_info.tangential_overlap) +
                                         damping_tangential_force;
                    }

                  // Calculation of torque
                  // Torque caused by tangential force (tangential_torque)
                  Tensor<1, 3> tangential_torque;
                  if (dim == 2)
                    {
                      Tensor<1, 3> tangential_force_3d;
                      tangential_force_3d[0] = tangential_force[0];
                      tangential_force_3d[1] = tangential_force[1];
                      tangential_torque =
                        cross_product_3d((effective_radius * normal),
                                         tangential_force_3d);
                    }
                  if (dim == 3)
                    {
                      Tensor<1, 3> tangential_force_3d;
                      tangential_force_3d[0] = tangential_force[0];
                      tangential_force_3d[1] = tangential_force[1];
                      tangential_force_3d[2] = tangential_force[2];
                      tangential_torque =
                        cross_product_3d((effective_radius * normal),
                                         tangential_force_3d);
                    }
                  // Rolling resistance torque using viscous rolling resistance
                  // model
                  auto omega_ij = particle_one.omega - particle_two.omega;
                  auto omega_ij_direction =
                    omega_ij / (omega_ij.norm() + DBL_MIN);

                  Tensor<1, 3> v_omega =
                    cross_product_3d(particle_one.omega,
                                     particle_one.radius * normal) -
                    cross_product_3d(particle_two.omega,
                                     particle_two.radius * -normal);

                  // Calculation of rolling resistance torque
                  Tensor<1, 3> rolling_resistance_torque =
                    -effective_coefficient_of_rolling_friction *
                    effective_radius * normal_force.norm() * v_omega.norm() *
                    omega_ij_direction;



                  contact_force[particle_one.particle_id] -=
                    (normal_force + tangential_force);
                  contact_force[particle_two.particle_id] +=
                    (normal_force + tangential_force);

                  contact_torque[particle_one.particle_id] =
                    contact_torque[particle_one.particle_id] -
                    tangential_torque + rolling_resistance_torque;
                  contact_torque[particle_two.particle_id] =
                    contact_torque[particle_two.particle_id] -
                    tangential_torque - rolling_resistance_torque;
                }

              else
                {
                  // if the adjacent pair is not in contact anymore
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_history.tangential_overlap[d]           = 0;
                      contact_history.tangential_relative_velocity[d] = 0;
                    }
                  pp_contact_map[particle_one.particle_id].erase(
                    particle_two.particle_id);
                }
            }
        }
    }
}
template <int dim>
void
IBParticlesDEM<dim>::update_particles_boundary_contact(std::vector<IBParticle<dim>>& particles,DoFHandler<dim>& dof_handler)
{

  for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
    {
      boundary_cells[p_i].clear();
      auto cells_at_boundary =
        LetheGridTools::find_boundary_cell_in_sphere(dof_handler,
                                                     particles[p_i].position,
                                                     particles[p_i].radius *
                                                       1.5);
      // Initialize a simple quadrature for on the system. This will be used to
      // obtain a single sample point on the boundary faces
      const FE_Q<dim> fe(1);
      for (unsigned int i = 0; i < cells_at_boundary.size(); ++i)
        {
          QGauss<dim - 1>   face_quadrature_formula(1);
          unsigned int      n_face_q_points = face_quadrature_formula.size();
          FEFaceValues<dim> fe_face_values(fe,
                                           face_quadrature_formula,
                                           update_values |
                                             update_quadrature_points |
                                             update_normal_vectors);
          for (int face_id = 0;
               face_id < int(GeometryInfo<dim>::faces_per_cell);
               ++face_id)
            {
              if (cells_at_boundary[i]->face(face_id)->at_boundary())
                {
                  fe_face_values.reinit(cells_at_boundary[i], face_id);

                  for (unsigned int f_q_point = 0; f_q_point < n_face_q_points;
                       ++f_q_point)
                    {
                      Tensor<1, dim> normal_vector =
                        -fe_face_values.normal_vector(f_q_point);
                      BoundaryCellsInfo boundary_information;
                      boundary_information.normal_vector = normal_vector;
                      boundary_information.point_on_boundary =
                        fe_face_values.quadrature_point(f_q_point);
                      boundary_cells[p_i].push_back(boundary_information);
                    }
                }
            }
        }


      auto
        global_boundary_cell=Utilities::MPI::all_gather(this->mpi_communicator,boundary_cells[p_i]);
      boundary_cells[p_i].clear();
      for(unsigned int i=0;i<global_boundary_cell.size();++i){
          boundary_cells[p_i].insert(boundary_cells[p_i].end(),
                                     global_boundary_cell[i].begin(), global_boundary_cell[i].end());
        }

    }
}


template <int dim>
void
IBParticlesDEM<dim>::calculate_pw_contact_force(
  const double                &dt_dem,
  std::vector<Tensor<1, dim>> &contact_force,
  std::vector<Tensor<1, 3>>   &contact_torque)
{
  double wall_youngs_modulus =
    parameters.particlesParameters->wall_youngs_modulus;
  double wall_poisson_ratio =
    parameters.particlesParameters->wall_poisson_ratio;
  double wall_rolling_friction_coefficient =
    parameters.particlesParameters
      ->wall_rolling_friction_coefficient;
  double wall_friction_coefficient =
    parameters.particlesParameters->wall_friction_coefficient;
  double wall_restitution_coefficient =
    parameters.particlesParameters
      ->wall_restitution_coefficient;


  for (auto &particle : dem_particles)
    {
      unsigned int boundary_index = 0;
      double       best_dist      = DBL_MAX;
      unsigned int best_index;

      for (auto &boundary_cell_iter : boundary_cells[particle.particle_id])
        {
          double dist =
            (boundary_cell_iter.point_on_boundary - particle.position).norm();
          if (dist < best_dist)
            {
              best_dist  = dist;
              best_index = boundary_index;
            }
          boundary_index += 1;
        }
      if (boundary_cells[particle.particle_id].size() > 0)
        {
          auto &boundary_cell =
            boundary_cells[particle.particle_id][best_index];


          auto                     boundary_cell_information = boundary_cell;
          ContactTangentialHistory contact_history;
          ContactTangentialHistory contact_info;
          try
            {
              contact_info =
                pw_contact_map[particle.particle_id][boundary_index];
            }
          catch (...)
            {
              for (int d = 0; d < dim; ++d)
                {
                  contact_history.tangential_overlap[d]           = 0;
                  contact_history.tangential_relative_velocity[d] = 0;
                }
              pw_contact_map[particle.particle_id][boundary_index] =
                contact_history;
            }

          auto normal_vector     = boundary_cell_information.normal_vector;
          auto point_on_boundary = boundary_cell_information.point_on_boundary;

          Tensor<1, 3> normal;
          if (dim == 2)
            {
              normal[0] = normal_vector[0];
              normal[1] = normal_vector[1];
            }
          if (dim == 3)
            {
              normal[0] = normal_vector[0];
              normal[1] = normal_vector[1];
              normal[2] = normal_vector[2];
            }

          // A vector (point_to_particle_vector) is defined which connects the
          // center of particle to the point_on_boundary. This vector will then
          // be projected on the normal vector of the boundary to obtain the
          // particle-wall distance
          Tensor<1, dim> point_to_particle_vector =
            particle.position - point_on_boundary;

          // Finding the projected vector on the normal vector of the boundary.
          // Here we have used the private function find_projection. Using this
          // projected vector, the particle-wall distance is calculated
          Tensor<1, dim> projected_vector =
            find_projection(point_to_particle_vector, normal_vector);

          double normal_overlap = particle.radius - (projected_vector.norm());

          if (normal_overlap > 0)
            {
              // Defining relative contact velocity
              Tensor<1, dim> contact_relative_velocity;
              if (dim == 3)
                {
                  Tensor<1, 3> rotational_velocity =
                    cross_product_3d((particle.radius * particle.omega),
                                     normal);
                  contact_relative_velocity[0] =
                    particle.velocity[0] + rotational_velocity[0];
                  contact_relative_velocity[1] =
                    particle.velocity[1] + rotational_velocity[1];
                  contact_relative_velocity[2] =
                    particle.velocity[2] + rotational_velocity[2];
                }
              if (dim == 2)
                {
                  // TODO correct this after dem_2d correction
                  contact_relative_velocity = particle.velocity;
                }

              // Calculation of normal relative velocity
              double normal_relative_velocity_value =
                contact_relative_velocity * normal_vector;
              Tensor<1, dim> normal_relative_velocity =
                normal_relative_velocity_value * normal_vector;

              // Calculation of tangential relative velocity
              Tensor<1, dim> tangential_relative_velocity =
                contact_relative_velocity - normal_relative_velocity;

              Tensor<1, dim> modified_tangential_overlap =
                contact_info.tangential_overlap +
                tangential_relative_velocity * dt_dem;


              const double effective_youngs_modulus =
                (particle.youngs_modulus * wall_youngs_modulus) /
                (wall_youngs_modulus *
                   (1 - particle.poisson_ratio * particle.poisson_ratio) +
                 particle.youngs_modulus *
                   (1 - wall_poisson_ratio * wall_poisson_ratio) +
                 DBL_MIN);

              const double effective_shear_modulus =
                (particle.youngs_modulus * wall_youngs_modulus) /
                ((2 * wall_youngs_modulus * (2 - particle.poisson_ratio) *
                  (1 + particle.poisson_ratio)) +
                 (2 * particle.youngs_modulus * (2 - wall_poisson_ratio) *
                  (1 + wall_poisson_ratio)) +
                 DBL_MIN);

              const double effective_coefficient_of_restitution =
                2 * particle.restitution_coefficient *
                wall_restitution_coefficient /
                (particle.restitution_coefficient +
                 wall_restitution_coefficient + DBL_MIN);

              const double effective_coefficient_of_friction =
                2 * particle.friction_coefficient * wall_friction_coefficient /
                (particle.friction_coefficient + wall_friction_coefficient +
                 DBL_MIN);

              const double effective_coefficient_of_rolling_friction =
                2 * particle.rolling_friction_coefficient *
                wall_rolling_friction_coefficient /
                (particle.rolling_friction_coefficient +
                 wall_rolling_friction_coefficient + DBL_MIN);

              const double radius_times_overlap_sqrt =
                sqrt(particle.radius * normal_overlap);
              const double log_coeff_restitution =
                log(effective_coefficient_of_restitution);
              const double model_parameter_beta =
                log_coeff_restitution /
                sqrt((log_coeff_restitution * log_coeff_restitution) + 9.8696);
              const double model_parameter_sn =
                2 * effective_youngs_modulus * radius_times_overlap_sqrt;

              // Calculation of normal and tangential spring and dashpot
              // constants using particle and wall properties
              const double normal_spring_constant =
                1.3333 * effective_youngs_modulus * radius_times_overlap_sqrt;
              const double normal_damping_constant =
                1.8257 * model_parameter_beta *
                sqrt(model_parameter_sn * particle.mass);
              const double tangential_spring_constant =
                -8 * effective_shear_modulus * radius_times_overlap_sqrt +
                DBL_MIN;

              // Calculation of normal force using spring and dashpot normal
              // forces

              Tensor<1, dim> normal_force =
                (normal_spring_constant * normal_overlap +
                 normal_damping_constant * normal_relative_velocity_value) *
                normal_vector;

              // Calculation of tangential force
              Tensor<1, dim> tangential_force =
                tangential_spring_constant * contact_info.tangential_overlap;

              const double coulomb_threshold =
                effective_coefficient_of_friction * normal_force.norm();

              // Check for gross sliding
              if (tangential_force.norm() > coulomb_threshold)
                {
                  // Gross sliding occurs and the tangential overlap and
                  // tangnetial force are limited to Coulumb's criterion
                  tangential_force =
                    coulomb_threshold *
                    (tangential_force / tangential_force.norm());

                  contact_info.tangential_overlap =
                    tangential_force / (tangential_spring_constant + DBL_MIN);
                }

              // Calculation of torque
              // Torque caused by tangential force (tangential_torque)
              Tensor<1, 3> tangential_torque;
              if (dim == 2)
                {
                  Tensor<1, 3> tangential_force_3d;
                  tangential_force_3d[0] = tangential_force[0];
                  tangential_force_3d[1] = tangential_force[1];
                  tangential_torque =
                    cross_product_3d((particle.radius * normal),
                                     tangential_force_3d);
                }
              if (dim == 3)
                {
                  Tensor<1, 3> tangential_force_3d;
                  tangential_force_3d[0] = tangential_force[0];
                  tangential_force_3d[1] = tangential_force[1];
                  tangential_force_3d[2] = tangential_force[2];
                  tangential_torque =
                    cross_product_3d((particle.radius * normal),
                                     tangential_force_3d);
                }

              // Rolling resistance torque
              Tensor<1, 3> angular_velocity = particle.omega;
              Tensor<1, 3> pw_angular_velocity;

              double omega_value = angular_velocity.norm();
              if (omega_value != 0)
                {
                  pw_angular_velocity = angular_velocity / omega_value;
                }

              Tensor<1, 3> v_omega =
                cross_product_3d(angular_velocity, particle.radius * normal);

              // Calculation of rolling resistance torque
              Tensor<1, 3> rolling_resistance_torque =
                -effective_coefficient_of_rolling_friction * particle.radius *
                normal_force.norm() * v_omega.norm() * pw_angular_velocity;


              // Updating the force of particles in the particle handler
              contact_force[particle.particle_id] +=
                normal_force + tangential_force;

              // Updating the torque acting on particles

              contact_torque[particle.particle_id] +=
                tangential_torque + rolling_resistance_torque;
            }
          else
            {
              for (int d = 0; d < dim; ++d)
                {
                  contact_info.tangential_overlap[d] = 0;
                }
            }
        }
    }
}
template <int dim>
void
IBParticlesDEM<dim>::particles_dem(double dt)
{ // add refilling containers
  using numbers::PI;


  double         rho =parameters.particlesParameters->density;
  double         dt_dem             = dt /parameters.particlesParameters->coupling_frequency ;

  std::vector<Tensor<1, dim>> contact_force(dem_particles.size());
  std::vector<Tensor<1, 3>>   contact_torque(dem_particles.size());
  std::vector<Tensor<1, 3>>   contact_wall_torque(dem_particles.size());
  std::vector<Tensor<1, dim>> contact_wall_force(dem_particles.size());

  std::vector<Tensor<1, dim>> current_fluid_force(dem_particles.size());
  std::vector<Tensor<1, 3>>   current_fluid_torque(dem_particles.size());

  std::vector<Tensor<1, dim>> velocity(dem_particles.size());
  std::vector<Point<dim>>     position(dem_particles.size());

  // local time for the dem step
  double         t = 0;
  Tensor<1, dim> g;
  this->parameters.particlesParameters->f_gravity->set_time(cfd_time);
  Tensor<1, dim> gravity;
  // initialized the particles
  for (unsigned int p_i = 0; p_i < dem_particles.size(); ++p_i)
    {
      dem_particles[p_i].position    = dem_particles[p_i].last_position[0];
      dem_particles[p_i].velocity    = dem_particles[p_i].last_velocity[0];
      dem_particles[p_i].omega       = dem_particles[p_i].last_omega[0];
      dem_particles[p_i].impulsion       = 0;
      dem_particles[p_i].omega_impulsion = 0;
      dem_particles[p_i].contact_impulsion = 0;
      g[0] = this->parameters.particlesParameters->f_gravity
               ->value(dem_particles[p_i].position, 0);
      g[1] = this->parameters.particlesParameters->f_gravity
               ->value(dem_particles[p_i].position, 1);
      if(dim==3)
        g[2] = this->parameters.particlesParameters->f_gravity
               ->value(dem_particles[p_i].position, 2);
    }

  // integrate on the with the sub_time_step
  while (t+dt_dem/2 < dt)
    {
      current_fluid_force.clear();
      current_fluid_force.resize(dem_particles.size());
      current_fluid_torque.clear();
      current_fluid_torque.resize(dem_particles.size());

      contact_torque.clear();
      contact_torque.resize(dem_particles.size());
      contact_force.clear();
      contact_force.resize(dem_particles.size());
      contact_wall_force.clear();
      contact_wall_force.resize(dem_particles.size());
      contact_wall_torque.clear();
      contact_wall_torque.resize(dem_particles.size());

      // Calculate particle-particle and particle-wall contact force
      calculate_pp_contact_force(dt_dem, contact_force, contact_torque);
      calculate_pw_contact_force(dt_dem,
                                 contact_wall_force,
                                 contact_wall_torque);


      for (unsigned int p_i = 0; p_i < dem_particles.size(); ++p_i)
        {
          auto inv_inertia = invert(dem_particles[p_i].inertia);
          if (dim == 2)
            gravity =
              g * (dem_particles[p_i].mass -
                   dem_particles[p_i].radius * dem_particles[p_i].radius * PI * rho);
          if (dim == 3)
            {
              gravity =
                g * (dem_particles[p_i].mass - 4.0 / 3.0 * dem_particles[p_i].radius *
                                             dem_particles[p_i].radius *
                                             dem_particles[p_i].radius * PI * rho);
            }

          // BDF 1 on the force
          current_fluid_force[p_i] =dem_particles[p_i].forces;
          current_fluid_torque[p_i] =dem_particles[p_i].torques;

          // Explicite Euler
          dem_particles[p_i].velocity =
            dem_particles[p_i].velocity +
            (current_fluid_force[p_i] + contact_force[p_i] +
             contact_wall_force[p_i] + gravity) *
              dt_dem / dem_particles[p_i].mass;
          dem_particles[p_i].position =
            dem_particles[p_i].position + dem_particles[p_i].velocity * dt_dem;

          dem_particles[p_i].omega =
            dem_particles[p_i].omega +
            inv_inertia *
              (current_fluid_torque[p_i] + contact_torque[p_i] +
               contact_wall_torque[p_i]) *
              dt_dem;

          // Integration of the impulsion
          dem_particles[p_i].impulsion +=
            (current_fluid_force[p_i] + gravity) *
            dt_dem;
          dem_particles[p_i].contact_impulsion+=(contact_wall_force[p_i]+ contact_force[p_i])*
                                              dt_dem;
          dem_particles[p_i].omega_impulsion +=
            (current_fluid_torque[p_i] + contact_torque[p_i] +
             contact_wall_torque[p_i]) *
            dt_dem;
        }
      t += dt_dem;
    }
}


template class IBParticlesDEM<2>;
template class IBParticlesDEM<3>;