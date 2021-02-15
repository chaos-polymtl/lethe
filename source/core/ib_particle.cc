#include <core/ib_particle.h>


template <int dim>
void
IBParticle<dim>::initialise_all() {
    mass=1;
    radius=1;
    local_alpha_torque=1;
    local_alpha_force=1;

    inertia[0][0]=0;
    inertia[1][1]=0;
    inertia[2][2]=0;

    forces[0] = 0;
    forces[1] = 0;

    velocity[0] = 0;
    velocity[1] = 0;

    position[0] = 0;
    position[1] = 0;

    torques[0] = 0;
    torques[1] = 0;
    torques[2] = 0;

    angular_position[0]=0;
    angular_position[1]=0;
    angular_position[2]=0;

    omega[0]=0;
    omega[1]=0;
    omega[2]=0;

    if (dim == 3) {
        forces[2] = 0;
        velocity[2] = 0;
        position[2] = 0;
    }

    last_forces=forces;
    last_position=position;
    last_velocity=velocity;
    velocity_iter=velocity;
    last_omega=omega;
    omega_iter=omega;
    last_angular_position=angular_position;
}
template <int dim>
void
IBParticle<dim>::initialise_last() {

    last_forces=forces;
    last_position=position;
    last_velocity=velocity;
    velocity_iter=velocity;
    last_omega=omega;
    omega_iter=omega;
    last_angular_position=angular_position;

}



extern template class IBParticle<2>;
extern template class IBParticle<3>;

