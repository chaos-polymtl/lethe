#include <dem/dem_properties.h>
#include <dem/rolling_resistance_torque_models.h>


using namespace DEM;


template <int dim>
Tensor<1, dim>
NoRollingResistanceTorque<dim>::calculate_rolling_resistance_torque(
  const double & /*effective_r*/,
  const ArrayView<const double> & /*particle_one_properties*/,
  const ArrayView<const double> & /*particle_two_properties*/,
  const double & /*effective_rolling_friction_coefficient*/,
  const double & /*normal_force_norm*/,
  const Tensor<1, dim> & /*normal_contact_vector*/)
{
  Tensor<1, dim> rolling_resistance;
  for (int d = 0; d < dim; ++d)
    rolling_resistance[d] = 0;

  return rolling_resistance;
}

template class NoRollingResistanceTorque<2>;
template class NoRollingResistanceTorque<3>;


template <int dim>
Tensor<1, dim>
ConstantRollingResistanceTorque<dim>::calculate_rolling_resistance_torque(
  const double &                 effective_r,
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties,
  const double &                 effective_rolling_friction_coefficient,
  const double &                 normal_force_norm,
  const Tensor<1, dim> & /*normal_contact_vector*/)
{
  // For calculation of rolling resistance torque, we need to obtain
  // omega_ij using rotational velocities of particles one and two
  Tensor<1, dim> particle_one_angular_velocity, particle_two_angular_velocity,
    omega_ij, omega_ij_direction;
  for (int d = 0; d < dim; ++d)
    {
      particle_one_angular_velocity[d] =
        particle_one_properties[DEM::PropertiesIndex::omega_x + d];
      particle_two_angular_velocity[d] =
        particle_two_properties[DEM::PropertiesIndex::omega_x + d];
    }

  omega_ij = particle_one_angular_velocity - particle_two_angular_velocity;
  omega_ij_direction = omega_ij / (omega_ij.norm() + DBL_MIN);

  // Calculation of rolling resistance torque
  Tensor<1, dim> rolling_resistance_torque =
    -effective_rolling_friction_coefficient * effective_r * normal_force_norm *
    omega_ij_direction;

  return rolling_resistance_torque;
}

template class ConstantRollingResistanceTorque<2>;
template class ConstantRollingResistanceTorque<3>;


template <int dim>
Tensor<1, dim>
ViscousRollingResistanceTorque<dim>::calculate_rolling_resistance_torque(
  const double &                 effective_r,
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties,
  const double &                 effective_rolling_friction_coefficient,
  const double &                 normal_force_norm,
  const Tensor<1, dim> &         normal_contact_vector)
{
  // For calculation of rolling resistance torque, we need to obtain
  // omega_ij using rotational velocities of particles one and two
  Tensor<1, dim> particle_one_angular_velocity, particle_two_angular_velocity,
    omega_ij, omega_ij_direction;
  for (int d = 0; d < dim; ++d)
    {
      particle_one_angular_velocity[d] =
        particle_one_properties[DEM::PropertiesIndex::omega_x + d];
      particle_two_angular_velocity[d] =
        particle_two_properties[DEM::PropertiesIndex::omega_x + d];
    }

  omega_ij = particle_one_angular_velocity - particle_two_angular_velocity;
  omega_ij_direction = omega_ij / (omega_ij.norm() + DBL_MIN);

  Tensor<1, dim> v_omega =
    cross_product_3d(particle_one_angular_velocity,
                     particle_one_properties[DEM::PropertiesIndex::dp] * 0.5 *
                       normal_contact_vector) -
    cross_product_3d(particle_two_angular_velocity,
                     particle_two_properties[DEM::PropertiesIndex::dp] * 0.5 *
                       -normal_contact_vector);

  // Calculation of rolling resistance torque
  Tensor<1, dim> rolling_resistance_torque =
    -effective_rolling_friction_coefficient * effective_r * normal_force_norm *
    v_omega.norm() * omega_ij_direction;

  return rolling_resistance_torque;
}

template class ViscousRollingResistanceTorque<2>;
template class ViscousRollingResistanceTorque<3>;
