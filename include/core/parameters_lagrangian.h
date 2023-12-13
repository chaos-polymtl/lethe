/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Shahab Golshan, Bruno Blais, Polytechnique Montreal, 2019-
 */

#ifndef lethe_parameters_lagrangian_h
#define lethe_parameters_lagrangian_h

// #include <deal.II/base/conditional_ostream.h>
// #include <deal.II/base/function.h>
#include <core/parameters.h>

#include <deal.II/base/parameter_handler.h>

#include <string>

using namespace dealii;
namespace Parameters
{
  namespace Lagrangian
  {
    enum ParticleParticleContactForceModel
    {
      linear,
      hertz_mindlin_limit_force,
      hertz_mindlin_limit_overlap,
      hertz,
      hertz_JKR
    };

    enum RollingResistanceMethod
    {
      no_resistance,
      constant_resistance,
      viscous_resistance
    };

    enum class SizeDistributionType
    {
      uniform,
      normal,
      custom
    };

    struct LagrangianPhysicalProperties
    {
    public:
      // Gravitational acceleration
      Tensor<1, 3> g;

      // Number of particle types
      unsigned int particle_type_number;

      // Average diameter of each particle type
      std::unordered_map<unsigned int, double> particle_average_diameter;

      // Size standard deviation of each particle type
      std::unordered_map<unsigned int, double> particle_size_std;

      // List of diameters for the custom distribution for each particle type
      std::unordered_map<unsigned int, std::vector<double>>
        particle_custom_diameter;

      // Probability of each diameter value based on volume fraction for the
      // custom distribution for each particle type
      std::unordered_map<unsigned int, std::vector<double>>
        particle_custom_probability;

      // Random seed for the size distribution
      std::vector<unsigned int> random_seed_for_distributions;

      // Distribution type of each particle type
      std::vector<SizeDistributionType> distribution_type;

      // Number of each particle type
      std::unordered_map<unsigned int, int> number;

      // Density of each particle type
      std::unordered_map<unsigned int, double> density_particle;

      // Young's modulus of each particle type
      std::unordered_map<unsigned int, double> youngs_modulus_particle;

      // Poisson's ratio of each particle type
      std::unordered_map<unsigned int, double> poisson_ratio_particle;

      // Surface energy of each particle type
      std::unordered_map<unsigned int, double> surface_energy_particle;

      // Coefficients of restitution of each particle type
      std::unordered_map<unsigned int, double> restitution_coefficient_particle;

      // Friction coefficient of each particle type
      std::unordered_map<unsigned int, double> friction_coefficient_particle;

      // Rolling friction coefficient of each particle type
      std::unordered_map<unsigned int, double>
        rolling_friction_coefficient_particle;

      // Young's modulus of wall
      double youngs_modulus_wall;

      // Poisson's ratio of wall
      double poisson_ratio_wall;

      // Coefficient of restitution of wall
      double restitution_coefficient_wall;

      // Friction coefficient of wall
      double friction_coefficient_wall;

      // Rolling friction coefficient wall
      double rolling_friction_wall;

      // Surface energy wall
      double surface_energy_wall;

      void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
      void
      declareDefaultEntry(ParameterHandler &prm);
      void
      parse_particle_properties(const unsigned int &particle_type,
                                ParameterHandler   &prm);

    private:
      unsigned int particle_type_maximum_number = 5;

      void
      initialize_containers(
        std::unordered_map<unsigned int, double> &particle_average_diameter,
        std::unordered_map<unsigned int, double> &particle_size_std,
        std::vector<SizeDistributionType>        &distribution_type,
        std::unordered_map<unsigned int, std::vector<double>>
          &particle_custom_diameter,
        std::unordered_map<unsigned int, std::vector<double>>
                                                 &particle_custom_probability,
        std::vector<unsigned int>                &random_seed_for_distributions,
        std::unordered_map<unsigned int, int>    &number,
        std::unordered_map<unsigned int, double> &density_particle,
        std::unordered_map<unsigned int, double> &youngs_modulus_particle,
        std::unordered_map<unsigned int, double> &poisson_ratio_particle,
        std::unordered_map<unsigned int, double>
          &restitution_coefficient_particle,
        std::unordered_map<unsigned int, double> &friction_coefficient_particle,
        std::unordered_map<unsigned int, double>
          &rolling_friction_coefficient_particle,
        std::unordered_map<unsigned int, double> &surface_energy_particle);
    };

    struct InsertionInfo
    {
      // Insertion method
      enum class InsertionMethod
      {
        volume,
        list,
        plane
      } insertion_method;

      // Inserted number of particles at each time step
      int inserted_this_step;

      // Insertion frequency
      int insertion_frequency;

      // Direction (axis) of insertion of particles (1st, 2nd, 3rd)
      unsigned int axis_0, axis_1, axis_2;

      // Insertion box info (xmin,xmax,ymin,ymax,zmin,zmax)
      double x_min, y_min, z_min, x_max, y_max, z_max;

      // Insertion initial conditions
      double vel_x, vel_y, vel_z, omega_x, omega_y, omega_z;

      // Insertion distance threshold
      double distance_threshold;

      // Insertion random number range
      double random_number_range;

      // Insertion random number seed
      int random_number_seed;

      std::vector<double> list_x, list_y, list_z, list_vx, list_vy, list_vz,
        list_wx, list_wy, list_wz, list_d;

      // Insertion plane definition
      Tensor<1, 3> insertion_plane_normal_vector;
      Point<3>     insertion_plane_point;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

    struct ModelParameters
    {
      // Load balance method
      enum class LoadBalanceMethod
      {
        none,
        once,
        frequent,
        dynamic,
        dynamic_with_disabling_contacts
      } load_balance_method;

      // Load balance step (for single step load-balancing)
      unsigned int load_balance_step;

      // Load balance frequency
      unsigned int load_balance_frequency;

      // Load balance threshold (for dynamic load-balancing)
      double load_balance_threshold;

      // Load balance check frequency (for dynamic load-balancing)
      unsigned int dynamic_load_balance_check_frequency;

      // Particle-particle, particle-wall broad and fine search frequency
      unsigned int contact_detection_frequency;

      // The particle weight based on a default cell weight of 1000
      unsigned int load_balance_particle_weight;

      // Factors applied on the particle weight in load balancing for active and
      // inactive cells (factor of mobile cells is always 1), only available
      // when dynamic disabling particle contacts is enable
      double active_load_balancing_factor;
      double inactive_load_balancing_factor;

      // Security factor for dynamic contact search
      double dynamic_contact_search_factor;

      // Contact detection method
      enum class ContactDetectionMethod
      {
        constant,
        dynamic
      } contact_detection_method;

      // Contact search neighborhood threshold (neighborhood diameter to
      // particle diameter)
      double neighborhood_threshold;

      // Particle-particle contact force model
      ParticleParticleContactForceModel particle_particle_contact_force_model;

      // Particle-wall contact force model
      enum class ParticleWallContactForceModel
      {
        linear,
        nonlinear,
        JKR
      } particle_wall_contact_force_method;

      // Rolling resistance torque method
      RollingResistanceMethod rolling_resistance_method;

      // Itegration method
      enum class IntegrationMethod
      {
        velocity_verlet,
        explicit_euler,
        gear3
      } integration_method;

      // Disable particle contacts to optimize performance
      bool disable_particle_contacts;

      // Minimal granular temperature value of cells where particle contacts
      // are considered
      double granular_temperature_threshold;

      // Maximal solid fraction value of cells where particle contacts are
      // considered no matter the granular temperature
      double solid_fraction_threshold;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

    /**
     * @brief ForceTorqueOnWall - Defines the parameters for the
     * force and torques calculation on boundaries of the domain.
     */
    template <int dim>
    class ForceTorqueOnWall
    {
    public:
      // Enable force post-processing
      bool calculate_force_torque;

      Parameters::Verbosity force_torque_verbosity;

      // Output frequency
      unsigned int output_frequency;

      // Prefix for simulation output
      std::string force_torque_output_name;

      // Center of mass
      Point<3> point_center_mass;

      void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

    template <int dim>
    class FloatingWalls
    {
    public:
      // Number of floating walls
      unsigned int floating_walls_number;

      // A point on each floating wall
      std::vector<Point<dim>> points_on_walls;

      // Normal vectors of the floating walls
      std::vector<Tensor<1, dim>> floating_walls_normal_vectors;

      // Beginning time
      std::vector<double> time_start;

      // Ending time
      std::vector<double> time_end;

      void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
      void
      declareDefaultEntry(ParameterHandler &prm);
      void
      parse_floating_wall(ParameterHandler &prm);

    private:
      unsigned int max_number_floating_walls = 9;
    };


    struct BCDEM
    {
    public:
      // Number of DEM boundary conditions
      unsigned int DEM_BC_number;

      // Boundary condition type
      enum class BoundaryType
      {
        fixed_wall,
        outlet,
        translational,
        rotational,
        periodic
      } BC_type;

      // Outlet boundary IDs
      std::vector<unsigned int> outlet_boundaries;

      // Translational velocities of moving boundaries
      std::unordered_map<unsigned int, Tensor<1, 3>>
        boundary_translational_velocity;

      // Rotational speeds of rotating boundaries in rad/s
      std::unordered_map<unsigned int, double> boundary_rotational_speed;

      // Rotational axes of rotating boundaries
      std::unordered_map<unsigned int, Tensor<1, 3>> boundary_rotational_vector;

      // Point on rotational axis
      std::unordered_map<unsigned int, Point<3>> point_on_rotation_axis;

      // Periodic boundary IDs
      types::boundary_id periodic_boundary_0;
      types::boundary_id periodic_boundary_1;
      types::boundary_id periodic_direction;


      void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
      void
      declareDefaultEntry(ParameterHandler &prm);
      void
      parse_boundary_conditions(ParameterHandler &prm);

    private:
      unsigned int DEM_BC_number_max = 10;
      void
      initialize_containers(
        std::unordered_map<unsigned int, Tensor<1, 3>>
          &boundary_translational_velocity,
        std::unordered_map<unsigned int, double> &boundary_rotational_speed,
        std::unordered_map<unsigned int, Tensor<1, 3>>
                                                   &boundary_rotational_vector,
        std::unordered_map<unsigned int, Point<3>> &point_on_rotation_axis,
        std::vector<unsigned int>                  &outlet_boundaries);
    };

    template <int dim>
    class GridMotion
    {
    public:
      // Grid motion type
      enum class MotionType
      {
        translational,
        rotational,
        none
      } motion_type;

      // Translational velocity of the moving grid
      Tensor<1, dim> grid_translational_velocity;

      // Rotational speed of rotating grid in rad/s
      double grid_rotational_speed;

      // Rotational axis of rotating grid. Similar to deal.II, we use 0=x axis,
      // 1=y axis, 2=z axis.
      unsigned int grid_rotational_axis;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

    /**
     * @brief Lagrangian Postprocessing - Defines the parameters
     * for the Lagrangian post-processing.
     *
     */
    struct LagrangianPostProcessing
    {
      // A bool variable which sets-up the Lagrangian post-processing
      bool Lagrangian_post_processing;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

    template <int dim>
    struct FloatingGrid
    {
      // Floating mesh motion information
      Parameters::Mesh mesh;

      // Floating grid motion information
      Parameters::Lagrangian::GridMotion<dim> motion;

      // Beginning time
      double time_start;
      // Ending time
      double time_end;

      void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

  } // namespace Lagrangian
} // namespace Parameters
#endif /* lethe_parameters_lagrangian_h */
