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

//#include <deal.II/base/conditional_ostream.h>
//#include <deal.II/base/function.h>
#include <core/parameters.h>

#include <deal.II/base/parameter_handler.h>

#include <string>

using namespace dealii;
namespace Parameters
{
  namespace Lagrangian
  {
    template <int dim>
    class LagrangianPhysicalProperties
    {
    public:
      // Gravitational acceleration
      Tensor<1, dim> g;

      // Choosing integration method
      enum class size_distribution_type
      {
        uniform,
        normal
      } size_distribution_type;

      // Number of particle types
      unsigned int particle_type_number;

      // Average diameter of each particle type
      std::unordered_map<unsigned int, double> particle_average_diameter;

      // Size standard deviation of each particle type
      std::unordered_map<unsigned int, double> particle_size_std;

      // Number of each particle type
      std::unordered_map<unsigned int, int> number;

      // Density of each particle type
      std::unordered_map<unsigned int, double> density_particle;

      // Young's modulus of each particle type
      std::unordered_map<unsigned int, double> youngs_modulus_particle;

      // Poisson's ratio of each particle type
      std::unordered_map<unsigned int, double> poisson_ratio_particle;

      // Coefficients of restituion of each particle type
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

      // Coefficient of restituion of wall
      double restitution_coefficient_wall;

      // Friction coefficient of wall
      double friction_coefficient_wall;

      // Rolling friction coefficient wall
      double rolling_friction_wall;

      void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
      void
      declareDefaultEntry(ParameterHandler &prm);
      void
      parse_particle_properties(const unsigned int &particle_type,
                                ParameterHandler &  prm);

    private:
      unsigned int particle_type_maximum_number = 5;

      void
      initialize_containers(
        std::unordered_map<unsigned int, double> &particle_average_diameter,
        std::unordered_map<unsigned int, double> &particle_size_std,
        std::unordered_map<unsigned int, int> &   number,
        std::unordered_map<unsigned int, double> &density_particle,
        std::unordered_map<unsigned int, double> &youngs_modulus_particle,
        std::unordered_map<unsigned int, double> &poisson_ratio_particle,
        std::unordered_map<unsigned int, double>
          &restitution_coefficient_particle,
        std::unordered_map<unsigned int, double> &friction_coefficient_particle,
        std::unordered_map<unsigned int, double>
          &rolling_friction_coefficient_particle);
    };

    struct InsertionInfo
    {
      // Choosing insertion method
      enum class InsertionMethod
      {
        uniform,
        non_uniform,
        list
      } insertion_method;

      // Inserted number of particles at each time step
      int inserted_this_step;

      // Insertion frequency
      int insertion_frequency;

      // Insertion box info (xmin,xmax,ymin,ymax,zmin,zmax)
      double x_min, y_min, z_min, x_max, y_max, z_max;

      // Insertion distance threshold
      double distance_threshold;

      // Insertion random number range
      double random_number_range;

      // Insertion random number seed
      int random_number_seed;

      std::vector<double> list_x, list_y, list_z;

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
        dynamic
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

      // Choosing particle-particle contact force model
      enum class PPContactForceModel
      {
        pp_linear,
        pp_nonlinear
      } pp_contact_force_method;

      // Choosing particle-wall contact force model
      enum class PWContactForceModel
      {
        pw_linear,
        pw_nonlinear
      } pw_contact_force_method;

      // Rolling resistance torque method
      enum class RollingResistanceMethod
      {
        no_resistance,
        constant_resistance,
        viscous_resistance
      } rolling_resistance_method;

      // Choosing integration method
      enum class IntegrationMethod
      {
        velocity_verlet,
        explicit_euler,
        gear3
      } integration_method;

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

      // Choosing how the outputs is gonna be displayed
      Parameters::Verbosity force_torque_verbosity;

      // Output frequency
      unsigned int output_frequency;

      // Prefix for simulation output
      std::string force_torque_output_name;

      // Center of mass
      Point<dim> point_center_mass;

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

    template <int dim>
    class BCDEM
    {
    public:
      // Number of DEM boundary conditions
      unsigned int DEM_BC_number;

      // Choosing BC type
      enum class BoundaryType
      {
        fixed_wall,
        outlet,
        translational,
        rotational
      } BC_type;

      // Outlet boundary IDs
      std::vector<unsigned int> outlet_boundaries;

      // Translational velocities of moving boundaries
      std::unordered_map<unsigned int, Tensor<1, dim>>
        boundary_translational_velocity;

      // Rotational speeds of rotating boundaries in rad/s
      std::unordered_map<unsigned int, double> boundary_rotational_speed;

      // Rotational axes of rotating boundaries
      std::unordered_map<unsigned int, Tensor<1, dim>>
        boundary_rotational_vector;

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
        std::unordered_map<unsigned int, Tensor<1, dim>>
          &boundary_translational_velocity,
        std::unordered_map<unsigned int, double> &boundary_rotational_speed,
        std::unordered_map<unsigned int, Tensor<1, dim>>
          &                        boundary_rotational_vector,
        std::vector<unsigned int> &outlet_boundaries);
    };

    template <int dim>
    class GridMotion
    {
    public:
      // Choosing grid motion type
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

      // Enable particles velocity post-processing
      bool calculate_particles_average_velocity;

      // Enable granular temperature post-processing
      bool calculate_granular_temperature;

      // Set initial step to start post-processing calculations
      unsigned int initial_step;

      // Set end step to finish post-processing calculations
      unsigned int end_step;

      // Set post-processing output frequency
      unsigned int output_frequency;

      // Prefix for particles velocity output
      std::string particles_velocity_name;

      // Prefix for granular temperature output
      std::string granular_temperature_name;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

  } // namespace Lagrangian
} // namespace Parameters

#endif /* PARAMETERS_H_ */
