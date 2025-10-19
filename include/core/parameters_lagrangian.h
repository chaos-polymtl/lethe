// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_parameters_lagrangian_h
#define lethe_parameters_lagrangian_h

#include <core/parameters.h>

#include <deal.II/base/parameter_handler.h>

#include <string>

using namespace dealii;
namespace Parameters
{
  namespace Lagrangian
  {
    enum class ParticleParticleContactForceModel : std::uint8_t
    {
      linear,
      hertz_mindlin_limit_force,
      hertz_mindlin_limit_overlap,
      hertz,
      hertz_JKR,
      DMT
    };

    enum class ParticleWallContactForceModel : std::uint8_t
    {
      linear,
      nonlinear,
      JKR,
      DMT
    };

    enum RollingResistanceMethod : std::uint8_t
    {
      none,
      constant,
      viscous,
      epsd
    };

    enum class SizeDistributionType : std::uint8_t
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
      std::vector<unsigned int> seed_for_distributions;

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

      // Hamaker constant of each particle type
      std::unordered_map<unsigned int, double> hamaker_constant_particle;

      // Coefficients of restitution of each particle type
      std::unordered_map<unsigned int, double> restitution_coefficient_particle;

      // Friction coefficient of each particle type
      std::unordered_map<unsigned int, double> friction_coefficient_particle;

      // Rolling viscous damping coefficient of each particle type
      std::unordered_map<unsigned int, double>
        rolling_viscous_damping_coefficient_particle;

      // Rolling friction coefficient of each particle type
      std::unordered_map<unsigned int, double>
        rolling_friction_coefficient_particle;

      // Thermal conductivity of each particle type
      std::unordered_map<unsigned int, double> thermal_conductivity_particle;

      // Specific heat of each particle type
      std::unordered_map<unsigned int, double> specific_heat_particle;

      // Microhardness of each particle type
      std::unordered_map<unsigned int, double> microhardness_particle;

      // Surface slope of each particle type
      std::unordered_map<unsigned int, double> surface_slope_particle;

      // Surface roughness of each particle type
      std::unordered_map<unsigned int, double> surface_roughness_particle;

      // Thermal accommodation coefficient of each particle type
      std::unordered_map<unsigned int, double> thermal_accommodation_particle;

      // Real Young's modulus of each particle type
      std::unordered_map<unsigned int, double> real_youngs_modulus_particle;

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

      // Rolling friction coefficient wall
      double rolling_viscous_damping_wall;

      // Surface energy wall
      double surface_energy_wall;

      // Hamaker constant wall
      double hamaker_constant_wall;

      // Thermal conductivity wall
      double thermal_conductivity_wall;

      // Microhardness wall
      double microhardness_wall;

      // Surface slope wall
      double surface_slope_wall;

      // Surface roughness wall
      double surface_roughness_wall;

      // Thermal accommodation wall
      double thermal_accommodation_wall;

      // Real Young's modulus of wall
      double real_youngs_modulus_wall;

      // Thermal conductivity of interstitial gas
      double thermal_conductivity_gas;

      // Specific heat of interstitial gas
      double specific_heat_gas;

      // Dynamic viscosity of interstitial gas
      double dynamic_viscosity_gas;

      // Specific heats ratio of interstitial gas
      double specific_heats_ratio_gas;

      // Molecular mean free path of interstitial gas
      double molecular_mean_free_path_gas;

      void
      declare_parameters(ParameterHandler &prm) const;
      void
      parse_parameters(ParameterHandler &prm);
      static void
      declareDefaultEntry(ParameterHandler &prm);
      void
      parse_particle_properties(const unsigned int     &particle_type,
                                const ParameterHandler &prm);

    private:
      unsigned int particle_type_maximum_number = 5;

      /**
       * @brief initialize_containers - Initialize the containers
       * used to store the physical properties of each particle type.
       *
       * @param[in,out] p_average_diameter Average diameter of each particle
       * type.
       * @param[in,out] p_size_std Size standard deviation of each particle
       * type.
       * @param[in,out] dist_types Distribution type of each particle type.
       * @param[in,out] p_custom_diameter List of diameters for the custom
       * distribution for each particle type.
       * @param[in,out] p_custom_probability Probability of each diameter value
       * based on volume fraction for the custom distribution for each particle
       * type.
       * @param[in,out] seed_for_dist Random seed for the size distribution.
       * @param[in,out] p_number Number of each particle type.
       * @param[in,out] p_density Density of each particle type.
       * @param[in,out] p_youngs_modulus Young's modulus of each particle type.
       * @param[in,out] p_poisson_ratio Poisson's ratio of each particle type.
       * @param[in,out] p_restitution_coefficient Coefficients of restitution
       * of each particle type.
       * @param[in,out] p_friction_coefficient Friction coefficient of each
       * particle type.
       * @param[in,out] p_rolling_viscous_damping_coefficient Rolling viscous
       * damping coefficient of each particle type.
       * @param[in,out] p_rolling_friction_coefficient Rolling friction
       * coefficient of each particle type.
       * @param[in,out] p_surface_energy Surface energy of each particle type.
       * @param[in,out] p_hamaker_constant Hamaker constant of each particle
       * type.
       * @param[in,out] p_thermal_conductivity Thermal conductivity of each
       * particle type.
       * @param[in,out] p_specific_heat Specific heat of each particle type.
       * @param[in,out] p_microhardness Microhardness of each particle type.
       * @param[in,out] p_surface_slope Surface slope of each particle type.
       * @param[in,out] p_surface_roughness Surface roughness of each particle
       * type.
       * @param[in,out] p_thermal_accommodation Thermal accommodation
       * coefficient of each particle type.
       * @param[in,out] p_real_youngs_modulus Real Young's modulus of each
       * particle type.
       */
      void
      initialize_containers(
        std::unordered_map<unsigned int, double> &p_average_diameter,
        std::unordered_map<unsigned int, double> &p_size_std,
        std::vector<SizeDistributionType>        &dist_types,
        std::unordered_map<unsigned int, std::vector<double>>
          &p_custom_diameter,
        std::unordered_map<unsigned int, std::vector<double>>
                                                 &p_custom_probability,
        std::vector<unsigned int>                &seed_for_dist,
        std::unordered_map<unsigned int, int>    &p_number,
        std::unordered_map<unsigned int, double> &p_density,
        std::unordered_map<unsigned int, double> &p_youngs_modulus,
        std::unordered_map<unsigned int, double> &p_poisson_ratio,
        std::unordered_map<unsigned int, double> &p_restitution_coefficient,
        std::unordered_map<unsigned int, double> &p_friction_coefficient,
        std::unordered_map<unsigned int, double>
          &p_rolling_viscous_damping_coefficient,
        std::unordered_map<unsigned int, double>
          &p_rolling_friction_coefficient,
        std::unordered_map<unsigned int, double> &p_surface_energy,
        std::unordered_map<unsigned int, double> &p_hamaker_constant,
        std::unordered_map<unsigned int, double> &p_thermal_conductivity,
        std::unordered_map<unsigned int, double> &p_specific_heat,
        std::unordered_map<unsigned int, double> &p_microhardness,
        std::unordered_map<unsigned int, double> &p_surface_slope,
        std::unordered_map<unsigned int, double> &p_surface_roughness,
        std::unordered_map<unsigned int, double> &p_thermal_accommodation,
        std::unordered_map<unsigned int, double> &p_real_youngs_modulus) const;
    };

    template <int dim>
    class InsertionInfo
    {
    public:
      // Insertion method
      enum class InsertionMethod
      {
        file,
        list,
        plane,
        volume
      } insertion_method;

      // Inserted number of particles at each time step
      int inserted_this_step;

      // Insertion frequency
      int insertion_frequency;

      /* Removal box: */
      bool removing_particles_in_region;

      // Clear box info (xmin,xmax,ymin,ymax,zmin,zmax)
      Point<3> clear_box_point_1, clear_box_point_2;

      /* File: */
      std::vector<std::string> list_of_input_files;

      /* Plane: */
      // Plane normal vector
      Tensor<1, 3> insertion_plane_normal_vector;
      // Plane point
      Point<3> insertion_plane_point;

      /* List */
      // Containers used for the list insertion method
      std::vector<double> list_x, list_y, list_z, list_vx, list_vy, list_vz,
        list_wx, list_wy, list_wz, list_d, list_T;

      /* Volume */
      // Direction sequence for the insertion of particles (1st, 2nd, 3rd)
      std::vector<unsigned int> direction_sequence;
      // Insertion box info (p_1 , p_2)
      Point<3> insertion_box_point_1, insertion_box_point_2;
      // Insertion initial velocity conditions
      Tensor<1, 3> initial_vel, initial_omega;
      // Function that returns the initial temperature of a particle based on
      // time or its position.
      std::shared_ptr<Function<dim>> initial_temperature_function;
      // Insertion distance threshold
      double distance_threshold;
      // Insertion random number range
      double insertion_maximum_offset;
      // Insertion random number seed
      int seed_for_insertion;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

    template <int dim>
    class ModelParameters
    {
    public:
      // Load balance method
      enum class LoadBalanceMethod
      {
        none,
        once,
        frequent,
        dynamic,
        dynamic_with_sparse_contacts
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

      // Function that returns the weight of a cell base on its barycenter
      // position.
      std::shared_ptr<Function<dim>> cell_weight_function;

      // The particle weight for load balancing
      unsigned int load_balance_particle_weight;

      // Factors applied on the particle weight in load balancing for active and
      // inactive cells (factor of mobile cells is always 1), only available
      // when adaptive sparse contacts is enable
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

      // Cut-off threshold where Van der Waals forces are ignored.
      double dmt_cut_off_threshold;

      // Particle-particle contact force model
      ParticleParticleContactForceModel particle_particle_contact_force_model;

      // Particle-wall contact force model
      ParticleWallContactForceModel particle_wall_contact_force_method;

      // Rolling resistance torque method
      RollingResistanceMethod rolling_resistance_method;

      // Model parameter for the EPSD rolling resistance model
      double f_coefficient_epsd;

      // Integration method
      enum class IntegrationMethod
      {
        velocity_verlet,
        explicit_euler
      } integration_method;

      // Solver type
      DEM::SolverType solver_type;

      // Disable particle contacts to optimize performance
      bool sparse_particle_contacts;

      // Enable advection of particles (applies cell average velocity and
      // acceleration to particles)
      bool advect_particles;

      // Minimal granular temperature value of cells where particle contacts
      // are considered
      double granular_temperature_threshold;

      // Maximal solid fraction value of cells where particle contacts are
      // considered no matter the granular temperature
      double solid_fraction_threshold;

      // Disable position integration
      bool disable_position_integration;

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

      static void
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
      declare_parameters(ParameterHandler &prm) const;
      void
      parse_parameters(ParameterHandler &prm);
      static void
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
      };

      // Vector of each boundary types
      std::vector<BoundaryType> bc_types;

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
      declare_parameters(ParameterHandler &prm) const;
      void
      parse_parameters(ParameterHandler &prm);
      static void
      declareDefaultEntry(ParameterHandler &prm);
      void
      parse_boundary_conditions(const ParameterHandler &prm);

    private:
      unsigned int DEM_BC_number_max = 10;
      void
      initialize_containers(
        std::unordered_map<unsigned int, Tensor<1, 3>> &boundary_trans_velocity,
        std::unordered_map<unsigned int, double>       &boundary_rot_speed,
        std::unordered_map<unsigned int, Tensor<1, 3>> &boundary_rot_vector,
        std::unordered_map<unsigned int, Point<3>>     &point_on_rot_axis,
        std::vector<unsigned int>                      &outlet_boundaries_id,
        std::vector<BoundaryType>                      &boundaries_types) const;
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
      bool lagrangian_post_processing_enabled;
      /// A bool variable which sets-up the force chains visualization
      bool force_chains;

      // Enable the logging of particle-wall contact statistics
      bool particle_wall_collision_statistics;

      // State whether collisions with all walls should be logged
      bool log_collisions_with_all_walls;

      // Boundary ids of the walls to log collisions with
      std::vector<int> particle_wall_collision_boundary_ids;

      Parameters::Verbosity collision_verbosity;

      // Exporting collision statistics csv filename
      std::string collision_stats_file_name;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

    template <int dim>
    class FloatingGrid
    {
    public:
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

    /**
     * @brief Defines parameters used to construct particle
     * ray tracing class.
     */
    template <int dim>
    struct ParticleRayTracing
    {
      // Location of the first photon to be inserted
      Point<3> starting_point;

      // In which directions will the photons be inserted relative to the first
      // photon.
      std::vector<Tensor<1, 3>> insertion_directions_units_vector;

      // How many photon will be inserted in each of those directions.
      std::vector<unsigned int> n_photons_each_directions;

      // What is the distance between each photon in each of those directions
      // considering an offset equal to 0.
      std::vector<double> step_between_photons_each_directions;

      // Reference displacement unit tensor
      Tensor<1, 3> ref_displacement_tensor_unit;

      // Related to the offset insertion position.
      double       max_insertion_offset;
      unsigned int prn_seed_photon_insertion;

      // Related to the offset in the displacement direction.
      double       max_angular_offset;
      unsigned int prn_seed_photon_displacement;

      // Declare and parse function
      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

  } // namespace Lagrangian
} // namespace Parameters
#endif /* lethe_parameters_lagrangian_h */
