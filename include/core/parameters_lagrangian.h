// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file parameters_lagrangian.h
 * @brief Parameter structures for Lagrangian (DEM) particle simulations.
 *
 * This file defines the parameter classes and enumerations used to configure
 * Lagrangian particle simulations in Lethe, including contact force models,
 * particle physical properties, insertion methods, model parameters,
 * floating walls, boundary conditions, grid motion, and post-processing.
 */

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
    /**
     * @brief Contact force model for particle-particle interactions.
     */
    enum class ParticleParticleContactForceModel : std::uint8_t
    {
      /// Linear model
      linear,
      /// Hertz-Mindlin model with a limitation on tangential force
      hertz_mindlin_limit_force,
      /// Hertz-Mindlin model with a limitation on tangential overlap
      hertz_mindlin_limit_overlap,
      /// Pure Hertz model
      hertz,
      /// Hertz model with Johnson-Kendall-Roberts cohesive force
      hertz_JKR,
      /// Derjaguin-Muller-Toporov model
      DMT
    };

    /**
     * @brief Contact force model for particle-wall interactions.
     */
    enum class ParticleWallContactForceModel : std::uint8_t
    {
      /// Linear model
      linear,
      /// Non-linear force model that corresponds to a Hertz-Mindlin with
      /// limitation on tangential overlap
      nonlinear,
      /// Johnson-Kendall-Roberts model
      JKR,
      /// Derjaguin-Muller-Toporov model
      DMT
    };

    /**
     * @brief Method used to compute the rolling resistance torque.
     */
    enum RollingResistanceMethod : std::uint8_t
    {
      /// No rolling resistance
      none,
      /// Constant rolling resistance
      constant,
      /// Viscous rolling resistance
      viscous,
      /// Elastic-plastic spring-dashpot
      epsd
    };

    /**
     * @brief Type of particle size distribution.
     */
    enum class SizeDistributionType : std::uint8_t
    {
      /// Uniform distribution
      uniform,
      /// Normal distribution
      normal,
      /// Log-normal distribution
      lognormal,
      /// Custom distribution
      custom
    };

    /**
     * @brief Weighting basis for particle size distributions.
     */
    enum class DistributionWeightingType : std::uint8_t
    {
      /// Number-based distribution weighting
      number_based,
      /// Volume-based distribution weighting
      volume_based
    };

    enum class ProbabilityFunctionType : std::uint8_t
    {
      /// Probability density function
      PDF,
      /// Cumulative density function
      CDF
    };

    /**
     * @brief Physical properties of particles, walls, and interstitial gas
     * for Lagrangian (DEM) simulations.
     *
     * This structure stores all material properties needed for contact force
     * computations, including per-particle-type mechanical and thermal
     * properties, wall properties, and interstitial gas properties.
     */
    struct LagrangianPhysicalProperties
    {
    public:
      /// Gravitational acceleration vector.
      Tensor<1, 3> g;

      /// Number of particle types.
      unsigned int particle_type_number;

      // Distribution type of each particle type (uniform, normal, lognormal,
      // custom)
      std::vector<SizeDistributionType> distribution_type;

      /// Average diameter of each particle type
      std::vector<double> particle_average_diameter;

      /// Size standard deviation of each particle type
      std::vector<double> particle_size_std;

      /// Indicate if the custom distribution is read from file
      std::vector<bool> custom_distribution_from_file;

      /// Filename for the custom distribution for each particle type
      std::vector<std::string> custom_distribution_filenames;

      /// Custom distribution function type (PDF or CDF)
      std::vector<ProbabilityFunctionType> custom_probability_function_type;

      /// Indicates whether the diameter values generated from the custom
      /// distribution are interpolated
      std::vector<bool> custom_distribution_interpolation;

      /// List of diameters for the custom distribution for each particle type
      std::vector<std::vector<double>> particle_custom_diameter;

      /// Probability of each diameter value based on volume fraction for the
      /// custom distribution for each particle type
      std::vector<std::vector<double>> particle_custom_probability;

      /// Distribution weighting type of each particle type (number-based,
      /// volume-based)
      std::vector<DistributionWeightingType> distribution_weighting_type;

      /// Random seed for the size distribution.
      std::vector<unsigned int> seed_for_distributions;

      /// Minimum diameter cutoff for lognormal distribution.
      std::vector<double> diameter_min_cutoff;

      /// Maximum diameter cutoff for lognormal distribution.
      std::vector<double> diameter_max_cutoff;

      /// Number particles of each particle type
      std::vector<int> number;

      /// Density of each particle type
      std::vector<double> density_particle;

      /// Young's modulus of each particle type
      std::vector<double> youngs_modulus_particle;

      /// Poisson's ratio of each particle type
      std::vector<double> poisson_ratio_particle;

      /// Coefficients of restitution of each particle type
      std::vector<double> restitution_coefficient_particle;

      /// Friction coefficient of each particle type
      std::vector<double> friction_coefficient_particle;

      /// Rolling viscous damping coefficient of each particle type
      std::vector<double> rolling_viscous_damping_coefficient_particle;

      /// Rolling friction coefficient of each particle type
      std::vector<double> rolling_friction_coefficient_particle;

      /// Surface energy of each particle type
      std::vector<double> surface_energy_particle;

      /// Hamaker constant of each particle type
      std::vector<double> hamaker_constant_particle;

      /// Thermal conductivity of each particle type
      std::vector<double> thermal_conductivity_particle;

      /// Specific heat of each particle type
      std::vector<double> specific_heat_particle;

      /// Microhardness of each particle type
      std::vector<double> microhardness_particle;

      /// Surface slope of each particle type
      std::vector<double> surface_slope_particle;

      /// Surface roughness of each particle type
      std::vector<double> surface_roughness_particle;

      /// Thermal accommodation coefficient of each particle type
      std::vector<double> thermal_accommodation_particle;

      /// Real Young's modulus of each particle type
      std::vector<double> real_youngs_modulus_particle;

      /// Young's modulus of the wall.
      double youngs_modulus_wall;

      /// Poisson's ratio of the wall.
      double poisson_ratio_wall;

      /// Coefficient of restitution of the wall.
      double restitution_coefficient_wall;

      /// Friction coefficient of the wall.
      double friction_coefficient_wall;

      /// Rolling friction coefficient of the wall.
      double rolling_friction_wall;

      /// Rolling viscous damping coefficient of the wall.
      double rolling_viscous_damping_wall;

      /// Surface energy of the wall.
      double surface_energy_wall;

      /// Hamaker constant of the wall.
      double hamaker_constant_wall;

      /// Thermal conductivity of the wall.
      double thermal_conductivity_wall;

      /// Microhardness of the wall.
      double microhardness_wall;

      /// Surface slope of the wall.
      double surface_slope_wall;

      /// Surface roughness of the wall.
      double surface_roughness_wall;

      /// Thermal accommodation coefficient of the wall.
      double thermal_accommodation_wall;

      /// Real Young's modulus of the wall.
      double real_youngs_modulus_wall;

      /// Thermal conductivity of the interstitial gas.
      double thermal_conductivity_gas;

      /// Specific heat of the interstitial gas.
      double specific_heat_gas;

      /// Dynamic viscosity of the interstitial gas.
      double dynamic_viscosity_gas;

      /// Specific heats ratio of the interstitial gas.
      double specific_heats_ratio_gas;

      /// Molecular mean free path of the interstitial gas.
      double molecular_mean_free_path_gas;

      /**
       * @brief Declare the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      declare_parameters(ParameterHandler &prm) const;

      /**
       * @brief Parse the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      parse_parameters(ParameterHandler &prm);

      /**
       * @brief Declare default parameter entries for a single particle type.
       *
       * @param[in,out] prm The parameter handler.
       */
      static void
      declareDefaultEntry(ParameterHandler &prm);

      /**
       * @brief Parse the physical properties for a given particle type.
       *
       * @param[in] particle_type Index of the particle type to parse.
       * @param[in, out] prm The parameter handler.
       */
      void
      parse_particle_properties(const unsigned int     &particle_type,
                                const ParameterHandler &prm);

    private:
      unsigned int particle_type_maximum_number = 5;

      /**
       * @brief initialize_containers - Initialize the containers
       * used to store the particle size distribution and physical properties
       * of each particle type.
       *
       * @param[in,out] dist_types Indicates the SizeDistributionType type.
       * @param[in,out] p_average_diameter Average diameter.
       * @param[in,out] p_size_std Diameter standard deviation for the
       * normal and lognormal distributions.
       * @param[in,out] custom_dist_read_from_file Indicates whether the
       * diameter and probability values are stored in a separate file.
       * @param[in,out] custom_dist_file_names Vector of filenames from which
       * the custom distributions are read.
       * @param[in,out] custom_function_type Indicates whether the values
       * defining the custom distribution correspond to the PDF or the CDF.
       * @param[in,out] custom_interpolation Indicates whether the diameter
       * values are interpolated from the input values.
       * @param[in,out] custom_diameter_values Lists of diameters for the custom
       * distribution.
       * @param[in,out] custom_probabilities_values Lists of probability values
       * corresponding to each diameter value in the custom distribution.
       * @param[in,out] distribution_weighting_basis_type Weighting basis for
       * the PSD.
       * @param[in,out] seed_for_dist Pseudo random seed for the PSD sampling.
       * @param[in,out] dia_min_cutoff Minimal cutoff value when sampling a PSD.
       * @param[in,out] dia_max_cutoff Maximal cutoff value when sampling a PSD.
       * @param[in,out] p_number Number of particles.
       * @param[in,out] p_density Density.
       * @param[in,out] p_youngs_modulus Young's modulus.
       * @param[in,out] p_poisson_ratio Poisson's ratio.
       * @param[in,out] p_restitution_coefficient Coefficient of restitution.
       * @param[in,out] p_friction_coefficient Sliding friction coefficient.
       * @param[in,out] p_rolling_viscous_damping_coefficient Rolling viscous
       * damping coefficient
       * @param[in,out] p_rolling_friction_coefficient Rolling friction
       * coefficient.
       * @param[in,out] p_surface_energy Surface energy.
       * @param[in,out] p_hamaker_constant Hamaker constant.
       * @param[in,out] p_thermal_conductivity Thermal conductivity.
       * @param[in,out] p_specific_heat Specific heat.
       * @param[in,out] p_microhardness Microhardness.
       * @param[in,out] p_surface_slope Surface slope.
       * @param[in,out] p_surface_roughness Surface roughness.
       * @param[in,out] p_thermal_accommodation Thermal accommodation
       * coefficient.
       * @param[in,out] p_real_youngs_modulus Real Young's modulus.
       */
      void
      initialize_containers(
        std::vector<SizeDistributionType>    &dist_types,
        std::vector<double>                  &p_average_diameter,
        std::vector<double>                  &p_size_std,
        std::vector<bool>                    &custom_dist_read_from_file,
        std::vector<std::string>             &custom_dist_file_names,
        std::vector<ProbabilityFunctionType> &custom_function_type,
        std::vector<bool>                    &custom_interpolation,
        std::vector<std::vector<double>>     &custom_diameter_values,
        std::vector<std::vector<double>>     &custom_probabilities_values,
        std::vector<DistributionWeightingType>
                                  &distribution_weighting_basis_type,
        std::vector<unsigned int> &seed_for_dist,
        std::vector<double>       &dia_min_cutoff,
        std::vector<double>       &dia_max_cutoff,
        std::vector<int>          &p_number,
        std::vector<double>       &p_density,
        std::vector<double>       &p_youngs_modulus,
        std::vector<double>       &p_poisson_ratio,
        std::vector<double>       &p_restitution_coefficient,
        std::vector<double>       &p_friction_coefficient,
        std::vector<double>       &p_rolling_viscous_damping_coefficient,
        std::vector<double>       &p_rolling_friction_coefficient,
        std::vector<double>       &p_surface_energy,
        std::vector<double>       &p_hamaker_constant,
        std::vector<double>       &p_thermal_conductivity,
        std::vector<double>       &p_specific_heat,
        std::vector<double>       &p_microhardness,
        std::vector<double>       &p_surface_slope,
        std::vector<double>       &p_surface_roughness,
        std::vector<double>       &p_thermal_accommodation,
        std::vector<double>       &p_real_youngs_modulus) const;
    };

    /**
     * @brief Parameters controlling particle insertion in DEM simulations.
     *
     * This class stores the insertion method, geometry, initial conditions,
     * and scheduling parameters used to introduce particles into the domain.
     * Four insertion methods are supported: file, list, plane, and volume.
     *
     * @tparam dim Number of spatial dimensions.
     */
    template <int dim>
    class InsertionInfo
    {
    public:
      /**
       * @brief Method used to insert particles.
       */
      enum class InsertionMethod
      {
        /// Insertion from a file listing particles
        file,
        /// Insertion from a list in the prm file
        list,
        /// Insertion from a plane defined by a point and a normal
        plane,
        /// Insertion within a volume defined by a box
        volume
      } insertion_method; ///< Method used to insert particles

      /// Number of particles inserted at each insertion step.
      int inserted_this_step;

      /// Frequency of insertion (in time steps).
      int insertion_frequency;

      /// Enable removal of particles in a specified region.
      bool removing_particles_in_region;

      /// First corner of the particle removal box.
      Point<3> clear_box_point_1;

      /// Second corner of the particle removal box.
      Point<3> clear_box_point_2;

      /// List of input files for the file insertion method.
      std::vector<std::string> list_of_input_files;

      /// Normal vector of the insertion plane (plane method).
      Tensor<1, 3> insertion_plane_normal_vector;

      /// Point on the insertion plane (plane method).
      Point<3> insertion_plane_point;

      // Position and velocity components for the list insertion method.
      std::vector<double> list_x, ///< x-position for list insertion.
        list_y,                   ///< y-position for list insertion.
        list_z,                   ///< z-position for list insertion.
        list_vx,                  ///< x-velocity for list insertion.
        list_vy,                  ///< y-velocity for list insertion.
        list_vz,                  ///< z-velocity for list insertion.
        list_wx,                  ///< x-angular velocity for list insertion.
        list_wy,                  ///< y-angular velocity for list insertion.
        list_wz,                  ///< z-angular velocity for list insertion.
        list_d,                   ///< Diameter for list insertion.
        list_T;                   ///< Temperature for list insertion

      /// Direction sequence for particle insertion (1st, 2nd, 3rd).
      std::vector<unsigned int> direction_sequence;

      /// First corner of the insertion box (volume method).
      Point<3> insertion_box_point_1;

      /// Second corner of the insertion box (volume method).
      Point<3> insertion_box_point_2;

      /// Initial translational velocity of inserted particles.
      Tensor<1, 3> initial_vel;

      /// Initial angular velocity of inserted particles.
      Tensor<1, 3> initial_omega;

      /// Function returning the initial temperature of a particle based on
      /// time or its position.
      std::shared_ptr<Function<dim>> initial_temperature_function;

      /// Minimum distance threshold between inserted particles.
      double distance_threshold;

      /// Maximum random offset applied to insertion positions.
      double insertion_maximum_offset;

      /// Random seed for particle insertion.
      int seed_for_insertion;

      /**
       * @brief Declare the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      static void
      declare_parameters(ParameterHandler &prm);

      /**
       * @brief Parse the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      parse_parameters(ParameterHandler &prm);
    };

    /**
     * @brief Numerical model parameters for DEM simulations.
     *
     * This class stores the contact detection, load balancing, force model
     * selection, time integration, and sparse contact optimization parameters.
     *
     * @tparam dim Number of spatial dimensions.
     */
    template <int dim>
    class ModelParameters
    {
    public:
      /**
       * @brief Load balancing strategy for parallel DEM simulations.
       */
      enum class LoadBalanceMethod
      {
        /// No load-balancing will be performed.
        none,
        /// Perform load-balancing once.
        once,
        /// Perform load-balancing at a given frequency.
        frequent,
        /// Perform load-balancing when the computational load among processes
        /// becomes too uneven.
        dynamic,
        /// Perform load-balancing in a similar manner to
        /// LoadBalanceMethod::dynamic but considering also the mobility status
        /// of the cells.
        dynamic_with_sparse_contacts
      } load_balance_method; ///< Load balancing strategy for parallel DEM
                             ///< simulations.

      /// Load balance step (for single-step load balancing).
      unsigned int load_balance_step;

      /// Load balance frequency (in time steps).
      unsigned int load_balance_frequency;

      /// Load balance threshold (for dynamic load balancing).
      double load_balance_threshold;

      /// Check frequency for dynamic load balancing.
      unsigned int dynamic_load_balance_check_frequency;

      /// Frequency of particle-particle and particle-wall contact detection.
      unsigned int contact_detection_frequency;

      /// Function returning the weight of a cell based on its barycenter
      /// position.
      std::shared_ptr<Function<dim>> cell_weight_function;

      /// Particle weight used for load balancing.
      unsigned int load_balance_particle_weight;

      /// Factor applied to particle weight for active cells in load balancing
      /// (only used with adaptive sparse contacts).
      double active_load_balancing_factor;

      /// Factor applied to particle weight for inactive cells in load balancing
      /// (only used with adaptive sparse contacts).
      double inactive_load_balancing_factor;

      /// Safety factor for dynamic contact search.
      double dynamic_contact_search_factor;

      /**
       * @brief Contact detection method used in the simulation.
       */
      enum class ContactDetectionMethod
      {
        /// Carry-out contact detection at a constant frequency.
        constant,
        /// Carry-out contact detection when the maximum displacement of a
        /// particle exceeds the smallest contact search criterion.
        dynamic
      } contact_detection_method; ///< Contact detection method used in the
                                  ///< simulation.

      /// Contact search neighborhood threshold (neighborhood diameter to
      /// particle diameter ratio).
      double neighborhood_threshold;

      /// Cut-off threshold beyond which Van der Waals forces are ignored.
      double dmt_cut_off_threshold;

      /// Particle-particle contact force model.
      ParticleParticleContactForceModel particle_particle_contact_force_model;

      /// Particle-wall contact force model.
      ParticleWallContactForceModel particle_wall_contact_force_method;

      /// Rolling resistance torque method.
      RollingResistanceMethod rolling_resistance_method;

      /// Model parameter for the EPSD rolling resistance model.
      double f_coefficient_epsd;

      /**
       * @brief Time integration method for particle motion.
       */
      enum class IntegrationMethod
      {
        /// Velocity Verlet second-order time integration scheme.
        velocity_verlet,
        /// Explicit Euler first-order time integration scheme.
        explicit_euler
      } integration_method; ///< Time integration method for particle motion.

      /// Solver type (DEM, CFD-DEM, or DEM multiphysics).
      DEM::SolverType solver_type;

      /// Enable sparse particle contacts to optimize performance.
      bool sparse_particle_contacts;

      /// Enable advection of particles using cell-averaged fluid velocity and
      /// acceleration.
      bool advect_particles;

      /// Minimum granular temperature for cells where particle contacts are
      /// evaluated.
      double granular_temperature_threshold;

      /// Maximum solid fraction for cells where particle contacts are always
      /// evaluated regardless of granular temperature.
      double solid_fraction_threshold;

      /// Disable position integration for particles.
      bool disable_position_integration;

      /**
       * @brief Declare the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      static void
      declare_parameters(ParameterHandler &prm);

      /**
       * @brief Parse the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
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
      /// Enable force and torque post-processing on wall boundaries.
      bool calculate_force_torque;

      /// Verbosity level for force and torque output.
      Parameters::Verbosity force_torque_verbosity;

      /// Output frequency (in time steps).
      unsigned int output_frequency;

      /// File name prefix for force and torque output.
      std::string force_torque_output_name;

      /// Center of mass used for torque computation.
      Point<3> point_center_mass;

      /**
       * @brief Declare the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      static void
      declare_parameters(ParameterHandler &prm);

      /**
       * @brief Parse the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      parse_parameters(ParameterHandler &prm);
    };

    /**
     * @brief Parameters for floating walls in DEM simulations.
     *
     * Floating walls are planar boundaries that can be activated and
     * deactivated at specified times during the simulation. Each wall
     * is defined by a point on its surface and a normal vector.
     *
     * @tparam dim Number of spatial dimensions.
     */
    template <int dim>
    class FloatingWalls
    {
    public:
      /// Number of floating walls.
      unsigned int floating_walls_number;

      /// A point on each floating wall surface.
      std::vector<Point<dim>> points_on_walls;

      /// Outward normal vector of each floating wall.
      std::vector<Tensor<1, dim>> floating_walls_normal_vectors;

      /// Activation time of each floating wall.
      std::vector<double> time_start;

      /// Deactivation time of each floating wall.
      std::vector<double> time_end;

      /**
       * @brief Declare the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      declare_parameters(ParameterHandler &prm) const;

      /**
       * @brief Parse the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      parse_parameters(ParameterHandler &prm);

      /**
       * @brief Declare default parameter entries for a single floating wall.
       *
       * @param[in,out] prm The parameter handler.
       */
      static void
      declareDefaultEntry(ParameterHandler &prm);

      /**
       * @brief Parse parameters for a single floating wall.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      parse_floating_wall(ParameterHandler &prm);

    private:
      unsigned int max_number_floating_walls = 9;
    };


    /**
     * @brief Boundary conditions for DEM simulations.
     *
     * This structure stores the boundary types, motion parameters, and
     * periodic boundary information for each boundary of the DEM domain.
     */
    struct BCDEM
    {
    public:
      /// Number of DEM boundary conditions.
      unsigned int DEM_BC_number;

      /**
       * @brief Type of boundary condition applied to a DEM domain boundary.
       */
      enum class BoundaryType
      {
        /// Static wall
        fixed_wall,
        /// Open boundary where particles may exit (outlet)
        outlet,
        /// Translating boundary at the velocity
        /// BCDEM::boundary_translational_velocity.
        translational,
        /// Rotating boundary described with BCDEM::boundary_rotational_vector,
        /// BCDEM::point_on_rotation_axis, and BCDEM::boundary_rotational_speed.
        rotational,
        /// Periodic boundary in the specified BCDEM::periodic_direction
        periodic
      };

      /// Boundary condition type for each boundary.
      std::vector<BoundaryType> bc_types;

      /// Boundary IDs designated as outlets.
      std::vector<unsigned int> outlet_boundaries;

      /// Translational velocity of each moving boundary.
      std::unordered_map<unsigned int, Tensor<1, 3>>
        boundary_translational_velocity;

      /// Rotational speed of each rotating boundary (rad/s).
      std::unordered_map<unsigned int, double> boundary_rotational_speed;

      /// Rotational axis vector of each rotating boundary.
      std::unordered_map<unsigned int, Tensor<1, 3>> boundary_rotational_vector;

      /// Point on the rotational axis of each rotating boundary.
      std::unordered_map<unsigned int, Point<3>> point_on_rotation_axis;

      /// First periodic boundary ID.
      types::boundary_id periodic_boundary_0;

      /// Second periodic boundary ID.
      types::boundary_id periodic_boundary_1;

      /// Direction of periodicity.
      types::boundary_id periodic_direction;

      /**
       * @brief Declare the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      declare_parameters(ParameterHandler &prm) const;

      /**
       * @brief Parse the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      parse_parameters(ParameterHandler &prm);

      /**
       * @brief Declare default parameter entries for a single boundary
       * condition.
       *
       * @param[in,out] prm The parameter handler.
       */
      static void
      declareDefaultEntry(ParameterHandler &prm);

      /**
       * @brief Parse boundary condition parameters for a single boundary.
       *
       * @param[in] prm The parameter handler.
       */
      void
      parse_boundary_conditions(const ParameterHandler &prm);

    private:
      unsigned int DEM_BC_number_max = 10;

      /**
       * @brief Initialize the containers used to store boundary condition
       * parameters.
       *
       * @param[in,out] boundary_trans_velocity Translational velocities.
       * @param[in,out] boundary_rot_speed Rotational speeds.
       * @param[in,out] boundary_rot_vector Rotational axis vectors.
       * @param[in,out] point_on_rot_axis Points on rotation axes.
       * @param[in,out] outlet_boundaries_id Outlet boundary IDs.
       * @param[in,out] boundaries_types Boundary types.
       */
      void
      initialize_containers(
        std::unordered_map<unsigned int, Tensor<1, 3>> &boundary_trans_velocity,
        std::unordered_map<unsigned int, double>       &boundary_rot_speed,
        std::unordered_map<unsigned int, Tensor<1, 3>> &boundary_rot_vector,
        std::unordered_map<unsigned int, Point<3>>     &point_on_rot_axis,
        std::vector<unsigned int>                      &outlet_boundaries_id,
        std::vector<BoundaryType>                      &boundaries_types) const;
    };

    /**
     * @brief Parameters for grid (mesh) motion in DEM simulations.
     *
     * This class defines the type and parameters of the grid motion, which
     * can be translational, rotational, or none.
     *
     * @tparam dim Number of spatial dimensions.
     */
    template <int dim>
    class GridMotion
    {
    public:
      /**
       * @brief Type of grid motion.
       */
      enum class MotionType
      {
        /// Translating grid at the velocity
        /// GridMotion::grid_translational_velocity.
        translational,
        /// Rotating grid at the speed GridMotion::grid_rotational_speed around
        /// the GridMotion::grid_rotational_axis.
        rotational,
        /// Static grid.
        none
      } motion_type; ///< Type of grid motion.

      /// Translational velocity of the moving grid.
      Tensor<1, dim> grid_translational_velocity;

      /// Rotational speed of the rotating grid (rad/s).
      double grid_rotational_speed;

      /// Rotational axis of the rotating grid (0=x, 1=y, 2=z).
      unsigned int grid_rotational_axis;

      /**
       * @brief Declare the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      static void
      declare_parameters(ParameterHandler &prm);

      /**
       * @brief Parse the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
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
      /// Enable Lagrangian post-processing.
      bool lagrangian_post_processing_enabled;

      /// Enable force chains visualization.
      bool force_chains;

      /// Enable logging of particle-wall contact statistics.
      bool particle_wall_collision_statistics;

      /// Log collisions with all walls (if false, only selected boundaries).
      bool log_collisions_with_all_walls;

      /// Boundary IDs of the walls for which collisions are logged.
      std::vector<int> particle_wall_collision_boundary_ids;

      /// Verbosity level for collision statistics output.
      Parameters::Verbosity collision_verbosity;

      /// File name for exporting collision statistics (CSV format).
      std::string collision_stats_file_name;

      /**
       * @brief Declare the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      static void
      declare_parameters(ParameterHandler &prm);

      /**
       * @brief Parse the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      parse_parameters(ParameterHandler &prm);
    };

    /**
     * @brief Parameters for a floating grid in DEM simulations.
     *
     * A floating grid is an auxiliary mesh that can move independently of
     * the main triangulation, activated between specified start and end times.
     *
     * @tparam dim Number of spatial dimensions.
     */
    template <int dim>
    class FloatingGrid
    {
    public:
      /// Mesh parameters for the floating grid.
      Parameters::Mesh mesh;

      /// Motion parameters for the floating grid.
      Parameters::Lagrangian::GridMotion<dim> motion;

      /// Activation time of the floating grid.
      double time_start;

      /// Deactivation time of the floating grid.
      double time_end;

      /**
       * @brief Declare the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      declare_parameters(ParameterHandler &prm);

      /**
       * @brief Parse the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
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
      /// Location of the first photon to be inserted.
      Point<3> starting_point;

      /// Unit direction vectors along which photons are inserted relative to
      /// the starting point.
      std::vector<Tensor<1, 3>> insertion_directions_units_vector;

      /// Number of photons to insert along each direction.
      std::vector<unsigned int> n_photons_each_directions;

      /// Spacing between consecutive photons along each direction (at zero
      /// offset).
      std::vector<double> step_between_photons_each_directions;

      /// Reference unit tensor defining the photon displacement direction.
      Tensor<1, 3> ref_displacement_tensor_unit;

      /// Maximum random offset applied to the photon insertion position.
      double max_insertion_offset;

      /// Random seed for photon insertion position offset.
      unsigned int prn_seed_photon_insertion;

      /// Maximum angular offset applied to the photon displacement direction.
      double max_angular_offset;

      /// Random seed for photon displacement angular offset.
      unsigned int prn_seed_photon_displacement;

      /**
       * @brief Declare the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      static void
      declare_parameters(ParameterHandler &prm);

      /**
       * @brief Parse the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      parse_parameters(ParameterHandler &prm);
    };

  } // namespace Lagrangian
} // namespace Parameters
#endif /* lethe_parameters_lagrangian_h */
