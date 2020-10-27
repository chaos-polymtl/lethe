/*
 * parametersdem.h
 *
 *  Created on: Dec 16, 2019
 *      Author: shahab
 */
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <string>

#ifndef PARAMETERS_LAGRANGIAN_H_
#  define PARAMETERS_LAGRANGIAN_H_

using namespace dealii;
namespace Parameters
{
  namespace Lagrangian
  {
    struct PhysicalProperties
    {
      // Gravitational acceleration
      double gx, gy, gz;

      // Particle diameter and density
      double diameter;
      double density;

      // Young's modulus of particle and wall
      double youngs_modulus_particle;
      double youngs_modulus_wall;

      // Poisson's ratios of particle and wall
      double poisson_ratio_particle;
      double poisson_ratio_wall;

      // Coefficients of restituion of particle and wall
      double restitution_coefficient_particle;
      double restitution_coefficient_wall;

      // Friction coefficients of particle and wall
      double friction_coefficient_particle;
      double friction_coefficient_wall;

      // Rollinrg friction coefficients of particle and wall
      double rolling_friction_particle;
      double rolling_friction_wall;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

    struct InsertionInfo
    {
      // Choosing insertion method
      enum class InsertionMethod
      {
        uniform,
        non_uniform
      } insertion_method;

      // Total number of particles
      unsigned int total_particle_number;

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

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

    struct ModelParameters
    {
      // Particle-particle, particle-wall broad and fine search frequency
      unsigned int repartition_frequency;

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

      // Choosing integration method
      enum class IntegrationMethod
      {
        velocity_verlet,
        explicit_euler
      } integration_method;

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
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
      void
      declareDefaultEntry(ParameterHandler &prm);
      void
      parse_floating_wall(ParameterHandler &prm);
    };

    template <int dim>
    class BoundaryMotions
    {
    public:
      // Number of moving boundaries
      unsigned int moving_boundary_number;

      // Translational velocities of moving boundaries
      std::unordered_map<int, Tensor<1, dim>> boundary_translational_velocity;

      // Rotational speeds of rotating boundaries
      std::unordered_map<int, double> boundary_rotational_speed;

      void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
      void
      declareDefaultEntry(ParameterHandler &prm);
      void
      parse_boundary_motions(ParameterHandler &  prm,
                             const unsigned int &boundary_id);
    };

  } // namespace Lagrangian
} // namespace Parameters

#endif /* PARAMETERS_H_ */
