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
#define PARAMETERS_LAGRANGIAN_H_

using namespace dealii;

namespace Parameters {
namespace Lagrangian {
struct SimulationControl {
  // Time step
  double dt;

  // End time step
  int final_time_step;

  // Total number of particles
  unsigned int total_particle_number;

  // Write frequency
  int write_frequency;

  static void declare_parameters(ParameterHandler &prm);
  void parse_parameters(ParameterHandler &prm);
};

struct PhysicalProperties {
  // Gravitational acceleration
  double gx, gy, gz;

  // Particle diameter and density
  double diameter;
  double density;

  // Young's modulus of particle and wall
  double Youngs_modulus_particle;
  double Youngs_modulus_wall;

  // Poisson's ratios of particle and wall
  double Poisson_ratio_particle;
  double Poisson_ratio_wall;

  // Coefficients of restituion of particle and wall
  double restitution_coefficient_particle;
  double restitution_coefficient_wall;

  // Friction coefficients of particle and wall
  double friction_coefficient_particle;
  double friction_coefficient_wall;

  // Rollinrg friction coefficients of particle and wall
  double rolling_friction_particle;
  double rolling_friction_wall;

  static void declare_parameters(ParameterHandler &prm);
  void parse_parameters(ParameterHandler &prm);
};

struct InsertionInfo {
  // Insertion time step
  int insertion_steps_number;

  // Inserted number of particles at each time step
  int inserted_this_step;

  // Insertion frequency
  int insertion_frequency;

  // Insertion box info (xmin,xmax,ymin,ymax,zmin,zmax)
  double x_min, y_min, z_min, x_max, y_max, z_max;

  // Insertion distance threshold
  double distance_threshold;

  // Insertion random number range
  double random_number_bin;

  static void declare_parameters(ParameterHandler &prm);
  void parse_parameters(ParameterHandler &prm);
};

struct OutputProperties {
  // Number of properties
  int properties_number;

  // Output directory
  std::string output_folder;

  // General information file (.pvtu) prefix
  std::string general_file_prefix;

  // Result (.vtu) name prefix
  std::string result_prefix;

  static void declare_parameters(ParameterHandler &prm);
  void parse_parameters(ParameterHandler &prm);
};

struct ModelParameters {
  // Particle-particle broad search frequency
  int pp_broad_search_frequency;

  // Particle-wall broad search frequency
  int pw_broad_search_frequency;

  // Print simulation info frequency
  int print_info_frequency;

  // Choosing particle-particle contact force model
  enum class PPContactForceModel {
    pp_linear,
    pp_nonlinear
  } pp_contact_force_method;

  // Choosing particle-wall contact force model
  enum class PWContactForceModel {
    pw_linear,
    pw_nonlinear
  } pw_contact_force_method;

  static void declare_parameters(ParameterHandler &prm);
  void parse_parameters(ParameterHandler &prm);
};

} // namespace Lagrangian
} // namespace Parameters

#endif /* PARAMETERS_H_ */
