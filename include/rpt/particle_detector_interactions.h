
#ifndef lethe_particle_detector_interactions_h
#define lethe_particle_detector_interactions_h

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <rpt/detector.h>
#include <rpt/parameters_rpt.h>
#include <rpt/radioactive_particle.h>
#include <rpt/rpt_calculating_parameters.h>

template <int dim>
class ParticleDetectorInteractions
{
public:
  ParticleDetectorInteractions(RadioParticle<dim> &      particle,
                               Detector<dim> &           detector,
                               RPTCalculatingParameters &rpt_parameters)
    : particle_position(particle.get_position())
    , detector_face_position(detector.get_face_position())
    , detector_middle_position(detector.get_middle_position())
    , detector_radius(detector.get_radius())
    , detector_length(detector.get_length())
    , fixed_parameters(rpt_parameters.rpt_param)
    , initial_parameters(rpt_parameters.initial_param)
  {}

  double
  calculate_count();

  double
  get_h();

  double
  get_rho();

  double
  get_alpha(double n_alpha, double n_theta);

  double
  get_theta(double n_alpha, double n_theta);

  double
  get_detector_path_length(double n_alpha, double n_theta);

  double
  get_reactor_path_length(double n_alpha, double n_theta);

private:
  void
  calculate_position_parameters();

  void
  calculate_solid_angle(double n_alpha, double n_theta);

  double
  calculate_detector_path_length();

  double
  calculate_reactor_path_length();

  double
  calculate_detector_interaction_probability(double &detector_path_length);

  double
  calculate_non_interaction_probability(double &reactor_path_length);

  void
  calculate_efficiency();

  std::vector<double>
  solve_t(Tensor<2, dim> e_inverse,
          Tensor<1, dim> detector_particle_origin,
          Tensor<1, dim> particle_position_rotation);

  double         efficiency;
  double         alpha;
  double         theta;
  double         theta_cri;
  double         weighting_factor_alpha;
  double         weighting_factor_theta;
  double         h;
  double         rho;
  double         OA_distance;
  double         OB_distance;
  Tensor<1, dim> detector_orientation_x;
  Tensor<1, dim> detector_orientation_y;
  Tensor<1, dim> detector_orientation_z;

  Point<dim> particle_position;
  Point<dim> detector_face_position;
  Point<dim> detector_middle_position;

  double                           detector_radius;
  double                           detector_length;
  Parameters::RPTParameters        fixed_parameters;
  Parameters::InitialRPTParameters initial_parameters;
};



#endif // lethe_particle_detector_interactions_h
