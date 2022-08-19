#include <rpt/particle_detector_interactions.h>

#include <cmath>
#include <cstdlib>

template <int dim>
void
ParticleDetectorInteractions<dim>::calculate_position_parameters()
{
  // Vector from the detector face center to the particle position (Xdp)
  Tensor<1, dim> face_detector_particle_distance =
    particle_position - detector_face_position;

  // Detector orientation vector parallel to z' in particle-detector coordinate
  // (e'z)
  Tensor<1, dim> face_middle_vector =
    detector_face_position - detector_middle_position;
  detector_orientation_z = face_middle_vector / face_middle_vector.norm();

  // Projection of the distance vector on the detector axis which is parallel to
  // the z' axis Particle position in particle-detector coordinate frame (h =
  // z'p)
  h = scalar_product(face_detector_particle_distance, detector_orientation_z);

  // x' coordinate of the particle-detector frame (e'x) (see equations 29 to 31)
  Tensor<1, dim> Z = h * detector_orientation_z;
  Tensor<1, dim> X = Z - face_detector_particle_distance;

  if (X.norm() < 1e-6) // When particle is aligned to the detector
    {
      detector_orientation_x[0] = 0;
      detector_orientation_x[1] = 0;
      detector_orientation_x[2] = 1;
    }
  else
    detector_orientation_x = X / X.norm();

  // y' coordinate of the particle-detector frame (e'y)
  detector_orientation_y =
    cross_product_3d(detector_orientation_z, detector_orientation_x);

  // Rho parameter : distance between the center of the detector and a line
  // parallel to the detector axis containing the particle position
  rho = X.norm();
}

template <int dim>
void
ParticleDetectorInteractions<dim>::calculate_solid_angle(double n_alpha,
                                                         double n_theta)
{
  double alpha_max;
  double theta_max = 0, theta_min = 0;
  double r = detector_radius;
  double l = detector_length;

  // Case 1 or 2: S1 or S2 (Tracer views the detector from both top and/or
  // lateral side of the detector)
  if (rho > detector_radius)
    {
      // Calculate alpha angles and its weigthting factor for the Monte Carlo
      // iteration
      alpha_max              = std::asin(r / rho);
      alpha                  = alpha_max * (2 * n_alpha - 1);
      weighting_factor_alpha = alpha_max / M_PI;

      // Calculate relevant distances prior next and other calculations
      OB_distance =
        rho * std::cos(alpha) -
        std::sqrt(std::pow(r, 2) - std::pow(rho * std::sin(alpha), 2));
      OA_distance =
        rho * std::cos(alpha) +
        std::sqrt(std::pow(r, 2) - std::pow(rho * std::sin(alpha), 2));

      if (h > 0) // Case 1 : S1 (h > 0)
        {
          theta_max = std::atan(OA_distance / h);
          theta_cri = std::atan(OB_distance / h);
          theta_min = std::atan(OB_distance / (h + l));
        }
      else if (std::abs(h) < 1e-6) // Case 1 : S1 (h = 0)
        {
          theta_max = M_PI_2; // pi/2
          theta_cri = M_PI_2; // pi/2
          theta_min = std::atan(OB_distance / l);
        }
      else // Cases where h < 0 happens are not implemented yet since they are
           // not usual
        {
          AssertThrow(
            h >= 0,
            ExcMessage(
              "For now, it is not possible to calculate count for particle positions on the side of the detector."));
        }

      // Calculate theta angle and its weighting factor for the Monte Carlo
      // iteration
      theta                  = std::acos(std::cos(theta_min) -
                        n_theta * (std::cos(theta_min) - std::cos(theta_max)));
      weighting_factor_theta = (std::cos(theta_min) - std::cos(theta_max)) / 2;
    }
  // Case 3 : S3 (Tracer views the detector only from the top surface)
  else // (rho < detector_radius)
    {
      theta_max = std::atan((r + rho) / h);
      theta_cri = std::atan((r - rho) / h);
      theta_min = 0;

      // Calculate theta angle and its weigthting factor for the Monte Carlo
      // iteration
      theta                  = std::acos(std::cos(theta_min) -
                        n_theta * (std::cos(theta_min) - std::cos(theta_max)));
      weighting_factor_theta = (std::cos(theta_min) - std::cos(theta_max)) / 2;

      if (theta <= theta_cri)
        {
          alpha_max              = M_PI;
          weighting_factor_alpha = 1;
        }
      else // theta > theta_cri
        {
          alpha_max =
            std::acos((std::pow(rho, 2) + std::pow(h * std::tan(theta), 2) -
                       std::pow(r, 2)) /
                      (2 * h * rho * std::tan(theta)));
          weighting_factor_alpha = alpha_max / M_PI;
        }
      alpha = alpha_max * (2 * n_alpha - 1);

      // Calculate relevant distances prior next and other calculations
      OB_distance =
        rho * std::cos(alpha) -
        std::sqrt(std::pow(r, 2) - std::pow(rho * std::sin(alpha), 2));
      OA_distance =
        rho * std::cos(alpha) +
        std::sqrt(std::pow(r, 2) - std::pow(rho * std::sin(alpha), 2));
    }
}

template <int dim>
double
ParticleDetectorInteractions<dim>::calculate_detector_path_length()
{
  double detector_path_length;
  double l       = detector_length;
  double theta_1 = std::atan(OA_distance / (h + l));

  // Case 1 : S1  (Tracer views the detector from both top and lateral side of
  // the detector)
  if (rho > detector_radius)
    {
      if (theta <= theta_cri && theta <= theta_1)
        detector_path_length =
          (h + l) / std::cos(theta) - OB_distance / std::sin(theta);
      else if (theta <= theta_cri && theta > theta_1)
        detector_path_length = (OA_distance - OB_distance) / std::sin(theta);
      else if (theta > theta_cri && theta <= theta_1)
        detector_path_length = l / std::cos(theta);
      else // (theta > theta_cri && theta > theta_1)
        detector_path_length =
          OA_distance / std::sin(theta) - h / std::cos(theta);
    }
  // Case 3 : S3 (Tracer views the detector only from the top surface)
  else
    {
      double theta_2 = std::atan(OA_distance / h);

      if (theta <= theta_1)
        detector_path_length = l / std::cos(theta);
      else if (theta > theta_1 && theta <= theta_2)
        detector_path_length =
          OA_distance / std::sin(theta) - h / std::cos(theta);
      else // Do not enters in the detector
        detector_path_length = 0;
    }

  return detector_path_length;
}

template <int dim>
std::vector<double> ParticleDetectorInteractions<dim>::solve_t(
  Tensor<2, dim> e_inverse,
  Tensor<1, dim> detector_particle_origin,
  Tensor<1, dim> particle_position_rotation)
{
  // Function value when evaluate with t (circle equation)
  auto F = [=](double t) {
    double R = parameters.reactor_radius;

    Tensor<1, dim> line_equations(
      {particle_position_rotation[0] + t * std::sin(theta) * std::cos(alpha),
       particle_position_rotation[1] + t * std::sin(theta) * std::sin(alpha),
       particle_position_rotation[2] + t * std::cos(M_PI - theta)});

    return std::pow((e_inverse * (line_equations) +
                     detector_particle_origin)[0],
                    2) +
           std::pow((e_inverse * (line_equations) +
                     detector_particle_origin)[1],
                    2) -
           std::pow(R, 2);
  };

  // Numerical derivative of the function value when evaluate with t
  auto dF = [=](double t) {
    double dt = 0.0001;

    return (F(t + dt) - F(t - dt)) / (2 * dt);
  };

  // Values for solving (positive/negative)
  std::vector<double> t0 = {-1, 1};
  double              tolerance, dt0, max_iteration;
  tolerance     = 1e-6;
  max_iteration = 1000;
  dt0           = 1;

  // Newton method
  unsigned int        i;
  double              t, dtn;
  std::vector<double> t_solutions(t0.size(), 0);

  for (unsigned int n = 0; n < t0.size(); n++)
    {
      i   = 0;
      dtn = dt0;
      t   = t0[n];

      while (std::abs(dtn) > tolerance && i < max_iteration)
        {
          dtn = -F(t) / dF(t);
          t += dtn;
          i++;

          if (i == max_iteration)
            std::cout << "Solving t didn't work" << std::endl;
        }

      t_solutions[n] = t;
    }

  return t_solutions;
}

template <int dim>
double
ParticleDetectorInteractions<dim>::calculate_reactor_path_length()
{
  Tensor<1, dim> h_vector, detector_particle_origin,
    particle_position_translation, particle_position_rotation;
  Point<dim> intersection_position;

  // Calculate translated particle position (eq 59 - 63)
  h_vector                      = h * detector_orientation_z;
  detector_particle_origin      = particle_position - h_vector;
  particle_position_translation = particle_position - detector_particle_origin;

  // Calculate rotated particle position (eq 64)
  Tensor<1, dim, Tensor<1, dim>> detector_orientation_vector(
    {detector_orientation_x, detector_orientation_y, detector_orientation_z});
  Tensor<2, dim> detector_orientation_matrix(detector_orientation_vector);
  particle_position_rotation =
    detector_orientation_matrix * particle_position_translation;

  // Find t value of the parametric equation for
  std::vector<double> t = solve_t(invert(detector_orientation_matrix),
                                  detector_particle_origin,
                                  particle_position_rotation);

  // Evaluate and determinate the intersection point closer to the detector
  std::vector<Point<dim>>     intersection_point(t.size());
  std::vector<Tensor<1, dim>> intersection_detector_distance_vector(t.size());
  std::vector<double>         intersection_detector_distance(t.size());

  for (unsigned int i = 0; i < t.size(); i++)
    {
      Tensor<1, dim> line_equations(
        {particle_position_rotation[0] +
           t[i] * std::sin(theta) * std::cos(alpha),
         particle_position_rotation[1] +
           t[i] * std::sin(theta) * std::sin(alpha),
         particle_position_rotation[2] + t[i] * std::cos(M_PI - theta)});

      intersection_point[i] =
        detector_particle_origin +
        invert(detector_orientation_matrix) * line_equations;

      intersection_detector_distance_vector[i] =
        detector_face_position - intersection_point[i];
      intersection_detector_distance[i] =
        std::abs(intersection_detector_distance_vector[i].norm());
    }

  double reactor_path_length =
    (intersection_detector_distance[0] < intersection_detector_distance[1]) ?
      std::abs((particle_position - intersection_point[0]).norm()) :
      std::abs((particle_position - intersection_point[1]).norm());

  return reactor_path_length;
}

template <int dim>
double
ParticleDetectorInteractions<dim>::calculate_detector_interaction_probability(
  double &detector_path_length)
{
  double mu_d = parameters.attenuation_coefficient_detector;
  double detector_interaction_probability =
    1 - std::exp(-mu_d * detector_path_length);

  return detector_interaction_probability;
}

template <int dim>
double
ParticleDetectorInteractions<dim>::calculate_non_interaction_probability(
  double &reactor_path_length)
{
  double mu_a                        = attenuation_coefficient_reactor;
  double non_interaction_probability = std::exp(-mu_a * reactor_path_length);

  return non_interaction_probability;
}

template <int dim>
void
ParticleDetectorInteractions<dim>::calculate_efficiency()
{
  // Initialize efficiency
  efficiency = 0;

  // Calculate position parameters for the particle position
  calculate_position_parameters();

  const double denum = 1. / RAND_MAX;

  for (unsigned int i = 0; i < parameters.n_monte_carlo_iteration; i++)
    {
      // Generate random values for Monte Carlo
      double n_alpha = (double)rand() * denum;
      double n_theta = (double)rand() * denum;

      // Initialize variables
      double detector_path_length, reactor_path_length,
        detector_interaction_probability, non_interaction_probability;

      // Calculate all parameters and distances for efficiency
      calculate_solid_angle(n_alpha, n_theta);
      detector_path_length = calculate_detector_path_length();
      detector_interaction_probability =
        calculate_detector_interaction_probability(detector_path_length);
      reactor_path_length = calculate_reactor_path_length();
      non_interaction_probability =
        calculate_non_interaction_probability(reactor_path_length);

      efficiency += weighting_factor_alpha * weighting_factor_theta *
                    detector_interaction_probability *
                    non_interaction_probability;
    }

  efficiency /= parameters.n_monte_carlo_iteration;
}

template <int dim>
double
ParticleDetectorInteractions<dim>::calculate_count()
{
  double T, nu, R, phi, tau, count;

  T   = parameters.sampling_time;
  nu  = parameters.gamma_rays_emitted;
  R   = activity;
  phi = parameters.peak_to_total_ratio;
  tau = dead_time;

  calculate_efficiency();

  // Calculate count for a particle position and one detector
  count =
    (T * nu * R * phi * efficiency) / (1 + tau * nu * R * phi * efficiency);

  return count;
}

template <int dim>
double
ParticleDetectorInteractions<dim>::get_h()
{
  calculate_position_parameters();

  return h;
}

template <int dim>
double
ParticleDetectorInteractions<dim>::get_rho()
{
  calculate_position_parameters();

  return rho;
}

template <int dim>
double
ParticleDetectorInteractions<dim>::get_alpha(double n_alpha, double n_theta)
{
  calculate_position_parameters();
  calculate_solid_angle(n_alpha, n_theta);

  return alpha;
}

template <int dim>
double
ParticleDetectorInteractions<dim>::get_theta(double n_alpha, double n_theta)
{
  calculate_position_parameters();
  calculate_solid_angle(n_alpha, n_theta);

  return theta;
}

template <int dim>
double
ParticleDetectorInteractions<dim>::get_detector_path_length(double n_alpha,
                                                            double n_theta)
{
  calculate_position_parameters();
  calculate_solid_angle(n_alpha, n_theta);
  double detector_path_length = calculate_detector_path_length();

  return detector_path_length;
}

template <int dim>
double
ParticleDetectorInteractions<dim>::get_reactor_path_length(double n_alpha,
                                                           double n_theta)
{
  calculate_position_parameters();
  calculate_solid_angle(n_alpha, n_theta);
  double reactor_path_length = calculate_reactor_path_length();

  return reactor_path_length;
}



template class ParticleDetectorInteractions<3>;
