#include <rpt/particle_detector_interactions.h>
#include <stdlib.h>
#include <time.h>

#include <cmath>

template <int dim>
void
ParticleDetectorInteractions<dim>::calculate_position_parameters()
{
  // Vector from the detector face center to the particle position (Xdp)
  face_detector_particle_distance = particle_position - detector_face_position;

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
  Tensor<1, dim> Z       = h * detector_orientation_z;
  Tensor<1, dim> X       = Z - face_detector_particle_distance;
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
  double theta_max, theta_min, theta_cri;
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

      // Calculate relevant distances prior other calculations
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
      else // Case 2 : S2 (h < 0)
        {
          theta_max = M_PI_2 + std::atan(std::abs(h) / OB_distance);
          theta_cri = M_PI_2 + std::atan(std::abs(h) / OB_distance);
          theta_min = std::atan(OB_distance / (l - std::abs(h)));
        }

      // Calculate theta angle and its weigthting factor for the Monte Carlo
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
    }
}

template <int dim>
void
ParticleDetectorInteractions<dim>::calculate_detector_path_length()
{}

template <int dim>
void
ParticleDetectorInteractions<dim>::calculate_reactor_path_length()
{}

template <int dim>
void
ParticleDetectorInteractions<dim>::calculate_detector_interaction_probability()
{}

template <int dim>
void
ParticleDetectorInteractions<dim>::calculate_non_interaction_probability()
{}

template <int dim>
void
ParticleDetectorInteractions<dim>::calculate_efficiency()
{}

template <int dim>
double
ParticleDetectorInteractions<dim>::calculate_count()
{
  // Loop for Monte Carlo here

  // Generate random values for Monte Carlo
  srand(time(NULL));
  double n_alpha = (double)rand() / RAND_MAX;
  double n_theta = (double)rand() / RAND_MAX;


  calculate_position_parameters();
  calculate_solid_angle(n_alpha, n_theta);

  double dummy = 0;
  return dummy;
}


template class ParticleDetectorInteractions<3>;