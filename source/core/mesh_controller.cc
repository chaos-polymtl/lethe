#include <core/mesh_controller.h>

double
MeshController::calculate_coarsening_factor(
  unsigned int current_number_of_elements)
{
  // Parameters of the controller
  double P = 0.5;
  double I = 0.3;
  double D = 0.1;

  // Evaluation of the error used to control the mesh refinement.
  double error = static_cast<signed>(target_number_of_elements -
                                     current_number_of_elements) /
                 static_cast<double>(target_number_of_elements);
  double previous_error = static_cast<signed>(target_number_of_elements -
                                              previous_number_of_elements) /
                          static_cast<double>(target_number_of_elements);

  previous_number_of_elements = current_number_of_elements;
  previous_mesh_control_error = previous_mesh_control_error + error;

  // Calculate the coarsening factor
  double coarsening_fraction_controled =
    -error * P - previous_mesh_control_error * I - (error - previous_error) * D;

  // Saturate the coarsening parameter.
  if (coarsening_fraction_controled < 0.0)
    {
      // Stop the integration if we are at the saturation point.
      coarsening_fraction_controled = 0;
      if (error > 0.0)
        previous_mesh_control_error = previous_mesh_control_error - error;
    }
  if (coarsening_fraction_controled > 1.0)
    {
      coarsening_fraction_controled = 1;
      // Stop the integration if we are at the saturation point.
      if (error < 1.0)
        previous_mesh_control_error = previous_mesh_control_error - error;
    }

  return coarsening_fraction_controled;
}
