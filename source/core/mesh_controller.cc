//
// Created by luckab on 12/07/23.
//

#include <core/mesh_controller.h>


using namespace std;
double
MeshController::calculate_coarsening_factor(int current_number_of_cell){
  // Parameter of the controller
  double P=0.5;
  double I=0.3;
  double D=0.1;

  // Evaluation of the error use to control the mesh refinement.
  double error=(target_number_of_elements-current_number_of_cell)/(double)target_number_of_elements;
  double previous_error=(target_number_of_elements-previous_number_of_elements)/(double)target_number_of_elements;
  previous_number_of_elements=current_number_of_cell;
  previous_mesh_control_error=previous_mesh_control_error+error;

  // Calculate the coarsening factor
  double coarsening_fraction_controled = -error * P -
                                      previous_mesh_control_error * I -(error-previous_error)*D ;

  // Saturate the coarsening parameter.
  if(coarsening_fraction_controled<0.0)
    {
      // Stop the integration if we are at the saturation point.
      coarsening_fraction_controled=0;
      if (error>0.0){
          previous_mesh_control_error=previous_mesh_control_error-error;
        }
    }
  if(coarsening_fraction_controled>1.0){
      coarsening_fraction_controled=1;
      // Stop the integration if we are at the saturation point.
      if (error<1.0){
          previous_mesh_control_error=previous_mesh_control_error-error;
        }
    }

  return coarsening_fraction_controled;
}


