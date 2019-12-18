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

#ifndef PARAMETERS_LAGRANGIAN_H_
#  define PARAMETERS_LAGRANGIAN_H_

using namespace dealii;

namespace Parameters
{
  namespace Lagrangian
  {
    struct SimulationControl
    {
      // Time step
      double dt;

      // End time step
      int tFinal;

      // Total number of particles
      int nTotal;

      // Write frequency
      int writeFrequency;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };


    struct PhysicalProperties
    {
      // Gravitational acceleration
      double gx, gy, gz;

      // Particle diameter and density
      double diameter;
      double density;

      // Spring and dashpot normal constants
      double kn;
      double ethan;

      // Spring and dashpot tangential constants
      double kt;
      double ethat;

      // Coefficient of friction
      double mu;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };


    struct InsertionInfo
    {
      // Insertion time step
      int tInsertion;

      // Inserted number of particles at each time step
      int nInsert;

      // Insertion frequency
      int insertFrequncy;

      // Insertion box info (xmin,xmax,ymin,ymax,zmin,zmax)
      double x_min, y_min, z_min, x_max, y_max, z_max;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };


    struct OutputProperties
    {
      // Number of properties
      int numProperties;

      // Number of fields
      int numFields;

      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

  } // namespace Lagrangian
} // namespace Parameters

#endif /* PARAMETERS_H_ */
