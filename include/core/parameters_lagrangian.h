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

      // Young's modulus of particle and wall
      double Yp;
      double Yw;

      // Poisson's ratios of particle and wall
      double vp;
      double vw;

      // Coefficients of restituion of particle and wall
      double ep;
      double ew;

      // Friction coefficients of particle and wall
      double mup;
      double muw;

      // Rollinrg friction coefficients of particle and wall
      double murp;
      double murw;

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
      int insertFrequency;

      // Insertion box info (xmin,xmax,ymin,ymax,zmin,zmax)
      double x_min, y_min, z_min, x_max, y_max, z_max;

      // Insertion distance threshold
      double distance_threshold;

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

    struct SimulationModel
    {
      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };

  } // namespace Lagrangian
} // namespace Parameters

#endif /* PARAMETERS_H_ */
