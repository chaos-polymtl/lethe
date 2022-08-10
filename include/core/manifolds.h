/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019 -
 */

#ifndef lethe_manifolds_h
#define lethe_manifolds_h

#include <deal.II/base/parameter_handler.h>

#include <deal.II/numerics/data_postprocessor.h>



using namespace dealii;

namespace Parameters
{
  /**
   * @brief Manifolds - This class manages attaching manifolds to
   * faces when a manual attachment is required (for example if the mesh
   * is coming from a GMSH file)
   */
  class Manifolds
  {
  public:
    enum class ManifoldType
    {
      none,
      spherical,
      iges
    };

    // ID of boundary condition
    std::vector<unsigned int> id;

    // List of boundary type for each number
    std::vector<ManifoldType> types;

    // Arguments of manifold
    std::vector<double> arg1;
    std::vector<double> arg2;
    std::vector<double> arg3;
    std::vector<double> arg4;
    std::vector<double> arg5;
    std::vector<double> arg6;

    // File names for cad manifolds
    std::vector<std::string> cad_files;

    // Number of boundary conditions
    unsigned int size;
    unsigned int max_size;

    void
    parse_boundary(ParameterHandler &prm, unsigned int i_bc);

    void
    declareDefaultEntry(ParameterHandler &prm, unsigned int i_bc);
    void
    declare_parameters(ParameterHandler &prm);

    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace Parameters

/**
 * @brief BounearyPostprocessor Post-processor class used to attach the boundary
 * id to the faces when outputting the surfaces within the domain.
 */
template <int dim>
class BoundaryPostprocessor : public DataPostprocessorScalar<dim>
{
public:
  BoundaryPostprocessor()
    : DataPostprocessorScalar<dim>("boundary_id", update_quadrature_points)
  {}
  virtual void
  evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &input_data,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    auto current_cell = input_data.template get_cell<dim>();

    for (unsigned int p = 0; p < input_data.evaluation_points.size(); ++p)
      {
        unsigned int boundary_id  = 0;
        double       min_distance = DBL_MAX;

        for (const auto face : current_cell->face_indices())
          {
            if (current_cell->face(face)->at_boundary())
              {
                for (const auto vertex :
                     current_cell->face(face)->vertex_indices())
                  {
                    double distance = input_data.evaluation_points[p].distance(
                      current_cell->face(face)->vertex(vertex));
                    if (distance < min_distance)
                      {
                        min_distance = distance;

                        boundary_id = current_cell->face(face)->boundary_id();
                      }
                  }
              }
          }
        computed_quantities[p][0] = boundary_id;
      }
  };

  virtual void
  evaluate_scalar_field(
    const DataPostprocessorInputs::Scalar<dim> &input_data,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    auto current_cell = input_data.template get_cell<dim>();

    for (unsigned int p = 0; p < input_data.evaluation_points.size(); ++p)
      {
        unsigned int boundary_id  = 0;
        double       min_distance = DBL_MAX;

        for (const auto face : current_cell->face_indices())
          {
            if (current_cell->face(face)->at_boundary())
              {
                for (const auto vertex :
                     current_cell->face(face)->vertex_indices())
                  {
                    double distance = input_data.evaluation_points[p].distance(
                      current_cell->face(face)->vertex(vertex));
                    if (distance < min_distance)
                      {
                        min_distance = distance;

                        boundary_id = current_cell->face(face)->boundary_id();
                      }
                  }
              }
          }
        computed_quantities[p][0] = boundary_id;
      }
  }
};


/**
 * @brief Attaches manifold to boundaries of the triangulation
 *
 * @param triangulation The triangulation to manifolds are attached
 *
 * @param manifolds_parameters The information about the type of manifolds attached to the faces
 */
template <int dim, int spacedim = dim>
void
attach_manifolds_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim, spacedim>>
                              triangulation,
  const Parameters::Manifolds manifolds);

/**
 * @brief Attaches CAD manifolds using IGES files to boundaries of the triangulation
 * This function attaches a CAD manifold to the faces of a triangulation. Note
 * that this only works in 3D.
 *
 * @param triangulation The triangulation to manifolds are attached
 *
 * @param std::string Filename of the cad file
 *
 * @param manifold_id Identifier of the manifold
 */
void attach_cad_to_manifold(
  std::shared_ptr<parallel::DistributedTriangulationBase<2>> triangulation,
  std::string                                                cad_name,
  unsigned int                                               manifold_id);

void attach_cad_to_manifold(
  std::shared_ptr<parallel::DistributedTriangulationBase<2, 3>> triangulation,
  std::string                                                   cad_name,
  unsigned int                                                  manifold_id);

void attach_cad_to_manifold(
  std::shared_ptr<parallel::DistributedTriangulationBase<3>> triangulation,
  std::string                                                cad_name,
  unsigned int                                               manifold_id);



#endif
