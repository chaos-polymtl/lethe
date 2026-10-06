// SPDX-FileCopyrightText: Copyright (c) 2019-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/manifolds.h>

#include <deal.II/base/point.h>

#include <deal.II/grid/manifold_lib.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

namespace Parameters
{
  void
  Manifolds::declareDefaultEntry(ParameterHandler &prm, const unsigned int i_bc)
  {
    prm.declare_entry("type",
                      enum_to_string(ManifoldType::none),
                      Patterns::Selection(
                        "none|spherical|cylindrical|iges|step"),
                      "Type of manifold description"
                      "Choices are <none|spherical|cylindrical|iges|step>.");

    prm.declare_entry(
      "id",
      Utilities::int_to_string(i_bc, 2),
      Patterns::List(Patterns::Integer()),
      "IDs of boundaries to which the manifold will be applied");

    prm.declare_entry("cad file",
                      "none",
                      Patterns::FileName(),
                      "CAD file name (IGES or STEP)");
    prm.declare_entry(
      "cad scale factor",
      "1e-3",
      Patterns::Double(),
      "Scale factor applied to the CAD file once read. OpenCASCADE reads CAD "
      "files in millimeters, so the default (1e-3) yields meters");
    prm.declare_entry(
      "point coordinates",
      "0,0,0",
      Patterns::List(Patterns::Double()),
      "Point coordinates describing the spherical or cylindrical manifold");
    prm.declare_entry("direction vector",
                      "0,1,0",
                      Patterns::List(Patterns::Double()),
                      "Central axis describing the cylindrical manifold");
  }

  void
  Manifolds::parse_boundary(const ParameterHandler &prm)
  {
    // Parse the list of boundary IDs for the current manifold type
    std::vector<unsigned int> ids_current =
      convert_string_to_vector<unsigned int>(prm, "id");

    for (const unsigned int id : ids_current)
      {
        this->types.emplace_back(string_to_enum<ManifoldType>(prm.get("type")));
        this->ids.emplace_back(id);
        this->manifold_point.emplace_back(prm.get("point coordinates"));
        this->manifold_direction.emplace_back(prm.get("direction vector"));
        this->cad_files.emplace_back(prm.get("cad file"));
        this->cad_scale_factors.emplace_back(
          prm.get_double("cad scale factor"));
      }
  }

  void
  Manifolds::declare_parameters(ParameterHandler &prm,
                                unsigned int      subsection_max_size)
  {
    prm.enter_subsection("manifolds");
    {
      prm.declare_entry("number",
                        "0",
                        Patterns::Integer(),
                        "Number of manifolds");
      for (unsigned int i = 0; i < subsection_max_size; i++)
        {
          prm.enter_subsection("manifold " + Utilities::int_to_string(i));
          declareDefaultEntry(prm, i);
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  void
  Manifolds::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("manifolds");
    {
      this->number_of_manifolds_subsections = prm.get_integer("number");

      for (unsigned int i = 0; i < this->number_of_manifolds_subsections; i++)
        {
          prm.enter_subsection("manifold " + Utilities::int_to_string(i));
          parse_boundary(prm);
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }
} // namespace Parameters

template <int dim, int spacedim>
void
attach_manifolds_to_triangulation(
  parallel::DistributedTriangulationBase<dim, spacedim> &triangulation,
  Parameters::Manifolds                                  manifolds)
{
  for (unsigned int i = 0; i < manifolds.types.size(); ++i)
    {
      if (manifolds.types[i] == Parameters::Manifolds::ManifoldType::spherical)
        {
          // Create a point using the parameter file input
          Point<spacedim> circle_center(
            value_string_to_tensor<spacedim>(manifolds.manifold_point[i]));

          SphericalManifold<dim, spacedim> manifold_description(circle_center);
          triangulation.set_manifold(manifolds.ids[i], manifold_description);
        }
      else if (manifolds.types[i] ==
               Parameters::Manifolds::ManifoldType::cylindrical)
        {
          if constexpr (spacedim == 3)
            {
              // Create a point using the parameter file input
              Point<spacedim> point_on_axis(
                value_string_to_tensor<spacedim>(manifolds.manifold_point[i]));

              // Create a tensor representing the direction of the length of the
              // cylinder
              Tensor<1, spacedim> cylinder_axis(
                value_string_to_tensor<spacedim>(
                  manifolds.manifold_direction[i]));

              CylindricalManifold<dim, spacedim> manifold_description(
                cylinder_axis, point_on_axis);
              triangulation.set_manifold(manifolds.ids[i],
                                         manifold_description);
            }
          else
            throw std::runtime_error(
              "Cylindrical manifolds are not supported in 2D");
        }
      else if (manifolds.types[i] ==
                 Parameters::Manifolds::ManifoldType::iges ||
               manifolds.types[i] == Parameters::Manifolds::ManifoldType::step)
        {
          attach_cad_to_manifold(triangulation,
                                 manifolds.cad_files[i],
                                 manifolds.ids[i],
                                 manifolds.types[i],
                                 manifolds.cad_scale_factors[i]);
        }
      else if (manifolds.types[i] == Parameters::Manifolds::ManifoldType::none)
        {
        }
      else
        throw std::runtime_error("Unsupported manifolds type");
    }
}

void
attach_cad_to_manifold(parallel::DistributedTriangulationBase<2> &,
                       const std::string &,
                       const unsigned int,
                       const Parameters::Manifolds::ManifoldType,
                       const double)
{
  throw std::runtime_error("CAD manifolds are not supported in 2D");
}

void
attach_cad_to_manifold(parallel::DistributedTriangulationBase<2, 3> &,
                       const std::string &,
                       const unsigned int,
                       const Parameters::Manifolds::ManifoldType,
                       const double)
{
  throw std::runtime_error("CAD manifolds are not supported in 2D/3D");
}

#ifdef DEAL_II_WITH_OPENCASCADE
void
attach_cad_to_manifold(parallel::DistributedTriangulationBase<3> &triangulation,
                       const std::string                         &cad_name,
                       const unsigned int                         manifold_id,
                       const Parameters::Manifolds::ManifoldType  cad_type,
                       const double cad_scale_factor)
{
  TopoDS_Shape cad_surface =
    (cad_type == Parameters::Manifolds::ManifoldType::step) ?
      OpenCASCADE::read_STEP(cad_name, cad_scale_factor) :
      OpenCASCADE::read_IGES(cad_name, cad_scale_factor);

  // Enforce manifold over boundary ID
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      for (const auto &face : cell->face_iterators())
        {
          if (face->boundary_id() == manifold_id)
            {
              face->set_all_manifold_ids(manifold_id);
            }
        }
    }

  // Define tolerance for interpretation of CAD file
  const double tolerance = OpenCASCADE::get_shape_tolerance(cad_surface) * 5;

  //  OpenCASCADE::NormalProjectionManifold<3,3> normal_projector(
  //        cad_surface, tolerance);
  OpenCASCADE::NormalToMeshProjectionManifold<3, 3> normal_projector(
    cad_surface, tolerance);

  triangulation.set_manifold(manifold_id, normal_projector);
}
#else
void
attach_cad_to_manifold(parallel::DistributedTriangulationBase<3> &,
                       const std::string &,
                       const unsigned int,
                       const Parameters::Manifolds::ManifoldType,
                       const double)
{
  throw std::runtime_error(
    "CAD manifolds require DEAL_II to be compiled with OPENCASCADE");
}
#endif // DEAL_II_WITH_OPENCASCADE


template void
attach_manifolds_to_triangulation(
  parallel::DistributedTriangulationBase<2> &triangulation,
  Parameters::Manifolds                      manifolds);

template void
attach_manifolds_to_triangulation(
  parallel::DistributedTriangulationBase<3> &triangulation,
  Parameters::Manifolds                      manifolds);

template void
attach_manifolds_to_triangulation(
  parallel::DistributedTriangulationBase<2, 3> &triangulation,
  Parameters::Manifolds                         manifolds);
