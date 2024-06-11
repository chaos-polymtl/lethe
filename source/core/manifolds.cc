#include <core/manifolds.h>

// Dealii Includes
#include <deal.II/base/point.h>

#include <deal.II/grid/manifold_lib.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

namespace Parameters
{
  void
  Manifolds::declareDefaultEntry(ParameterHandler &prm, unsigned int i_bc)
  {
    prm.declare_entry("type",
                      "none",
                      Patterns::Selection("none|spherical|cylindrical|iges"),
                      "Type of manifold description"
                      "Choices are <none|spherical|cylindrical|iges>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(i_bc, 2),
                      Patterns::Integer(),
                      "Mesh id for boundary conditions");

    prm.declare_entry("cad file",
                      "none",
                      Patterns::FileName(),
                      "IGES file name");
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
  Manifolds::parse_boundary(ParameterHandler &prm, unsigned int i_bc)
  {
    const std::string op = prm.get("type");
    if (op == "none")
      types[i_bc] = ManifoldType::none;
    else if (op == "spherical")
      types[i_bc] = ManifoldType::spherical;
    else if (op == "cylindrical")
      types[i_bc] = ManifoldType::cylindrical;
    else if (op == "iges")
      types[i_bc] = ManifoldType::iges;

    id[i_bc]                 = prm.get_integer("id");
    manifold_point[i_bc]     = prm.get("point coordinates");
    manifold_direction[i_bc] = prm.get("direction vector");
    cad_files[i_bc]          = prm.get("cad file");
  }

  void
  Manifolds::declare_parameters(ParameterHandler &prm)
  {
    max_size = 14;
    manifold_point.resize(max_size);
    manifold_direction.resize(max_size);
    cad_files.resize(max_size);

    prm.enter_subsection("manifolds");
    {
      prm.declare_entry("number",
                        "0",
                        Patterns::Integer(),
                        "Number of boundary conditions");
      id.resize(max_size);
      types.resize(max_size);

      prm.enter_subsection("manifold 0");
      declareDefaultEntry(prm, 0);
      prm.leave_subsection();

      prm.enter_subsection("manifold 1");
      declareDefaultEntry(prm, 1);
      prm.leave_subsection();

      prm.enter_subsection("manifold 2");
      declareDefaultEntry(prm, 2);
      prm.leave_subsection();

      prm.enter_subsection("manifold 3");
      declareDefaultEntry(prm, 3);
      prm.leave_subsection();

      prm.enter_subsection("manifold 4");
      declareDefaultEntry(prm, 4);
      prm.leave_subsection();

      prm.enter_subsection("manifold 5");
      declareDefaultEntry(prm, 5);
      prm.leave_subsection();

      prm.enter_subsection("manifold 6");
      declareDefaultEntry(prm, 6);
      prm.leave_subsection();

      prm.enter_subsection("manifold 7");
      declareDefaultEntry(prm, 7);
      prm.leave_subsection();

      prm.enter_subsection("manifold 8");
      declareDefaultEntry(prm, 8);
      prm.leave_subsection();

      prm.enter_subsection("manifold 9");
      declareDefaultEntry(prm, 9);
      prm.leave_subsection();

      prm.enter_subsection("manifold 10");
      declareDefaultEntry(prm, 10);
      prm.leave_subsection();

      prm.enter_subsection("manifold 11");
      declareDefaultEntry(prm, 11);
      prm.leave_subsection();

      prm.enter_subsection("manifold 12");
      declareDefaultEntry(prm, 12);
      prm.leave_subsection();

      prm.enter_subsection("manifold 13");
      declareDefaultEntry(prm, 13);
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  void
  Manifolds::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("manifolds");
    {
      size = prm.get_integer("number");
      types.resize(size);
      id.resize(size);
      if (size >= 1)
        {
          prm.enter_subsection("manifold 0");
          parse_boundary(prm, 0);
          prm.leave_subsection();
        }
      if (size >= 2)
        {
          prm.enter_subsection("manifold 1");
          parse_boundary(prm, 1);
          prm.leave_subsection();
        }
      if (size >= 3)
        {
          prm.enter_subsection("manifold 2");
          parse_boundary(prm, 2);
          prm.leave_subsection();
        }
      if (size >= 4)
        {
          prm.enter_subsection("manifold 3");
          parse_boundary(prm, 3);
          prm.leave_subsection();
        }
      if (size >= 5)
        {
          prm.enter_subsection("manifold 4");
          parse_boundary(prm, 4);
          prm.leave_subsection();
        }

      if (size >= 6)
        {
          prm.enter_subsection("manifold 5");
          parse_boundary(prm, 5);
          prm.leave_subsection();
        }
      if (size >= 7)
        {
          prm.enter_subsection("manifold 6");
          parse_boundary(prm, 6);
          prm.leave_subsection();
        }
      if (size >= 8)
        {
          prm.enter_subsection("manifold 7");
          parse_boundary(prm, 7);
          prm.leave_subsection();
        }
      if (size >= 9)
        {
          prm.enter_subsection("manifold 8");
          parse_boundary(prm, 8);
          prm.leave_subsection();
        }
      if (size >= 10)
        {
          prm.enter_subsection("manifold 9");
          parse_boundary(prm, 9);
          prm.leave_subsection();
        }
      if (size >= 11)
        {
          prm.enter_subsection("manifold 10");
          parse_boundary(prm, 10);
          prm.leave_subsection();
        }
      if (size >= 12)
        {
          prm.enter_subsection("manifold 11");
          parse_boundary(prm, 11);
          prm.leave_subsection();
        }
      if (size >= 13)
        {
          prm.enter_subsection("manifold 12");
          parse_boundary(prm, 12);
          prm.leave_subsection();
        }
      if (size >= 14)
        {
          prm.enter_subsection("manifold 13");
          parse_boundary(prm, 13);
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
            entry_string_to_tensor<spacedim>(manifolds.manifold_point[i]));

          static const SphericalManifold<dim, spacedim> manifold_description(
            circle_center);
          triangulation.set_manifold(manifolds.id[i], manifold_description);
        }
      else if (manifolds.types[i] ==
               Parameters::Manifolds::ManifoldType::cylindrical)
        {
          if constexpr (spacedim == 3)
            {
              // Create a point using the parameter file input
              Point<spacedim> point_on_axis(
                entry_string_to_tensor<spacedim>(manifolds.manifold_point[i]));

              // Create a tensor representing the direction of the length of the
              // cylinder
              Tensor<1, spacedim> cylinder_axis(
                entry_string_to_tensor<spacedim>(
                  manifolds.manifold_direction[i]));

              static const CylindricalManifold<dim, spacedim>
                manifold_description(cylinder_axis, point_on_axis);
              triangulation.set_manifold(manifolds.id[i], manifold_description);
            }
          else
            throw std::runtime_error(
              "Cylindrical manifolds are not supported in 2D");
        }
      else if (manifolds.types[i] == Parameters::Manifolds::ManifoldType::iges)
        {
          attach_cad_to_manifold(triangulation,
                                 manifolds.cad_files[i],
                                 manifolds.id[i]);
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
                       std::string,
                       unsigned int)
{
  throw std::runtime_error("IGES manifolds are not supported in 2D");
}

void
attach_cad_to_manifold(parallel::DistributedTriangulationBase<2, 3> &,
                       std::string,
                       unsigned int)
{
  throw std::runtime_error("IGES manifolds are not supported in 2D/3D");
}

#ifdef DEAL_II_WITH_OPENCASCADE
void
attach_cad_to_manifold(parallel::DistributedTriangulationBase<3> &triangulation,
                       std::string                                cad_name,
                       unsigned int                               manifold_id)
{
  TopoDS_Shape cad_surface = OpenCASCADE::read_IGES(cad_name, 1e-3);

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
                       std::string,
                       unsigned int)
{
  throw std::runtime_error(
    "IGES manifolds require DEAL_II to be compiled with OPENCASCADE");
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
