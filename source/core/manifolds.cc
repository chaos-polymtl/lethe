// Dealii Includes
// Base
#include <deal.II/base/point.h>

// Grid
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include "core/manifolds.h"

namespace Parameters
{
  void
  Manifolds::declareDefaultEntry(ParameterHandler &prm, unsigned int i_bc)
  {
    prm.declare_entry("type",
                      "none",
                      Patterns::Selection("none|spherical|iges"),
                      "Type of manifold description"
                      "Choices are <none|spherical|iges>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(i_bc, 2),
                      Patterns::Integer(),
                      "Mesh id for boundary conditions");

    prm.declare_entry("cad file",
                      "none",
                      Patterns::FileName(),
                      "IGES file name");

    prm.declare_entry("arg1",
                      "0",
                      Patterns::Double(),
                      "Argument of construction no. 1");
    prm.declare_entry("arg2",
                      "0",
                      Patterns::Double(),
                      "Argument of construction no. 2");
    prm.declare_entry("arg3",
                      "0",
                      Patterns::Double(),
                      "Argument of construction no. 3");
    prm.declare_entry("arg4",
                      "0",
                      Patterns::Double(),
                      "Argument of construction no. 4");
    prm.declare_entry("arg5",
                      "0",
                      Patterns::Double(),
                      "Argument of construction no. 5");
    prm.declare_entry("arg6",
                      "0",
                      Patterns::Double(),
                      "Argument of construction no. 6");
  }

  void
  Manifolds::parse_boundary(ParameterHandler &prm, unsigned int i_bc)
  {
    const std::string op = prm.get("type");
    if (op == "none")
      types[i_bc] = ManifoldType::none;
    else if (op == "spherical")
      types[i_bc] = ManifoldType::spherical;
    else if (op == "iges")
      types[i_bc] = ManifoldType::iges;

    id[i_bc]        = prm.get_integer("id");
    arg1[i_bc]      = prm.get_double("arg1");
    arg2[i_bc]      = prm.get_double("arg2");
    arg3[i_bc]      = prm.get_double("arg3");
    arg4[i_bc]      = prm.get_double("arg4");
    arg5[i_bc]      = prm.get_double("arg5");
    arg6[i_bc]      = prm.get_double("arg6");
    cad_files[i_bc] = prm.get("cad file");
  }

  void
  Manifolds::declare_parameters(ParameterHandler &prm)
  {
    max_size = 7;
    arg1.resize(max_size);
    arg2.resize(max_size);
    arg3.resize(max_size);
    arg4.resize(max_size);
    arg5.resize(max_size);
    arg6.resize(max_size);
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
    }
    prm.leave_subsection();
  }
} // namespace Parameters

template <int dim, int spacedim>
void
attach_manifolds_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim, spacedim>>
                        triangulation,
  Parameters::Manifolds manifolds)
{
  for (unsigned int i = 0; i < manifolds.types.size(); ++i)
    {
      if (manifolds.types[i] == Parameters::Manifolds::ManifoldType::spherical)
        {
          Point<spacedim> circleCenter;
          circleCenter = Point<spacedim>(manifolds.arg1[i], manifolds.arg2[i]);
          static const SphericalManifold<dim, spacedim> manifold_description(
            circleCenter);
          triangulation->set_manifold(manifolds.id[i], manifold_description);
          triangulation->set_all_manifold_ids_on_boundary(manifolds.id[i],
                                                          manifolds.id[i]);
        }
      else if (manifolds.types[i] == Parameters::Manifolds::ManifoldType::iges)
        {
          attach_cad_to_manifold(triangulation,
                                 manifolds.cad_files[i],
                                 manifolds.id[i]);
        }
      else if (manifolds.types[i] == Parameters::Manifolds::ManifoldType::none)
        {}
      else
        throw std::runtime_error("Unsupported manifolds type");
    }
}

void
attach_cad_to_manifold(
  std::shared_ptr<parallel::DistributedTriangulationBase<2>>,
  std::string,
  unsigned int)
{
  throw std::runtime_error("IGES manifolds are not supported in 2D");
}

void
attach_cad_to_manifold(
  std::shared_ptr<parallel::DistributedTriangulationBase<2, 3>>,
  std::string,
  unsigned int)
{
  throw std::runtime_error("IGES manifolds are not supported in 2D/3D");
}

void
attach_cad_to_manifold(
  std::shared_ptr<parallel::DistributedTriangulationBase<3>> triangulation,
  std::string                                                cad_name,
  unsigned int                                               manifold_id)
{
#ifdef DEAL_II_WITH_OPENCASCADE

  TopoDS_Shape cad_surface = OpenCASCADE::read_IGES(cad_name, 1e-3);

  // Enforce manifold over boundary ID
  for (const auto &cell : triangulation->active_cell_iterators())
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

  triangulation->set_manifold(manifold_id, normal_projector);
#else
  throw std::runtime_error(
    "IGES manifolds require DEAL_II to be compiled with OPENCASCADE");

#endif // DEAL_II_WITH_OPENCASCADE
}


template void
attach_manifolds_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<2>> triangulation,
  Parameters::Manifolds                                      manifolds);

template void
attach_manifolds_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<3>> triangulation,
  Parameters::Manifolds                                      manifolds);

template void
attach_manifolds_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<2, 3>> triangulation,
  Parameters::Manifolds                                         manifolds);
