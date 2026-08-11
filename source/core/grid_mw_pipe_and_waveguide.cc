#include <core/grid_mw_pipe_and_waveguide.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace dealii;

/*
This file creates a mesh that represents a rectangular volume entering a
cylinder laterally, perpendicular to its axis.

The program is divided into two parts:

1)The function "detailed_circle" is first defined, inspired by
hyper_ball in the deal.II library. This function creates a 2D coarse mesh of a
circle which is adapted to be connected with a rectangle on its side.

2)The core of the file creates the 3D pipe assembled with the waveguide
*/

// See the doxygen documentation to learn more about this function
template <int dim, int spacedim>
void
GridMwPipeAndWaveguide<dim, spacedim>::detailed_circle(
  Triangulation<2, 2> &tria,
  const Point<2>      &p,
  const double        &a,
  const double        &b,
  const double         r,
  const double         R,
  const bool           shape)
{
  const auto embed_point = [](const double x, const double y) -> Point<2> {
    return Point<2>(x, y);
  };

  const auto vertices = [&]() -> std::vector<Point<2>> {
    if (shape)
      // Inner circle
      return {
        p + embed_point(-1., -1.) * (r / std::sqrt(2.0)),
        p + embed_point(+1., -1.) * (r / std::sqrt(2.0)),
        p + embed_point(-1., +1.) * (r / std::sqrt(2.0)),
        p + embed_point(+1., +1.) * (r / std::sqrt(2.0)),
        // The four points which make the inner square
        p + embed_point(-1., +1.) * (R / std::sqrt(2.0)),
        p + embed_point(+1., +1.) * (R / std::sqrt(2.0)),
        p + embed_point(-1., -1.) * (R / std::sqrt(2.0)),
        p + embed_point(+1., -1.) * (R / std::sqrt(2.0)),
        // Two points linked with the rectangle on is left
        embed_point(b, 0),
        embed_point(b, a),
        // The two points (0, r) and (0, -r) from the center of the circle
        p + embed_point(0, -r),
        p + embed_point(0, +r),
        // The two points on the right side of the outer circle that lie along
        // the width of the rectangle
        p + embed_point(std::sqrt(r * r - a * a / 4), -a / 2),
        p + embed_point(std::sqrt(r * r - a * a / 4), +a / 2),
        // The other points where the width of the rectangle intersects the
        // central square
        p + embed_point(-(std::sqrt(R * R - a * a / 4)), +a / 2),
        p + embed_point(-(std::sqrt(R * R - a * a / 4)), -a / 2),
        p + embed_point(+(std::sqrt(R * R - a * a / 4)), -a / 2),
        p + embed_point(+(std::sqrt(R * R - a * a / 4)), +a / 2),
        // The two points of the central square on the circle's axis of symmetry
        p + embed_point(0, -R),
        p + embed_point(0, +R),
        // The two points where the width of the rectangle intersects the axis
        // of symmetry of the circle
        p + embed_point(0, +a / 2),
        p + embed_point(0, -a / 2)};

    else
      // Inner square
      return {
        p + embed_point(-1., -1.) * (r / std::sqrt(2.0)),
        p + embed_point(+1., -1.) * (r / std::sqrt(2.0)),
        p + embed_point(-1., +1.) * (r / std::sqrt(2.0)),
        p + embed_point(+1., +1.) * (r / std::sqrt(2.0)),
        // The other points that form the inner square
        p + embed_point(-1., -1.) * R,
        p + embed_point(+1., -1.) * R,
        p + embed_point(-1., +1.) * R,
        p + embed_point(+1., +1.) * R,
        // The two points attached to the rectangle
        embed_point(b, 0),
        embed_point(b, a),
        // The two points (0, r) and (0, -r) from the center of the circle
        p + embed_point(0, -r),
        p + embed_point(0, +r),
        // The two points on the right side of the outer circle that lie along
        // the width of the rectangle
        p + embed_point(std::sqrt(r * r - a * a / 4), -a / 2),
        p + embed_point(std::sqrt(r * r - a * a / 4), +a / 2),
        // The other points where the width of the rectangle intersects the
        // central square
        p + embed_point(-R, +a / 2),
        p + embed_point(-R, -a / 2),
        p + embed_point(+R, -a / 2),
        p + embed_point(+R, +a / 2),
        // The two points of the central square on the circle's axis of symmetry
        p + embed_point(0, -R),
        p + embed_point(0, +R),
        // The two points where the width of the rectangle intersects the axis
        // of symmetry of the circle
        p + embed_point(0, +a / 2),
        p + embed_point(0, -a / 2)};
  }();

  std::vector<CellData<2>> cells(16, CellData<2>());
  static constexpr int     circle_cell_vertices[16][4] = {
    // We start with the 6 squares in the center square
    {6, 18, 15, 21},
    {15, 21, 14, 20},
    {14, 20, 4, 19},
    {18, 7, 21, 16},
    {21, 16, 20, 17},
    {20, 17, 19, 5},
    // We fill in the outer cells, starting from the top left
    {9, 14, 2, 4},
    {8, 15, 9, 14},
    {0, 6, 8, 15},
    {0, 10, 6, 18},
    {10, 1, 18, 7},
    {7, 1, 16, 12},
    {16, 12, 17, 13},
    {17, 13, 5, 3},
    {19, 5, 11, 3},
    {4, 19, 2, 11}};

  for (unsigned int i = 0; i < 16; ++i)
    {
      for (unsigned int j = 0; j < 4; ++j)
        cells[i].vertices[j] = circle_cell_vertices[i][j];
    }

  tria.create_triangulation(std::vector<Point<2>>(std::begin(vertices),
                                                  std::end(vertices)),
                            cells,
                            SubCellData()); // no boundary information
  tria.set_all_manifold_ids_on_boundary(0);
  tria.set_manifold(0, SphericalManifold<2, 2>(p));
  tria.set_manifold(1, FlatManifold<2, 2>());
}



// The core of the file begins here

template <int dim, int spacedim>
GridMwPipeAndWaveguide<dim, spacedim>::GridMwPipeAndWaveguide(
  const std::string &grid_arguments)
{
  const std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');

  rectangle_width   = std::stod(arguments[0]);
  rectangle_length  = std::stod(arguments[1]);
  outer_radius      = std::stod(arguments[2]);
  inner_length      = std::stod(arguments[3]);
  bottom_height     = std::stod(arguments[4]);
  middle_height     = std::stod(arguments[5]);
  top_height        = std::stod(arguments[6]);
  bottom_resolution = std::stod(arguments[7]);
  middle_resolution = std::stod(arguments[8]);
  top_resolution    = std::stod(arguments[9]);
  shape             = (arguments[10] == "true");
}

template <int dim, int spacedim>
void
GridMwPipeAndWaveguide<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> &triangulation)
{
  /*
  Two geometric conditions that need to be verified to prevent the generation of
  inverted cells.

  The first condition concerns the inner shape. To improve mesh quality, a
  central rectangle of width equal to the waveguide width is introduced. After
  extrusion, the cells within this rectangle are nearly hexahedral, which
  minimizes aspect ratio. However, this rectangle imposes a geometric
  constraint: the inner shape must be large enough to fully contain it,
  otherwise inverted cells are generated.

  Then, the inner shape have to be strictly inscribed in the outer shape.
  */

  // Tests when the inner shape is a circle
  if (shape)
    {
      AssertThrow(
        inner_length > rectangle_width / 2,
        ExcMessage(
          "The inner shape is not bigger than the rectangle that crosses the cylinder along the x-axis. This situations leads to the generation of inverted cells"));

      AssertThrow(
        inner_length < outer_radius,
        ExcMessage(
          "The internal radius is greater than the external radius. This situations leads to the generation of inverted cells"));
    }

  // Tests when the inner shape is a square
  else
    {
      AssertThrow(
        inner_length > rectangle_width / 2,
        ExcMessage(
          "The square side is too small compared to the width of the waveguide. This situations leads to the generation of inverted cells"));

      AssertThrow(
        inner_length < outer_radius / (std::sqrt(2)),
        ExcMessage(
          "The inner square length is greater than the external radius. This situations leads to the generation of inverted cells"));
    }

  // Dimension Test
  if constexpr (!(dim == 3 && spacedim == 3))
    {
      AssertThrow(false, ExcMessage("GridCavityMw is only supported in 3D."));
      return;
    }

  else
    {
      // COARSE 2D PARTS OF THE MESH

      // PART -I-:Creating the Shape of the Medium: The paving stone embedded in
      // the Cylinder

      // Creating the first shape: a rectangle (after extrusion, this rectangle
      // will represent the paving stone). This rectangle is the reason why we
      // divide the volume into three different pieces.
      Triangulation<2> rect;
      GridGenerator::hyper_rectangle(
        rect, Point<2>(0.0, 0.0), Point<2>(rectangle_length, rectangle_width));

      // Creating the second shape: a rough circle with two points located
      // exactly at the right corners of the rectangle
      Triangulation<2> disk;
      // Position of the center of the circle
      const Point<2> center(rectangle_length +
                              sqrt(outer_radius * outer_radius -
                                   rectangle_width * rectangle_width / 4),
                            rectangle_width / 2);
      detailed_circle(disk,
                      center,
                      rectangle_width,
                      rectangle_length,
                      outer_radius,
                      inner_length,
                      shape);

      Triangulation<2> merged;
      // The last argument of the following function is
      // false since it is useless to copy the manifolds ids
      // because we will set them after
      GridGenerator::merge_triangulations(rect, disk, merged, 1e-12, false);

      // PART -II-   Creating the bottom cylinder
      Triangulation<dim - 1, spacedim - 1> base;
      detailed_circle(base,
                      center,
                      rectangle_width,
                      rectangle_length,
                      outer_radius,
                      inner_length,
                      shape);

      // PART -III-  Creating the top cylinder
      Triangulation<2, 2> top;
      detailed_circle(top,
                      center,
                      rectangle_width,
                      rectangle_length,
                      outer_radius,
                      inner_length,
                      shape);

      // EXTRUDE

      Triangulation<3> extruded_base;
      Triangulation<3> extruded_middle;
      Triangulation<3> extruded_top;
      GridGenerator::extrude_triangulation(base,
                                           bottom_resolution,
                                           bottom_height,
                                           extruded_base);
      GridGenerator::extrude_triangulation(merged,
                                           middle_resolution,
                                           middle_height,
                                           extruded_middle);
      GridGenerator::extrude_triangulation(top,
                                           top_resolution,
                                           top_height,
                                           extruded_top);

      // MANIFOLD

      // Lambda function that places the id of the cylindrical manifold in the
      // inside of the cylinder
      auto attach_cylindrical_manifold_id = [&](Triangulation<3> &tria,
                                                double            height) {
        for (auto &cell : tria.active_cell_iterators())
          if (cell->center()[0] >= rectangle_length - 1e-6)
            // Concerns only the cylinder because every cell of the
            // waveguide verifies the condition : x >= rectangle_length
            for (const auto f : cell->face_indices())
              if (cell->face(f)->at_boundary())
                {
                  const Point<3> face_center = cell->face(f)->center();
                  if (std::abs(face_center[2]) > 1e-10 &&
                      std::abs(face_center[2] - height) > 1e-10)
                    // Does not take into account the the flat faces
                    // at z=0 and z=height because we do not want a manifold
                    // applied there.
                    cell->set_all_manifold_ids(1);
                }
      };

      attach_cylindrical_manifold_id(extruded_base, bottom_height);
      attach_cylindrical_manifold_id(extruded_middle, middle_height);
      attach_cylindrical_manifold_id(extruded_top, top_height);

      // Now that we have the IDs in the right place, we can install the
      // manifolds
      const Point<3>     center_3d(center[0], center[1], 0.0);
      const Tensor<1, 3> axis({0.0, 0.0, 1.0});
      extruded_base.set_manifold(1, CylindricalManifold<3>(axis, center_3d));
      extruded_middle.set_manifold(1, CylindricalManifold<3>(axis, center_3d));
      extruded_top.set_manifold(1, CylindricalManifold<3>(axis, center_3d));

      // MERGING THE 3D MESHES

      // Move the cylinders in the right positions to prepare the merge
      Tensor<1, 3> shift_down;
      shift_down[0] = 0.0;
      shift_down[1] = 0.0;
      shift_down[2] = -bottom_height; // The base_cylinder must have its base at
                                      // -bottom_height
      GridTools::shift(shift_down, extruded_base);


      Tensor<1, 3> shift_up;
      shift_up[0] = 0.0;
      shift_up[1] = 0.0;
      shift_up[2] = +middle_height; // The top_cylinder  must have its base at
                                    // +middle_height
      GridTools::shift(shift_up, extruded_top);


      // Merge the 3 parts of the mesh
      std::vector<const Triangulation<3, 3> *> parts = {&extruded_base,
                                                        &extruded_middle,
                                                        &extruded_top};

      GridGenerator::merge_triangulations(
        parts,
        triangulation,
        1e-12,
        true); // “true” is used here to preserve the cell IDs; this is useful
               // for attaching the manifold right after

      // Attach the manifolds that we previously set on each part
      triangulation.set_manifold(1, CylindricalManifold<3>(axis, center_3d));


      // BOUNDARY AND MATERIAL IDs

      /*
      //Boundary
      0: MW Inlet (x=0)
      1: Fluid Inlet  (y = -bottom_height)
      2: Fluid Outlet (y = middle_height + top_height)
      3: PEC walls (every other wall)

      //Material
      0: Fluid part (x> rectangle_length)
      1: Waveguide paving stone (x<=  rectangle_length)
      */

      for (auto &cell : triangulation.active_cell_iterators())
        {
          // Boundary ids
          for (const auto f : cell->face_indices())
            {
              if (cell->face(f)->at_boundary())
                {
                  const Point<3> face_center = cell->face(f)->center();
                  if (face_center[0] <= 1e-6)
                    cell->face(f)->set_boundary_id(0); // MW Inlet
                  else if (face_center[2] <= -bottom_height + 1e-6)
                    cell->face(f)->set_boundary_id(1); // Fluid Inlet
                  else if (face_center[2] >= middle_height + top_height - 1e-6)
                    cell->face(f)->set_boundary_id(2); // Fluid Outlet
                  else
                    cell->face(f)->set_boundary_id(
                      3); // PEC Walls (Perfect Electric Conductors)
                }
            }
          // Material ids
          if (cell->center()[0] <= rectangle_length)
            cell->set_material_id(1); // waveguide (Paving stone)
          else
            cell->set_material_id(0); // Fluid (cylinder)
        }
    }
}

// Explicit instantiations
template class GridMwPipeAndWaveguide<2, 2>;
template class GridMwPipeAndWaveguide<2, 3>;
template class GridMwPipeAndWaveguide<3, 3>;
