#include <core/grid_pipe_and_waveguide.h>

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
#include <numbers>

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
GridPipeAndWaveguide<dim, spacedim>::detailed_circle(Triangulation<2, 2> &tria,
                                                     const Point<2>      &p,
                                                     const double        &a,
                                                     const double        &b,
                                                     const double         r,
                                                     const double         R,
                                                     const bool inner_circle)
{
  const auto embed_point = [](const double x, const double y) -> Point<2> {
    return Point<2>(x, y);
  };

  const auto vertices = [&]() -> std::vector<Point<2>> {
    if (inner_circle)
      // Inner circle
      return {
        p + embed_point(-1., -1.) * (r / std::numbers::sqrt2),
        p + embed_point(+1., -1.) * (r / std::numbers::sqrt2),
        p + embed_point(-1., +1.) * (r / std::numbers::sqrt2),
        p + embed_point(+1., +1.) * (r / std::numbers::sqrt2),
        // The four points which make the inner square
        p + embed_point(-1., +1.) * (R / std::numbers::sqrt2),
        p + embed_point(+1., +1.) * (R / std::numbers::sqrt2),
        p + embed_point(-1., -1.) * (R / std::numbers::sqrt2),
        p + embed_point(+1., -1.) * (R / std::numbers::sqrt2),
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
        p + embed_point(-1., -1.) * (r / std::numbers::sqrt2),
        p + embed_point(+1., -1.) * (r / std::numbers::sqrt2),
        p + embed_point(-1., +1.) * (r / std::numbers::sqrt2),
        p + embed_point(+1., +1.) * (r / std::numbers::sqrt2),
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
}



// The core of the file begins here

template <int dim, int spacedim>
GridPipeAndWaveguide<dim, spacedim>::GridPipeAndWaveguide(
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
  bottom_resolution = std::stoi(arguments[7]);
  middle_resolution = std::stoi(arguments[8]);
  top_resolution    = std::stoi(arguments[9]);
  inner_circle      = (arguments[10] == "true");
}

template <int dim, int spacedim>
void
GridPipeAndWaveguide<dim, spacedim>::make_grid(
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
  if (inner_circle)
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
      AssertThrow(false,
                  ExcMessage("GridPipeAndWaveguide is only supported in 3D."));
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
                      inner_circle);

      Triangulation<2> merged;
      // The last argument of the following function is
      // false since it is useless to copy the manifolds ids
      // because we will set them after
      GridGenerator::merge_triangulations(rect, disk, merged, 1e-12, false);

      // PART -II-   Creating the bottom cylinder
      Triangulation<2, 2> base;
      detailed_circle(base,
                      center,
                      rectangle_width,
                      rectangle_length,
                      outer_radius,
                      inner_length,
                      inner_circle);

      // PART -III-  Creating the top cylinder
      Triangulation<2, 2> top;
      detailed_circle(top,
                      center,
                      rectangle_width,
                      rectangle_length,
                      outer_radius,
                      inner_length,
                      inner_circle);

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

      // The manifold ids are not copied since they are set on the merged mesh
      GridGenerator::merge_triangulations(parts, triangulation, 1e-12, false);


      // BOUNDARY AND MATERIAL IDs

      /*
      Material ids:
      - cylinder_material_id: cells of the cylinder (x > rectangle_length)
      - waveguide_material_id: cells of the waveguide (x <= rectangle_length)

      Boundary ids:
      - waveguide_port_boundary_id: end face of the waveguide (x = 0)
      - waveguide_wall_boundary_id: other boundary faces of the waveguide
      - bottom_boundary_id: bottom surface of the cylinder
        (z = -bottom_height)
      - top_boundary_id: top surface of the cylinder
        (z = middle_height + top_height)
      - cylinder_wall_boundary_id: lateral wall of the cylinder

      The material id of a cell is set first, since it is used to distinguish
      the walls of the waveguide from the walls of the cylinder.
      */

      for (auto &cell : triangulation.active_cell_iterators())
        {
          const bool cell_is_in_waveguide =
            cell->center()[0] <= rectangle_length;
          cell->set_material_id(cell_is_in_waveguide ? waveguide_material_id :
                                                       cylinder_material_id);

          for (const auto f : cell->face_indices())
            {
              if (!cell->face(f)->at_boundary())
                continue;

              const Point<3> face_center = cell->face(f)->center();
              if (face_center[0] <= 1e-6)
                cell->face(f)->set_boundary_id(waveguide_port_boundary_id);
              else if (cell_is_in_waveguide)
                cell->face(f)->set_boundary_id(waveguide_wall_boundary_id);
              else if (face_center[2] <= -bottom_height + 1e-6)
                cell->face(f)->set_boundary_id(bottom_boundary_id);
              else if (face_center[2] >= middle_height + top_height - 1e-6)
                cell->face(f)->set_boundary_id(top_boundary_id);
              else
                cell->face(f)->set_boundary_id(cylinder_wall_boundary_id);
            }
        }


      // MANIFOLD

      /*
      The cylindrical manifold is attached to every cell of the cylinder that
      has a face on its lateral wall, that is a face whose vertices are all at a
      distance outer_radius from the axis of the cylinder. The manifold id is
      set on the whole cell, so that its interior lines and faces also follow
      the cylinder when the mesh is refined. This includes the cells of the
      middle section that are next to the waveguide: their face shared with the
      waveguide is an interior face of the mesh, but it lies on the lateral
      wall of the cylinder. The manifold ids are set on the merged mesh so that
      they are the same in the three sections. The cells of the waveguide keep
      the flat manifold, so that the walls of the waveguide stay planar.
      */
      const double distance_to_axis_tolerance = 1e-8 * outer_radius;
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->material_id() != cylinder_material_id)
            continue;

          for (const auto f : cell->face_indices())
            {
              bool face_is_on_lateral_wall = true;
              for (const auto v : cell->face(f)->vertex_indices())
                {
                  const Point<3> &vertex = cell->face(f)->vertex(v);
                  const double    distance_to_axis =
                    center.distance(Point<2>(vertex[0], vertex[1]));
                  if (std::abs(distance_to_axis - outer_radius) >
                      distance_to_axis_tolerance)
                    face_is_on_lateral_wall = false;
                }

              if (face_is_on_lateral_wall)
                {
                  cell->set_all_manifold_ids(cylindrical_manifold_id);
                  break;
                }
            }
        }

      const Point<3>     center_3d(center[0], center[1], 0.0);
      const Tensor<1, 3> axis({0.0, 0.0, 1.0});
      triangulation.set_manifold(cylindrical_manifold_id,
                                 CylindricalManifold<3>(axis, center_3d));
    }
}

// Explicit instantiations
template class GridPipeAndWaveguide<2, 2>;
template class GridPipeAndWaveguide<2, 3>;
template class GridPipeAndWaveguide<3, 3>;
