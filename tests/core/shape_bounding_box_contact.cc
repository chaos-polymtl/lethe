// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test the bounding box check of the contact between shapes, just
 * before and just after they touch, in both orders of the shapes.
 *
 * - Plane and rotated box, for an axis-aligned and for an oblique plane. They
 *   touch when the center of the box is at the sum over its axes of
 *   h_j |n . B_j| from the plane.
 * - Boxes of the same orientation, along each of their axes k. They touch at
 *   2 h_k.
 * - A small box in front of each face k of a rotated box. They touch at
 *   h_k + sum over j of s_j |A_k . B_j|.
 * - Two unit cubes whose edges cross, in a rotated frame. They touch at
 *   sqrt(2), and only an axis which is the cross product of two edges separates
 *   them.
 *
 * A_k and B_j are the columns of the rotation matrices of the boxes.
 */

// Deal.II includes
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/physics/transformations.h>

// Lethe
#include <core/shape.h>
#include <core/tensors_and_points_dimension_manipulation.h>

// Tests
#include <../tests/tests.h>

#include <cmath>
#include <numbers>
#include <string>

using namespace dealii;

// Relative distance from the position at which the shapes touch
constexpr double relative_offset = 1e-6;

/**
 * @brief Build a box with its rotation matrix.
 *
 * @param[in] half_lengths Half lengths of the box.
 * @param[in] position Center of the box.
 * @param[in] orientation XYZ angles of the box.
 *
 * @return The box.
 */
template <int dim>
HyperRectangle<dim>
make_box(const Tensor<1, dim> &half_lengths,
         const Point<dim>     &position,
         const Tensor<1, 3>   &orientation)
{
  HyperRectangle<dim> box(half_lengths, position, orientation);
  box.set_position(position);
  box.set_orientation(orientation);
  return box;
}

/**
 * @brief Get an axis of a shape in the global frame.
 *
 * @param[in] shape The shape.
 * @param[in] k Index of the axis.
 *
 * @return The unit vector along the axis.
 */
template <int dim>
Tensor<1, dim>
get_axis(Shape<dim> &shape, const unsigned int k)
{
  const Tensor<2, 3> rotation_matrix = shape.get_rotation_matrix();
  Tensor<1, dim>     axis;
  for (int d = 0; d < dim; ++d)
    axis[d] = rotation_matrix[d][k];
  return axis;
}

/**
 * @brief Check that two shapes are in contact just before the distance at which
 * they touch, and not just after it, in both orders.
 *
 * @param[in,out] fixed_shape Shape which stays in place.
 * @param[in,out] moved_shape Shape which is moved.
 * @param[in] origin Point from which the distance is measured.
 * @param[in] direction Unit vector along which the shape is moved.
 * @param[in] touching_distance Distance at which the shapes touch.
 *
 * @return true if the four checks are right.
 */
template <int dim>
bool
changes_when_touching(Shape<dim>           &fixed_shape,
                      Shape<dim>           &moved_shape,
                      const Point<dim>     &origin,
                      const Tensor<1, dim> &direction,
                      const double          touching_distance)
{
  moved_shape.set_position(origin + (1. - relative_offset) * touching_distance *
                                      direction);
  const bool contact_before = fixed_shape.bounding_box_contact(moved_shape) &&
                              moved_shape.bounding_box_contact(fixed_shape);

  moved_shape.set_position(origin + (1. + relative_offset) * touching_distance *
                                      direction);
  const bool contact_after = fixed_shape.bounding_box_contact(moved_shape) ||
                             moved_shape.bounding_box_contact(fixed_shape);

  return contact_before && !contact_after;
}

template <int dim>
void
test()
{
  auto log_check = [](const std::string &description, const bool is_right) {
    deallog << description << ": " << (is_right ? "true" : "false")
            << std::endl;
  };

  Point<dim>     center;
  Tensor<1, dim> half_lengths;
  Tensor<1, dim> small_half_lengths;
  center[0]             = 1.;
  center[1]             = 2.;
  half_lengths[0]       = 0.4;
  half_lengths[1]       = 0.15;
  small_half_lengths[0] = 0.03;
  small_half_lengths[1] = 0.05;

  Tensor<1, 3> orientation;
  Tensor<1, 3> other_orientation;
  if constexpr (dim == 2)
    {
      orientation[2]       = 0.5;
      other_orientation[2] = -0.7;
    }
  else
    {
      center[2]             = 3.;
      half_lengths[2]       = 0.25;
      small_half_lengths[2] = 0.02;
      orientation           = Tensor<1, 3>({0.3, 0.4, 0.5});
      other_orientation     = Tensor<1, 3>({-0.6, 0.2, 1.1});
    }

  deallog << "Dimension " << dim << std::endl;

  // Plane and box. Without rotation, the normal of a plane is y in 2D and z in
  // 3D.
  {
    auto       box      = make_box<dim>(half_lengths, center, orientation);
    const auto is_right = [&](const Tensor<1, 3> &plane_orientation) {
      Plane<dim> plane(center, plane_orientation);
      plane.set_position(center);
      plane.set_orientation(plane_orientation);
      const Tensor<1, dim> normal = get_axis<dim>(plane, dim - 1);

      double touching_distance = 0.;
      for (int j = 0; j < dim; ++j)
        touching_distance +=
          half_lengths[j] * std::abs(normal * get_axis<dim>(box, j));
      return changes_when_touching<dim>(
        plane, box, center, normal, touching_distance);
    };
    log_check("Axis-aligned plane and box", is_right(Tensor<1, 3>()));
    log_check("Oblique plane and box", is_right(other_orientation));
  }

  // Boxes of the same orientation
  {
    auto box   = make_box<dim>(half_lengths, center, orientation);
    auto other = make_box<dim>(half_lengths, center, orientation);

    bool is_right = true;
    for (int k = 0; k < dim; ++k)
      is_right = changes_when_touching<dim>(box,
                                            other,
                                            center,
                                            get_axis<dim>(box, k),
                                            2. * half_lengths[k]) &&
                 is_right;
    log_check("Boxes of the same orientation", is_right);
  }

  // Small box in front of each face of a rotated box
  {
    auto box   = make_box<dim>(half_lengths, center, orientation);
    auto small = make_box<dim>(small_half_lengths, center, other_orientation);

    bool is_right = true;
    for (int k = 0; k < dim; ++k)
      {
        const Tensor<1, dim> face_normal = get_axis<dim>(box, k);

        double touching_distance = half_lengths[k];
        for (int j = 0; j < dim; ++j)
          touching_distance += small_half_lengths[j] *
                               std::abs(face_normal * get_axis<dim>(small, j));
        is_right = changes_when_touching<dim>(
                     box, small, center, face_normal, touching_distance) &&
                   is_right;
      }
    log_check("Boxes of different orientations", is_right);
  }

  // Crossed edges
  if constexpr (dim == 3)
    {
      using Physics::Transformations::Rotations::rotation_matrix_3d;
      const double       angle = 0.25 * std::numbers::pi;
      const Tensor<2, 3> frame =
        make_box<dim>(half_lengths, center, orientation).get_rotation_matrix();

      Tensor<2, 3> rotation_one =
        frame * rotation_matrix_3d(Tensor<1, 3>({1., 0., 0.}), angle);
      Tensor<2, 3> rotation_two =
        frame * rotation_matrix_3d(Tensor<1, 3>({0., 1., 0.}), angle);

      const Tensor<1, 3> unit_cube({0.5, 0.5, 0.5});

      auto cube_one =
        make_box<dim>(unit_cube,
                      center,
                      Shape<dim>::rotation_matrix_to_xyz_angles(rotation_one));
      auto cube_two =
        make_box<dim>(unit_cube,
                      center,
                      Shape<dim>::rotation_matrix_to_xyz_angles(rotation_two));

      const Tensor<1, 3> direction = frame * Tensor<1, 3>({0., 0., 1.});
      log_check("Boxes with crossed edges",
                changes_when_touching<dim>(
                  cube_one, cube_two, center, direction, std::numbers::sqrt2));
    }
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      initlog();
      test<2>();
      test<3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
