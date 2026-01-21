// SPDX-FileCopyrightText: Copyright (c) 2021 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/config.h>

using namespace dealii;


/*-------------- Mortar manager base -------------------------------*/
/**
 * @brief Base class for the mortar manager
 */
template <int dim>
class MortarManagerBase
{
public:
  template <int dim2>
  MortarManagerBase(unsigned int            n_subdivisions,
                    double                  radius,
                    const Quadrature<dim2> &quadrature,
                    const double            rotation_angle);

  template <int dim2>
  MortarManagerBase(const std::vector<unsigned int> &n_subdivisions,
                    const std::vector<double>       &radius,
                    const Quadrature<dim2>          &quadrature,
                    const double                     rotation_angle);

  /**
   * @brief Default destructor.
   */
  virtual ~MortarManagerBase() = default;

  /**
   * @brief Verify if cells of the inner and outer domains are aligned
   */
  bool
  is_mesh_aligned() const
  {
    AssertThrow(dim != 1, ExcInternalError());

    constexpr double tolerance = 1e-8;
    const double     delta_0   = 2 * numbers::PI / n_subdivisions[0];

    return std::abs(rotation_angle / delta_0 -
                    std::round(rotation_angle / delta_0)) < tolerance;
  }

  /**
   * @brief Returns the total number of mortars
   */
  unsigned int
  get_n_total_mortars() const
  {
    if constexpr (dim == 1)
      return 1;

    unsigned int n_total_subdivisions = n_subdivisions[0];

    // In 3D, besides the subdivisions in the plane perpendicular to the
    // rotation axis, we also need to account for the number of subdivisions
    // along the rotation axis direction
    if constexpr (dim == 3)
      n_total_subdivisions *= n_subdivisions[1];

    if (this->is_mesh_aligned()) // aligned
      {
        return n_total_subdivisions;
      }
    else // inside/outside
      {
        return 2 * n_total_subdivisions;
      }
  }

  /**
   * @brief Returns the number of mortars per face
   */
  unsigned int
  get_n_mortars() const
  {
    if constexpr (dim == 1)
      return 1;

    if (this->is_mesh_aligned()) // aligned
      {
        return 1;
      }
    else // inside/outside
      {
        return 2;
      }
  }

  /**
   * @brief Returns the indices of all mortars at both sides of the interface
   *
   * @param[in] face_center Face center
   */
  std::vector<unsigned int>
  get_mortar_indices(const Point<dim> &face_center, const bool is_inner) const
  {
    if constexpr (dim == 1)
      return std::vector<unsigned int>{0};

    // Mesh alignment type and cell indexes
    const auto [type, id_in_plane, id_out_plane] =
      get_config(face_center, is_inner);

    if (type == 0) // aligned
      {
        std::vector<unsigned int> indices;

        const unsigned int index = id_in_plane;

        AssertIndexRange(index, n_subdivisions[0]);

        indices.emplace_back(index + n_subdivisions[0] * id_out_plane);

        return indices;
      }
    else if (type == 1) // inside
      {
        std::vector<unsigned int> indices;

        for (unsigned int q = 0; q < 2; ++q)
          {
            const unsigned int index =
              (id_in_plane * 2 + 1 + q) % (n_subdivisions[0] * 2);

            AssertIndexRange(index, n_subdivisions[0] * 2);

            indices.emplace_back(index + 2 * n_subdivisions[0] * id_out_plane);
          }

        return indices;
      }
    else // outside
      {
        std::vector<unsigned int> indices;

        for (unsigned int q = 0; q < 2; ++q)
          {
            const unsigned int index = id_in_plane * 2 + q;

            AssertIndexRange(index, n_subdivisions[0] * 2);

            indices.emplace_back(index + 2 * n_subdivisions[0] * id_out_plane);
          }

        return indices;
      }
  }


  /**
   * @brief Returns the total number of quadrature points at the inner/outer boundary interface
   */
  unsigned int
  get_n_total_points() const
  {
    if constexpr (dim == 1)
      return 1;

    return get_n_total_mortars() * n_quadrature_points;
  }

  /**
   * @brief Returns the coordinates of the quadrature points at both sides of the interface
   */
  unsigned int
  get_n_points() const
  {
    if constexpr (dim == 1)
      return 1;

    return get_n_mortars() * n_quadrature_points;
  }

  /**
   * @brief Returns the coordinates of the quadrature points at both sides of the interface
   *
   * @param[in] face_center Face center
   *
   * @return points Coordinate of quadrature points of the cell
   */
  std::vector<Point<dim>>
  get_points(const Point<dim> &face_center, const bool is_inner) const
  {
    if constexpr (dim == 1)
      return std::vector<Point<dim>>{face_center};

    // Mesh alignment type and cell index
    const auto [type, id_in_plane, id_out_plane] =
      get_config(face_center, is_inner);
    // Angle variation within each cell
    const double delta_0 = 2 * numbers::PI / n_subdivisions[0];
    double       delta_1 = 1.0;

    if constexpr (dim == 3)
      delta_1 = radius[1] / n_subdivisions[1];

    if (type == 0) // aligned
      {
        std::vector<Point<dim>> points;

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          {
            const auto x =
              from_1D((id_in_plane + quadrature.point(q)[0]) * delta_0);

            if constexpr (dim == 3)
              points.emplace_back(
                x[0], x[1], (id_out_plane + quadrature.point(q)[1]) * delta_1);
            else
              points.emplace_back(x);
          }

        return points;
      }
    else // Point at the inner boundary lies somewhere in the face of the outer
         // boundary cell
      {
        // rad_0: first cell vertex (fixed)
        // rad_1: shifted vertex
        // rad_2: last cell vertex (fixed)
        double rad_0, rad_1, rad_2;
        // Minimum rotation angle
        double rot_min =
          rotation_angle - std::floor(rotation_angle / delta_0) * delta_0;

        if (type == 2) // outside
          {
            rad_0 = id_in_plane * delta_0;
            rad_1 = id_in_plane * delta_0 + rot_min;
            rad_2 = (id_in_plane + 1) * delta_0;
          }
        else // inside
          {
            rad_0 = id_in_plane * delta_0 + rot_min;
            rad_1 = (id_in_plane + 1) * delta_0;
            rad_2 = (id_in_plane + 1) * delta_0 + rot_min;
          }

        std::vector<Point<dim>> points;

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          {
            const auto x =
              from_1D(rad_0 + quadrature.point(q)[0] * (rad_1 - rad_0));

            if constexpr (dim == 3)
              points.emplace_back(
                x[0], x[1], (id_out_plane + quadrature.point(q)[1]) * delta_1);
            else
              points.emplace_back(x);
          }

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          {
            const auto x =
              from_1D(rad_1 + quadrature.point(q)[0] * (rad_2 - rad_1));

            if constexpr (dim == 3)
              points.emplace_back(
                x[0], x[1], (id_out_plane + quadrature.point(q)[1]) * delta_1);
            else
              points.emplace_back(x);
          }

        return points;
      }
  }


  /**
   * @brief Returns the coordinates of the quadrature points at the interface
   *
   * @param[in] face_center Face center
   *
   * @return points Coordinate of quadrature points of the cell
   */
  std::vector<Point<std::max(1, dim - 1)>>
  get_points_ref(const Point<dim> &face_center, const bool is_inner) const
  {
    if (dim == 1)
      return std::vector<Point<std::max(1, dim - 1)>>{
        Point<std::max(1, dim - 1)>()};

    const auto [type, _, __] = get_config(face_center, is_inner);

    const double delta_0 = 2 * numbers::PI / n_subdivisions[0];

    if (type == 0) // aligned
      {
        std::vector<Point<std::max(1, dim - 1)>> points;
        points.reserve(n_quadrature_points);

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          points.emplace_back(quadrature.point(q));

        return points;
      }
    else // inside/outside
      {
        double rad_0, rad_1, rad_2;

        double rot_min =
          (rotation_angle - std::floor(rotation_angle / delta_0) * delta_0) /
          delta_0;

        if (type == 2) // outside
          {
            rad_0 = 0.0;
            rad_1 = rot_min;
            rad_2 = 1.0;
          }
        else // inside
          {
            rad_0 = 0.0;
            rad_1 = 1.0 - rot_min;
            rad_2 = 1.0;
          }

        std::vector<Point<std::max(1, dim - 1)>> points;
        points.reserve(2 * n_quadrature_points);

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          {
            const double x = rad_0 + quadrature.point(q)[0] * (rad_1 - rad_0);

            if (dim == 2)
              points.emplace_back(x);
            else
              points.emplace_back(x, quadrature.point(q)[1]);
          }

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          {
            const double x = rad_1 + quadrature.point(q)[0] * (rad_2 - rad_1);

            if (dim == 2)
              points.emplace_back(x);
            else
              points.emplace_back(x, quadrature.point(q)[1]);
          }

        return points;
      }
  }


  /**
   * @brief Returns the weights of the quadrature points at both sides of the interface
   *
   * @param[in] face_center Face center
   *
   * @return points Angular weights of quadrature points of the cell
   */
  std::vector<double>
  get_weights(const Point<dim> &face_center, const bool is_inner) const
  {
    if (dim == 1)
      return std::vector<double>{1.0};

    // Mesh alignment type and cell index
    const auto [type, id_in_plane, _] = get_config(face_center, is_inner);
    // Angle variation within each cell
    const double delta_0 = 2 * numbers::PI / n_subdivisions[0];
    double       delta_1 = 1.0;

    if (dim == 3)
      delta_1 = radius[1] / n_subdivisions[1];

    if (type == 0) // aligned
      {
        std::vector<double> weights;
        weights.reserve(n_quadrature_points);

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          weights.emplace_back(radius[0] * quadrature.weight(q) * delta_0 *
                               delta_1);

        return weights;
      }
    else // inside/outside
      {
        double rad_0, rad_1, rad_2;

        double rot_min =
          rotation_angle - std::floor(rotation_angle / delta_0) * delta_0;

        if (type == 2) // outside
          {
            rad_0 = id_in_plane * delta_0;
            rad_1 = id_in_plane * delta_0 + rot_min;
            rad_2 = (id_in_plane + 1) * delta_0;
          }
        else // inside
          {
            rad_0 = id_in_plane * delta_0 + rot_min;
            rad_1 = (id_in_plane + 1) * delta_0;
            rad_2 = (id_in_plane + 1) * delta_0 + rot_min;
          }

        std::vector<double> weights;
        weights.reserve(2 * n_quadrature_points);

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          weights.emplace_back(radius[0] * quadrature.weight(q) *
                               (rad_1 - rad_0) * delta_1);

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          weights.emplace_back(radius[0] * quadrature.weight(q) *
                               (rad_2 - rad_1) * delta_1);

        return weights;
      }
  }

  /**
   * @brief Returns the normal vector for the quadrature points
   *
   * @param[in] face_center Face center
   *
   * @return result Normal vectors of the cell quadrature points
   */
  std::vector<Tensor<1, dim, double>>
  get_normals(const Point<dim> &face_center, const bool is_inner) const
  {
    // Coordinates of cell quadrature points
    const auto points = get_points(face_center, is_inner);

    std::vector<Tensor<1, dim, double>> result;

    result.reserve(points.size());

    for (const auto &point : points)
      result.emplace_back(get_normal(point));

    return result;
  }

  /// Number of cells at the interface between inner and outer domains
  std::vector<unsigned int> n_subdivisions;
  /// Vector containing the radius at the mortar interface and the domain length
  /// in the direction of the rotation axis
  std::vector<double> radius;

protected:
  /**
   * @brief Returns the mesh alignment type and cell index
   *
   * @param[in] face_center Face center
   *
   * @return type Cell configuration type at the interface
   * type = 0: mesh aligned
   * type = 1: mesh not aligned, inner domain (allows rotation)
   * type = 2: mesh not aligned, outer domain (fixed)
   * @return id_in_plane Index of the cell in which lies the rotated cell center
   * @return id_out_plane Second index of the cell in which lies the rotated cell center.
   *
   * Note that the id_out_plane corresponds to indexes along the rotation axis,
   * and it is necessary only for 3D problems.
   */
  std::tuple<unsigned int, unsigned int, unsigned int>
  get_config(const Point<dim> &face_center, const bool is_inner) const
  {
    const auto angle_cell_center = to_1D(face_center);

    // Angular variation in each cell
    const double delta_0 = 2 * numbers::PI / n_subdivisions[0];
    // Minimum rotation angle
    const double rot_min =
      rotation_angle - std::floor(rotation_angle / delta_0) * delta_0;

    AssertThrow(rot_min <= delta_0, ExcInternalError());

    // Point position in the cell
    const double segment = (angle_cell_center - delta_0 / 2) / delta_0;
    // Point position after rotation
    const double segment_rot =
      (angle_cell_center - delta_0 / 2 - rot_min) / delta_0;
    // Cell index in the direction of the rotation axis
    unsigned int id_out_plane = 0;

    if constexpr (dim == 3)
      {
        const double delta_1 = radius[1] / n_subdivisions[1];
        id_out_plane         = static_cast<unsigned int>(
          std::round((face_center[2] - delta_1 / 2) / delta_1));
      }

    if (this->is_mesh_aligned())
      {
        // Case 1: mesh is aligned
        return {0,
                static_cast<unsigned int>(std::round(segment)),
                id_out_plane};
      }
    else
      {
        // Case 2: mesh is not aligned
        if (!is_inner)
          // outer (fixed) domain
          return {2,
                  static_cast<unsigned int>(std::round(segment)),
                  id_out_plane};
        else
          // inner (rotated) domain
          return {1,
                  (static_cast<unsigned int>(std::round(segment_rot)) +
                   2 * n_subdivisions[0]) %
                    (2 * n_subdivisions[0]),
                  id_out_plane};
      }
  }


  /**
   * @brief Convert radiant to quadrature point in real space.
   */
  virtual Point<dim>
  from_1D(const double radiant) const = 0;

  /**
   * @brief Convert quadrature point in real space to radiant.
   */
  virtual double
  to_1D(const Point<dim> &point) const = 0;

  /**
   * @brief Return the normal for a given quadrature point.
   */
  virtual Tensor<1, dim, double>
  get_normal(const Point<dim> &point) const = 0;

  /// Mortar quadrature
  Quadrature<std::max(1, dim - 1)> quadrature;
  /// Number of quadrature points per cell
  const unsigned int n_quadrature_points;
  /// Rotation angle for the inner domain
  const double rotation_angle;
};


/*-------------- Mortar manager circle -------------------------------*/
template <int dim>
class MortarManagerCircle : public MortarManagerBase<dim>
{
public:
  template <int dim2>
  MortarManagerCircle(unsigned int            n_subdivisions,
                      double                  radius,
                      const Quadrature<dim2> &quadrature,
                      const double            rotation_angle);

  template <int dim2>
  MortarManagerCircle(std::vector<unsigned int> n_subdivisions,
                      std::vector<double>       radius,
                      const Quadrature<dim2>   &quadrature,
                      const double              rotation_angle);

protected:
  Point<dim>
  from_1D(const double rad) const override;

  double
  to_1D(const Point<dim> &point) const override;

  Tensor<1, dim, double>
  get_normal(const Point<dim> &point) const override;
};

template <int dim>
template <int dim2>
MortarManagerCircle<dim>::MortarManagerCircle(
  unsigned int            n_subdivisions,
  double                  radius,
  const Quadrature<dim2> &quadrature,
  const double            rotation_angle)
  : MortarManagerBase<dim>(n_subdivisions, radius, quadrature, rotation_angle)
{}

template <int dim>
template <int dim2>
MortarManagerCircle<dim>::MortarManagerCircle(
  std::vector<unsigned int> n_subdivisions,
  std::vector<double>       radius,
  const Quadrature<dim2>   &quadrature,
  const double              rotation_angle)
  : MortarManagerBase<dim>(n_subdivisions, radius, quadrature, rotation_angle)
{}

template <int dim>
Point<dim>
MortarManagerCircle<dim>::from_1D(const double radiant) const
{
  Point<dim> point;

  point[0] = this->radius[0] * std::cos(radiant);
  point[1] = this->radius[0] * std::sin(radiant);

  return point;
}

template <int dim>
double
MortarManagerCircle<dim>::to_1D(const Point<dim> &point) const
{
  return std::fmod(std::atan2(point[1], point[0]) + 2 * numbers::PI,
                   2 * numbers::PI);
}

template <int dim>
Tensor<1, dim, double>
MortarManagerCircle<dim>::get_normal(const Point<dim> &point_in) const
{
  Point<dim> point = point_in;

  if (dim == 3)
    point[2] = 0.0;

  return point / point.norm();
}

template <int dim>
template <int dim2>
MortarManagerBase<dim>::MortarManagerBase(unsigned int n_subdivisions,
                                          double       radius,
                                          const Quadrature<dim2> &quadrature_in,
                                          const double rotation_angle)
  : MortarManagerBase(std::vector<unsigned int>{n_subdivisions, 1},
                      std::vector<double>{radius, 1.0},
                      quadrature_in,
                      rotation_angle)
{}


template <int dim>
template <int dim2>
MortarManagerBase<dim>::MortarManagerBase(
  const std::vector<unsigned int> &n_subdivisions,
  const std::vector<double>       &radius,
  const Quadrature<dim2>          &quadrature_in,
  const double                     rotation_angle)
  : n_subdivisions(n_subdivisions)
  , radius(radius)
  , quadrature(quadrature_in.get_tensor_basis()[0])
  , n_quadrature_points(quadrature.size())
  , rotation_angle(rotation_angle)
{}
