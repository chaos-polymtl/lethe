// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

struct Settings
{
  bool
  parse(const std::string &prm_filename);

  int          dimension;
  unsigned int max_refinement_level;
  unsigned int min_refinement_level;
  unsigned int initial_refinement;
  unsigned int initial_refinement_steps;

  double max_reinitialization_distance;
  double tanh_thickness;

  double time_step;
  double time_end;

  std::string output_name;
  std::string output_path;
};

bool
Settings::parse(const std::string &prm_filename)
{
  ParameterHandler prm;
  prm.declare_entry("dim",
                    "2",
                    Patterns::Integer(),
                    "The problem dimension <2|3>");
  prm.declare_entry("max refinement level",
                    "6",
                    Patterns::Integer(),
                    "Maximum number of refinment level");
  prm.declare_entry("min refinement level",
                    "4",
                    Patterns::Integer(),
                    "Minimum number of refinment level");
  prm.declare_entry("initial refinement level",
                    "1",
                    Patterns::Integer(),
                    "Initial global refinement level");
  prm.declare_entry("initial refinement steps",
                    "1",
                    Patterns::Integer(),
                    "Numeber of initial adaptive refinement steps");
  prm.declare_entry("max reinitialization distance",
                    "1.",
                    Patterns::Double(),
                    "Maximum reinitialization distance value");
  prm.declare_entry("tanh thickness",
                    "1.",
                    Patterns::Double(),
                    "Interface thickness for the tanh transformation");
  prm.declare_entry("time step", "1.", Patterns::Double(), "Time step value");
  prm.declare_entry("time end",
                    "1",
                    Patterns::Double(),
                    "Simulation time value");
  prm.declare_entry("output name",
                    "solution",
                    Patterns::FileName(),
                    "Name for vtu files");
  prm.declare_entry("output path",
                    "./",
                    Patterns::FileName(),
                    "Path for vtu output files");

  if (prm_filename.size() == 0)
    {
      std::cout
        << "****  Error: No input file provided!\n"
        << "****  Error: Call this program as './vof_advection input.prm\n"
        << '\n'
        << "****  You may want to use one of the input files in this\n"
        << "****  directory, or use the following default values\n"
        << "****  to create an input file:\n";
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
#if DEAL_II_VERSION_GTE(9, 7, 0)
        prm.print_parameters(std::cout, ParameterHandler::DefaultStyle);
#else
        prm.print_parameters(std::cout, ParameterHandler::Text);
#endif
      return false;
    }

  try
    {
      prm.parse_input(prm_filename);
    }
  catch (std::exception &e)
    {
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        std::cerr << e.what() << std::endl;
      return false;
    }

  this->dimension                = prm.get_integer("dim");
  this->max_refinement_level     = prm.get_integer("max refinement level");
  this->min_refinement_level     = prm.get_integer("min refinement level");
  this->initial_refinement       = prm.get_integer("initial refinement level");
  this->initial_refinement_steps = prm.get_integer("initial refinement steps");

  this->max_reinitialization_distance =
    prm.get_double("max reinitialization distance");
  this->tanh_thickness = prm.get_double("tanh thickness");

  this->time_step   = prm.get_double("time step");
  this->time_end    = prm.get_double("time end");
  this->output_name = prm.get("output name");
  this->output_path = prm.get("output path");

  return true;
}

template <typename T>
T
sgn(T val)
{
  return (static_cast<T>(0) < val) - (val < static_cast<T>(0));
}

namespace PrototypeGridTools
{
  // Already implemented in LetheGridTools
  template <int dim>
  inline double
  compute_point_2_interface_min_distance(
    const std::vector<Point<dim>> &triangle,
    const Point<dim>              &point)
  {
    double            D;
    const Point<dim> &point_0 = triangle[0];
    const Point<dim> &point_1 = triangle[1];

    if constexpr (dim == 3)
      {
        const Point<dim> &point_2 = triangle[2];

        Tensor<1, dim> vector_to_plane;
        Point<dim>     pt_in_triangle;

        vector_to_plane = point_0 - point;

        const Tensor<1, dim> e_0 = point_1 - point_0;
        const Tensor<1, dim> e_1 = point_2 - point_0;

        const double a = e_0.norm_square();
        const double b = scalar_product(e_0, e_1);
        const double c = e_1.norm_square();
        const double d = scalar_product(e_0, vector_to_plane);
        const double e = scalar_product(e_1, vector_to_plane);

        const double det = a * c - b * b;

        double s = b * e - c * d;
        double t = b * d - a * e;

        if (s + t <= det)
          {
            if (s < 0)
              {
                if (t < 0)
                  {
                    // Region 4
                    if (d < 0)
                      {
                        t = 0;
                        if (-d >= a)
                          s = 1;
                        else
                          s = -d / a;
                      }
                    else
                      {
                        s = 0;
                        if (e >= 0)
                          t = 0;
                        else if (-e >= c)
                          t = 1;
                        else
                          t = e / c;
                      }
                  }
                else
                  {
                    // Region 3
                    s = 0;
                    if (e >= 0)
                      t = 0;
                    else if (-e >= c)
                      t = 1;
                    else
                      t = -e / c;
                  }
              }
            else if (t < 0)
              {
                // Region 5
                t = 0;
                if (d >= 0)
                  s = 0;
                else if (-d >= a)
                  s = 1;
                else
                  s = -d / a;
              }
            else
              {
                // Region 0
                const double inv_det = 1. / det;
                s *= inv_det;
                t *= inv_det;
              }
          }
        else
          {
            if (s < 0)
              {
                // Region 2
                const double tmp0 = b + d;
                const double tmp1 = c + e;
                if (tmp1 > tmp0)
                  {
                    const double numer = tmp1 - tmp0;
                    const double denom = a - 2 * b + c;
                    if (numer >= denom)
                      s = 1;
                    else
                      s = numer / denom;

                    t = 1 - s;
                  }
                else
                  {
                    s = 0;
                    if (tmp1 <= 0)
                      t = 1;
                    else if (e >= 0)
                      t = 0;
                    else
                      t = -e / c;
                  }
              }
            else if (t < 0)
              {
                // Region 6
                const double tmp0 = b + e;
                const double tmp1 = a + d;
                if (tmp1 > tmp0)
                  {
                    const double numer = tmp1 - tmp0;
                    const double denom = a - 2 * b + c;
                    if (numer >= denom)
                      t = 1;
                    else
                      t = numer / denom;
                    s = 1 - t;
                  }
                else
                  {
                    t = 0;
                    if (tmp1 <= 0)
                      s = 1;
                    else if (d >= 0)
                      s = 0;
                    else
                      s = -d / a;
                  }
              }
            else
              {
                // Region 1
                const double numer = (c + e) - (b + d);
                if (numer <= 0)
                  s = 0;
                else
                  {
                    const double denom = a - 2 * b + c;
                    if (numer >= denom)
                      s = 1;
                    else
                      s = numer / denom;
                  }
                t = 1 - s;
              }
          }

        pt_in_triangle = point_0 + s * e_0 + t * e_1;

        D = pt_in_triangle.distance(point);
      }

    if constexpr (dim == 2)
      {
        const Tensor<1, dim> d = point_1 - point_0;

        const double t_bar = d * (point - point_0) / (d.norm() * d.norm());

        if (t_bar <= 0.0)
          {
            const Tensor<1, dim> point_minus_p0 = point - point_0;
            D                                   = point_minus_p0.norm();
          }
        else if (t_bar >= 1.0)
          {
            const Tensor<1, dim> point_minus_p1 = point - point_1;
            D                                   = point_minus_p1.norm();
          }
        else
          {
            const Tensor<1, dim> projection = point - (point_0 + t_bar * d);
            D                               = projection.norm();
          }
      }

    return D;
  }
} // namespace PrototypeGridTools

template <int dim>
class Visualization : public dealii::DataOutInterface<0, dim>
{
public:
  Visualization();

  void
  build_patches(
    const std::map<types::global_cell_index, std::vector<Point<dim>>>
      &interface_reconstruction_vertices);

private:
  /**
   * Implementation of the corresponding function of the base class.
   */
  const std::vector<DataOutBase::Patch<0, dim>> &
  get_patches() const override;

  std::vector<std::string>
  get_dataset_names() const override;

  /**
   * Output information that is filled by build_patches() and
   * written by the write function of the base class.
   */
  std::vector<DataOutBase::Patch<0, dim>> patches;

  /**
   * A list of field names for all data components stored in patches.
   */
  std::vector<std::string> dataset_names;
};

template <int dim>
Visualization<dim>::Visualization()
{}

template <int dim>
void
Visualization<dim>::build_patches(
  const std::map<types::global_cell_index, std::vector<Point<dim>>>
    &interface_reconstruction_vertices)
{
  for (auto const &cell : interface_reconstruction_vertices)
    {
      std::vector<Point<dim>> vertices = cell.second;
      for (const Point<dim> &vertex : vertices)
        {
          DataOutBase::Patch<0, dim> temp;
          temp.vertices[0] = vertex;
          patches.push_back(temp);
        }
    }
}

template <int dim>
const std::vector<DataOutBase::Patch<0, dim>> &
Visualization<dim>::get_patches() const
{
  return patches;
}

template <int dim>
std::vector<std::string>
Visualization<dim>::get_dataset_names() const
{
  return dataset_names;
}

namespace InterfaceTools
{

  /**
   * @brief Scalar function defined by the DOF values of a single cell. Based on the CellWiseFunction and RefSpaceFEFieldFunction of dealii.
   */
  template <int dim, typename VectorType = Vector<double>>
  class LocalCellWiseFunction : public Function<dim>
  {
  public:
    LocalCellWiseFunction(unsigned int fe_degree);

    void
    set_active_cell(const VectorType &in_local_dof_values);

    double
    value(const Point<dim>  &point,
          const unsigned int component = 0) const override;

    Tensor<1, dim>
    gradient(const Point<dim>  &point,
             const unsigned int component = 0) const override;

    SymmetricTensor<2, dim>
    hessian(const Point<dim>  &point,
            const unsigned int component = 0) const override;

  private:
    FE_Q<dim> fe;

    unsigned int n_local_dofs;

    Vector<typename VectorType::value_type> local_dof_values;
  };

  template <int dim, typename VectorType>
  LocalCellWiseFunction<dim, VectorType>::LocalCellWiseFunction(
    unsigned int fe_degree)
    : fe(fe_degree)
  {
    n_local_dofs = fe.dofs_per_cell;
  }

  template <int dim, typename VectorType>
  void
  LocalCellWiseFunction<dim, VectorType>::set_active_cell(
    const VectorType &in_local_dof_values)
  {
    local_dof_values = in_local_dof_values;
  }

  template <int dim, typename VectorType>
  double
  LocalCellWiseFunction<dim, VectorType>::value(
    const Point<dim>  &point,
    const unsigned int component) const
  {
    double value = 0;
    for (unsigned int i = 0; i < n_local_dofs; ++i)
      value +=
        local_dof_values[i] * fe.shape_value_component(i, point, component);

    return value;
  }

  template <int dim, typename VectorType>
  Tensor<1, dim>
  LocalCellWiseFunction<dim, VectorType>::gradient(
    const Point<dim>  &point,
    const unsigned int component) const
  {
    Tensor<1, dim> gradient;
    for (unsigned int i = 0; i < n_local_dofs; ++i)
      gradient +=
        local_dof_values[i] * fe.shape_grad_component(i, point, component);

    return gradient;
  }

  template <int dim, typename VectorType>
  SymmetricTensor<2, dim>
  LocalCellWiseFunction<dim, VectorType>::hessian(
    const Point<dim>  &point,
    const unsigned int component) const
  {
    Tensor<2, dim> hessian;
    for (unsigned int i = 0; i < n_local_dofs; ++i)
      hessian +=
        local_dof_values[i] * fe.shape_grad_grad_component(i, point, component);

    return symmetrize(hessian);
  }

  /**
   * @brief
   * Compute the volume enclosed by the 0 level of a level set field
   * inside a cell.
   *
   * @param[in] fe_point_evaluation FePointEvaluation
   *
   * @param[in] cell Cell for which the volume is computed
   *
   * @param[in] cell_dof_values cell DOFs value of the level set field
   *
   * @param[in] corr correction to apply to the DOF values (constant for all
   * DOFs)
   *
   * @param n_quad_points[in] number of quadrature points for the volume integration
   * faces
   *
   * @return cell-wise volume
   */
  template <int dim>
  double
  compute_cell_wise_volume(
    FEPointEvaluation<1, dim>                            &fe_point_evaluation,
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    Vector<double>     cell_dof_level_set_values,
    const double       corr,
    const unsigned int n_quad_points)
  {
    const unsigned int n_dofs = cell_dof_level_set_values.size();

    // Initialize required variables to compute local volume
    const BoundingBox<dim>     unit_box = create_unit_bounding_box<dim>();
    LocalCellWiseFunction<dim> signed_distance_function =
      LocalCellWiseFunction<dim>(cell->get_fe().degree);

    hp::QCollection<1> q_collection;
    q_collection.push_back(QGauss<1>(n_quad_points));

    NonMatching::QuadratureGenerator<dim> quadrature_generator =
      NonMatching::QuadratureGenerator<dim>(q_collection);

    for (unsigned int j = 0; j < n_dofs; j++)
      {
        cell_dof_level_set_values[j] += corr;
      }

    signed_distance_function.set_active_cell(cell_dof_level_set_values);
    quadrature_generator.generate(signed_distance_function, unit_box);

    const Quadrature<dim> inside_quadrature =
      quadrature_generator.get_inside_quadrature();

    if (inside_quadrature.size() == 0)
      return 0.0;

    fe_point_evaluation.reinit(cell, inside_quadrature.get_points());

    double inside_cell_volume = 0.0;
    for (unsigned int q = 0; q < inside_quadrature.size(); q++)
      {
        /* Compute the volume int 1*JxW*dOmega. FEPointEvaluation.JxW() does not
        return the right thing.*/
        inside_cell_volume += fe_point_evaluation.jacobian(q).determinant() *
                              inside_quadrature.weight(q);
      }

    return inside_cell_volume;
  }

  /**
   * @brief
   * Compute the volume enclosed by the 0 level of a level set field
   * in the domain.
   *
   * @param[in] mapping Mapping of the domain
   *
   * @param[in] dof_handler DofHandler associated to the triangulation on which
   * the volume is computed
   *
   * @param[in] level_set_vector Level-set vector
   *
   * @param[in] mpi_communicator MPI communicator
   *
   * @return Volume enclosed by the 0 level
   */
  template <int dim, typename VectorType>
  double
  compute_volume(const Mapping<dim>       &mapping,
                 const DoFHandler<dim>    &dof_handler,
                 const FiniteElement<dim> &fe,
                 const VectorType         &level_set_vector,
                 const MPI_Comm           &mpi_communicator)
  {
    FEPointEvaluation<1, dim> fe_point_evaluation(
      mapping, fe, update_jacobians | update_JxW_values);

    double volume = 0.0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            const unsigned int n_dofs_per_cell =
              cell->get_fe().n_dofs_per_cell();
            Vector<double> cell_dof_level_set_values(n_dofs_per_cell);

            cell->get_dof_values(level_set_vector,
                                 cell_dof_level_set_values.begin(),
                                 cell_dof_level_set_values.end());

            volume += compute_cell_wise_volume(fe_point_evaluation,
                                               cell,
                                               cell_dof_level_set_values,
                                               0.0,
                                               cell->get_fe().degree + 1);
          }
      }
    volume = Utilities::MPI::sum(volume, mpi_communicator);

    return volume;
  }

  /**
   * @brief
   * Reconstruct the interface defined by the 0 level of a level set field
   * in the domain.
   *
   * @param[in] mapping Mapping of the domain
   *
   * @param[in] dof_handler DofHandler associated to the triangulation for which
   * the interface is reconstructed
   *
   * @param[in] fe Finite element
   *
   * @param[in] level_set_vector Level-set vector
   *
   * @param[in,out] interface_reconstruction_vertices Cell-wise map of the
   * reconstructed surface vertices. The map contains vectors storing the
   * vertices of the reconstructed surface for each intersected volume cell
   * (dim).
   *
   * @param[in,out] interface_reconstruction_cells Cell-wise map of the
   * reconstructed surface cells. The map contains vectors storing the cell
   * (dim-1) of the reconstructed surface for each intersected volume cell
   * (dim).
   *
   * @param[in,out] intersected_dofs Set of DOFs that belong to intersected
   * volume cell (dim).
   *
   */
  template <int dim, typename VectorType>
  void
  reconstruct_interface(
    const Mapping<dim>       &mapping,
    const DoFHandler<dim>    &dof_handler,
    const FiniteElement<dim> &fe,
    const VectorType         &level_set_vector,
    std::map<types::global_cell_index, std::vector<Point<dim>>>
      &interface_reconstruction_vertices,
    std::map<types::global_cell_index, std::vector<CellData<dim - 1>>>
                                      &interface_reconstruction_cells,
    std::set<types::global_dof_index> &intersected_dofs)
  {
    // Warning: for fe.degree=2, at least 2 subdivisions should be used for the
    // marching-cube algorithm
    GridTools::MarchingCubeAlgorithm<dim, VectorType> marching_cube(mapping,
                                                                    fe,
                                                                    1,
                                                                    1e-10);
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            std::vector<Point<dim>>        surface_vertices;
            std::vector<CellData<dim - 1>> surface_cells;

            marching_cube.process_cell(
              cell, level_set_vector, 0.0, surface_vertices, surface_cells);

            // If the cell is intersected, reconstruct the interface in it
            if (surface_vertices.size() != 0)
              {
                const unsigned int cell_index =
                  cell->global_active_cell_index();

                // Store the interface reconstruction vertices and cells
                interface_reconstruction_vertices[cell_index] =
                  surface_vertices;
                interface_reconstruction_cells[cell_index] = surface_cells;

                // Store the dofs of the intersected volume cell
                std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
                cell->get_dof_indices(dof_indices);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    intersected_dofs.insert(dof_indices[i]);
                  }
              }
          }
      }
  }

  template <int dim>
  class SignedDistanceSolver
  {
  public:
    SignedDistanceSolver(const parallel::DistributedTriangulationBase<dim>
                                         &background_triangulation,
                         const FE_Q<dim> &background_fe,
                         const double     p_max_distance)
      : dof_handler(background_triangulation)
      , fe(background_fe)
      , mapping(fe.degree)
      , max_distance(p_max_distance)
      , pcout(std::cout,
              (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
    {
      set_dof_opposite_faces_map();
    }

    void
    setup_dofs(const MPI_Comm &mpi_communicator);
    void
    set_level_set_from_background_mesh(
      const DoFHandler<dim>               &background_dof_handler,
      const TrilinosWrappers::MPI::Vector &background_level_set_vector,
      const MPI_Comm                      &mpi_communicator);
    void
    solve(const MPI_Comm &mpi_communicator);

    TrilinosWrappers::MPI::Vector
    get_level_set(const MPI_Comm &mpi_communicator);

    void
    output_interface_reconstruction(const std::string  output_name,
                                    const std::string  output_path,
                                    const unsigned int it) const;


    DoFHandler<dim> dof_handler;

  private:
    using VectorType = TrilinosWrappers::MPI::Vector;
    using MatrixType = TrilinosWrappers::SparseMatrix;

    void
    zero_out_ghost_values();
    void
    update_ghost_values();
    void
    exchange_distance();

    void
    initialize_local_distance();

    void
    compute_first_neighbors_distance();
    void
    compute_second_neighbors_distance(const MPI_Comm &mpi_communicator);

    void
    compute_signed_distance_from_distance();

    void
    compute_cell_wise_volume_correction(const MPI_Comm &mpi_communicator);

    void
    compute_volume_correction_L2_projection(Vector<double> &eta_cell,
                                            const MPI_Comm &mpi_communicator);

    void
    conserve_global_volume(const MPI_Comm &mpi_communicator);

    inline void
    get_dof_opposite_faces(unsigned int               local_dof_id,
                           std::vector<unsigned int> &local_opposite_faces);

    inline unsigned int
    get_n_opposite_faces_per_dof(unsigned int local_dof_id);

    inline void
    set_dof_opposite_faces_map();

    inline void
    get_face_transformation_jacobian(
      const DerivativeForm<1, dim, dim> &cell_transformation_jac,
      const unsigned int                 local_face_id,
      DerivativeForm<1, dim - 1, dim>   &face_transformation_jac);

    inline Point<dim>
    transform_ref_face_point_to_ref_cell(const Point<dim - 1> &x_ref_face,
                                         const unsigned int    local_face_id);

    inline void
    compute_residual(const Tensor<1, dim>                  &x_n_to_x_I_real,
                     const Tensor<1, dim>                  &distance_gradient,
                     const DerivativeForm<1, dim - 1, dim> &transformation_jac,
                     Tensor<1, dim - 1>                    &residual_ref);

    inline std::vector<Point<dim>>
    compute_numerical_jacobian_stencil(const Point<dim>  &x_ref,
                                       const unsigned int local_face_id,
                                       const double       perturbation);

    inline Tensor<1, dim>
    transform_ref_face_correction_to_ref_cell(
      const Vector<double> &correction_ref_face,
      const unsigned int    local_face_id);

    inline void
    compute_numerical_jacobian(
      const std::vector<Point<dim>>     &stencil_real,
      const Point<dim>                  &x_I_real,
      const std::vector<Tensor<1, dim>> &distance_gradients,
      const std::vector<DerivativeForm<1, dim - 1, dim>>
                               &transformation_jacobians,
      const double              perturbation,
      LAPACKFullMatrix<double> &jacobian_matrix);

    inline double
    compute_distance(const Tensor<1, dim> &x_n_to_x_I_real,
                     const double          distance);


    FE_Q<dim>     fe;
    MappingQ<dim> mapping;

    const double max_distance;

    ConditionalOStream pcout;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    IndexSet locally_active_dofs;

    TrilinosWrappers::MPI::Vector level_set;

    LinearAlgebra::distributed::Vector<double> signed_distance;
    LinearAlgebra::distributed::Vector<double> signed_distance_with_ghost;
    LinearAlgebra::distributed::Vector<double> distance;
    LinearAlgebra::distributed::Vector<double> distance_with_ghost;
    LinearAlgebra::distributed::Vector<double> volume_correction;

    AffineConstraints<double> constraints;
    MatrixType                system_matrix_volume_correction;
    VectorType                system_rhs_volume_correction;

    std::map<types::global_cell_index, std::vector<Point<dim>>>
      interface_reconstruction_vertices;
    std::map<types::global_cell_index, std::vector<CellData<dim - 1>>>
      interface_reconstruction_cells;

    std::set<types::global_dof_index>                 intersected_dofs;
    std::map<unsigned int, std::vector<unsigned int>> dof_opposite_faces_map;
  };

  template <int dim>
  void
  SignedDistanceSolver<dim>::setup_dofs(const MPI_Comm &mpi_communicator)
  {
    dof_handler.distribute_dofs(fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    locally_active_dofs = DoFTools::extract_locally_active_dofs(dof_handler);

    level_set.reinit(locally_owned_dofs,
                     locally_relevant_dofs,
                     mpi_communicator);

    signed_distance.reinit(locally_owned_dofs,
                           locally_active_dofs,
                           mpi_communicator);
    signed_distance_with_ghost.reinit(locally_owned_dofs,
                                      locally_active_dofs,
                                      mpi_communicator);

    distance.reinit(locally_owned_dofs, locally_active_dofs, mpi_communicator);
    distance_with_ghost.reinit(locally_owned_dofs,
                               locally_active_dofs,
                               mpi_communicator);

    volume_correction.reinit(locally_owned_dofs,
                             locally_active_dofs,
                             mpi_communicator);

    constraints.clear();
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.close();

    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs =*/false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               dof_handler.locally_owned_dofs(),
                                               mpi_communicator,
                                               locally_relevant_dofs);
    system_matrix_volume_correction.reinit(locally_owned_dofs,
                                           locally_owned_dofs,
                                           dsp,
                                           mpi_communicator);
    system_rhs_volume_correction.reinit(locally_owned_dofs, mpi_communicator);
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::set_level_set_from_background_mesh(
    const DoFHandler<dim>               &background_dof_handler,
    const TrilinosWrappers::MPI::Vector &background_level_set_vector,
    const MPI_Comm                      &mpi_communicator)
  {
    TrilinosWrappers::MPI::Vector tmp_local_level_set(this->locally_owned_dofs,
                                                      mpi_communicator);

    VectorTools::interpolate_to_different_mesh(background_dof_handler,
                                               background_level_set_vector,
                                               dof_handler,
                                               tmp_local_level_set);

    level_set = tmp_local_level_set;
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::solve(const MPI_Comm &mpi_communicator)
  {
    // Gain the writing right.
    zero_out_ghost_values();

    // Clear maps and sets
    interface_reconstruction_vertices.clear();
    interface_reconstruction_cells.clear();
    intersected_dofs.clear();

    // Initialize local distance vectors.
    initialize_local_distance();

    // Identify intersected cells and compute the interface reconstruction.
    reconstruct_interface(mapping,
                          dof_handler,
                          fe,
                          level_set,
                          interface_reconstruction_vertices,
                          interface_reconstruction_cells,
                          intersected_dofs);

    /* Compute the distance for the dofs of the intersected cells (the ones in
    the intersected_dofs set). They correspond to the first neighbor dofs.*/
    compute_first_neighbors_distance();

    /* Compute signed distance from distance (only first neighbors have an
    updated value)*/
    compute_signed_distance_from_distance();

    // Conserve local and global volume
    compute_cell_wise_volume_correction(mpi_communicator);
    conserve_global_volume(mpi_communicator);

    /* Compute the distance for the dofs of the rest of the mesh. They
    correspond to the second neighbors dofs. */
    compute_second_neighbors_distance(mpi_communicator);

    // Compute signed distance from distance (all DOFs have updated value)
    compute_signed_distance_from_distance();

    // Update ghost values to regain reading ability.
    update_ghost_values();
  }

  template <int dim>
  TrilinosWrappers::MPI::Vector
  SignedDistanceSolver<dim>::get_level_set(const MPI_Comm &mpi_communicator)
  {
    TrilinosWrappers::MPI::Vector tmp_local_level_set(this->locally_owned_dofs,
                                                      mpi_communicator);

    for (auto p : this->locally_owned_dofs)
      {
        tmp_local_level_set(p) = signed_distance(p);
      }

    level_set = tmp_local_level_set;

    return level_set;
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::zero_out_ghost_values()
  {
    /* To have the right to write in a LinearAlgebra::distributed::Vector, we
    have to zero out the ghost values.*/
    signed_distance.zero_out_ghost_values();
    signed_distance_with_ghost.zero_out_ghost_values();
    distance.zero_out_ghost_values();
    distance_with_ghost.zero_out_ghost_values();
    volume_correction.zero_out_ghost_values();
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::update_ghost_values()
  {
    /* To have the right to read a LinearAlgebra::distributed::Vector, we have
    to update the ghost values.*/
    signed_distance.update_ghost_values();
    signed_distance_with_ghost.update_ghost_values();
    distance.update_ghost_values();
    distance_with_ghost.update_ghost_values();
    volume_correction.update_ghost_values();
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::exchange_distance()
  {
    // Exchange and "select" min value between the processes
    distance.compress(VectorOperation::min);

    // Update local ghost (distance becomes read only)
    distance.update_ghost_values();

    /* Copy distance to distance_with_ghost to keep the knowledge of local ghost
    values and to have a read only version of the vector*/
    distance_with_ghost = distance;

    /* Zero out ghost DOFs to regain write functionalities in distance (it
    becomes write only, that is why we need distance_with_ghost - to read the
    ghost values in it).*/
    distance.zero_out_ghost_values();

    /* Copy the ghost values back in distance (zero_out_ghost_values() puts
    zeros in ghost DOFs)*/
    for (auto p : this->locally_active_dofs)
      {
        /* We need to have the ghost values in distance for future
        compress(VectorOperation::min) operation*/
        distance(p) = distance_with_ghost(p);
      }
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::initialize_local_distance()
  {
    /* Initialization of the active dofs to the max distance we want to
     redistanciate. It requires a loop on the active dofs to initialize also
     the distance value of the ghost dofs. Otherwise, they are set to 0.0 by
     default. */
    for (auto p : this->locally_active_dofs)
      {
        distance(p)            = max_distance;
        distance_with_ghost(p) = max_distance;

        const double level_set_value  = level_set(p);
        signed_distance(p)            = max_distance * sgn(level_set_value);
        signed_distance_with_ghost(p) = max_distance * sgn(level_set_value);
      }
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::compute_first_neighbors_distance()
  {
    /* The signed distance for the first neighbors (the Dofs belonging to the
     * cells intersected by the reconstructed interface. This is a brute force
     * distance computation, meaning the distance is computed geometrically, as
     * presented by Ausas (2010). */

    // DoF coordinates
    std::map<types::global_dof_index, Point<dim>> dof_support_points =
      DoFTools::map_dofs_to_support_points(mapping, dof_handler);

    // Loop over the intersected cells (volume cells)
    for (auto &intersected_cell : interface_reconstruction_cells)
      {
        const unsigned int cell_index = intersected_cell.first;

        // Create interface recontruction triangulation (surface triangulation)
        // in the intersected volume cell
        std::vector<Point<dim>> surface_vertices =
          interface_reconstruction_vertices.at(cell_index);
        std::vector<CellData<dim - 1>> surface_cells = intersected_cell.second;

        Triangulation<dim - 1, dim> surface_triangulation;
        surface_triangulation.create_triangulation(surface_vertices,
                                                   surface_cells,
                                                   {});

        /* Loop over all DoFs of the volume mesh belonging to a intersected
         * volume cell. This is more expensive, but it is required to have the
         * the right signed distance approximation for the first neighbors.*/
        for (const unsigned int &intersected_dof : intersected_dofs)
          {
            const Point<dim> y = dof_support_points.at(intersected_dof);

            /* Loop over the surface cells of the interface reconstruction in
             * the volume cell. In 2D, there is only 1 surface cell (line),
             * while in 3D, it can vary from 1 to 4 or 5 (triangles), depending
             * on the marching cube algorithm.*/
            for (const auto &surface_cell :
                 surface_triangulation.active_cell_iterators())
              {
                // Store the current surface cell vertex coordinates
                unsigned int surface_cell_n_vertices =
                  surface_cell->n_vertices();
                std::vector<Point<dim>> surface_cell_vertices(
                  surface_cell_n_vertices);
                for (unsigned int p = 0; p < surface_cell_n_vertices; p++)
                  {
                    surface_cell_vertices[p] = surface_cell->vertex(p);
                  }

                // Compute the geometrical distance between the surface cell
                // (line in 2D, triangle in 3D) and the DoF
                double D =
                  PrototypeGridTools::compute_point_2_interface_min_distance(
                    surface_cell_vertices, y);

                // Select the minimum distance
                distance(intersected_dof) =
                  std::min(std::abs(distance(intersected_dof)), std::abs(D));
              }
          }
      }
    exchange_distance();
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::compute_second_neighbors_distance(
    const MPI_Comm &mpi_communicator)
  {
    /* The signed distance for the second neighbors (the cells not intersected
     * by the interface is resolved according the minimization problem
     * presented by Ausas (2010). The method looks for the point in the opposite
     * faces of each second neighbor DoFs that minimizes the distance to the
     * interface. It works in a similar manner as a marching algorithm from the
     * knowledge of the signed distance for the interface first neighbors. */

    std::vector<unsigned int> dof_opposite_faces;
    const unsigned int        dofs_per_cell = fe.n_dofs_per_cell();

    std::map<types::global_dof_index, Point<dim>> dof_support_points =
      DoFTools::map_dofs_to_support_points(mapping, dof_handler);

    Point<dim - 1> ref_face_center_point = Point<dim - 1>();
    ref_face_center_point(0)             = 0.5;
    if constexpr (dim == 3)
      ref_face_center_point(1) = 0.5;

    FEPointEvaluation<1, dim> fe_point_evaluation(
      mapping, fe, update_values | update_gradients | update_jacobians);

    /* The method is iterative, hence, we solve as long as the distance
    approximation changes for at least one dof. We use the flag change to
    track this change. */
    bool change = true;

    /* The count corresponds to how many times we iterate. In fact, it
    correspond to the number of cell layers (starting from the interface)
    that the approximation of the distance is known. */
    int count = 0;
    while (change)
      {
        pcout << "Redistanciating layer " << count << std::endl;
        change = false;

        for (const auto &cell : dof_handler.active_cell_iterators())
          {
            if (cell->is_locally_owned())
              {
                const unsigned int cell_index =
                  cell->global_active_cell_index();

                // If the cell is intersected, the distance is already computed.
                if (interface_reconstruction_vertices.find(cell_index) !=
                    interface_reconstruction_vertices.end())
                  {
                    continue;
                  }

                std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

                cell->get_dof_indices(dof_indices);
                std::vector<double> cell_dof_values(dofs_per_cell);
                cell->get_dof_values(distance_with_ghost,
                                     cell_dof_values.begin(),
                                     cell_dof_values.end());

                // Loop over the cell's Dofs
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    /* If the dof belongs to an intersected cell, the distance
                    is already computed */
                    if (intersected_dofs.find(dof_indices[i]) !=
                        intersected_dofs.end())
                      {
                        continue;
                      }

                    // Get opposite faces
                    const unsigned int n_opposite_faces_per_dof =
                      get_n_opposite_faces_per_dof(i);
                    dof_opposite_faces.resize(n_opposite_faces_per_dof);
                    get_dof_opposite_faces(i, dof_opposite_faces);

                    // Get the real coordinates of the current DoF I
                    const Point<dim> x_I_real =
                      dof_support_points.at(dof_indices[i]);

                    // Loop on opposite faces F_J
                    for (unsigned int j = 0; j < n_opposite_faces_per_dof; ++j)
                      {
                        /* The minimization problem is: Find x in the face F_J
                        (opposite to the DoF of interest I) such that:

                          |d|_{x_I} = min(phi(x) +|x_I - x|)

                        where x_I is the coord of the DoF I, phi(x) is the
                        distance (not signed) at the point x, belonging to the
                        face F_J. Here, we solve the problem in the reference
                        space (dim - 1).
                        */

                        // Initialize required variables
                        Point<dim> x_n_ref =
                          transform_ref_face_point_to_ref_cell(
                            ref_face_center_point, dof_opposite_faces[j]);
                        Point<dim> x_n_real;

                        double correction_norm = 1.0;
                        int    newton_it       = 0;

                        // Check to constrain the solution in the face F_J
                        int outside_check = 0;

                        // Solve the minimization problem with Newton method
                        // using a numerical jacobian
                        while (correction_norm > 1e-10 && outside_check < 3 &&
                               newton_it < 100)
                          {
                            /* Set stencil for numerical jacobian computation.
                             The entries of the vector are the following:
                                      4

                                 1    0    2

                                      3
                            The entry 0 is the current evaluation point. */

                            const double            perturbation = 0.01;
                            std::vector<Point<dim>> stencil_ref =
                              compute_numerical_jacobian_stencil(
                                x_n_ref, dof_opposite_faces[j], perturbation);

                            std::vector<Point<dim>> stencil_real(2 * dim - 1);
                            std::vector<Tensor<1, dim>> distance_gradients(
                              2 * dim - 1);
                            std::vector<DerivativeForm<1, dim, dim>>
                              cell_transformation_jacobians(2 * dim - 1);
                            std::vector<DerivativeForm<1, dim - 1, dim>>
                              face_transformation_jacobians(2 * dim - 1);


                            /* Prepare FEPointEvaluation to compute value and
                            gradient at the stencil points*/
                            fe_point_evaluation.reinit(cell, stencil_ref);
                            fe_point_evaluation.evaluate(
                              cell_dof_values, EvaluationFlags::gradients);

                            // Get the required values at each stencil point
                            for (unsigned int k = 0; k < 2 * dim - 1; k++)
                              {
                                stencil_real[k] =
                                  fe_point_evaluation.quadrature_point(k);
                                distance_gradients[k] =
                                  fe_point_evaluation.get_gradient(k);
                                cell_transformation_jacobians[k] =
                                  fe_point_evaluation.jacobian(k);
                                get_face_transformation_jacobian(
                                  cell_transformation_jacobians[k],
                                  dof_opposite_faces[j],
                                  face_transformation_jacobians[k]);
                              }

                            /* Compute the jacobian matrix. The Ax=b system is
                            formulated as the dim-1 system. We solve for the
                            correction in the reference face. */
                            LAPACKFullMatrix<double> jacobian_matrix(dim - 1,
                                                                     dim - 1);
                            compute_numerical_jacobian(
                              stencil_real,
                              x_I_real,
                              distance_gradients,
                              face_transformation_jacobians,
                              perturbation,
                              jacobian_matrix);

                            const Tensor<1, dim> x_n_to_x_I_real =
                              x_I_real - stencil_real[0];

                            // Compute the right hand side.
                            Tensor<1, dim - 1> residual_n;
                            compute_residual(x_n_to_x_I_real,
                                             distance_gradients[0],
                                             face_transformation_jacobians[0],
                                             residual_n);

                            // Convert the right hand side to the right format
                            // for the linear solver
                            Vector<double> residual_n_vec(dim - 1);
                            residual_n.unroll(residual_n_vec.begin(),
                                              residual_n_vec.end());
                            residual_n_vec *= -1.0;

                            jacobian_matrix.set_property(
                              LAPACKSupport::general);

                            /* Factorize and solve the matrix. The correction is
                            put back in residual_n_vec. */
                            jacobian_matrix.compute_lu_factorization();
                            jacobian_matrix.solve(residual_n_vec);

                            // Compute the norm of the correction
                            correction_norm = residual_n_vec.l2_norm();

                            /* Transform the dim-1 correction (in the reference
                            face) to dim (in the reference cell) */
                            Tensor<1, dim> correction =
                              transform_ref_face_correction_to_ref_cell(
                                residual_n_vec, dof_opposite_faces[j]);

                            /* Compute the solution (the point x_n_ref on the
                            face minimizing the distance)*/
                            Point<dim> x_n_p1_ref = stencil_ref[0] + correction;

                            /* Relax the correction if it brings us outside of
                            the cell */
                            double relaxation = 1.0;

                            /* Check if the Newton method results in a solution
                            outside the face. For example in 3D, we could have:
                                 _____________
                                |             |     solution
                                |             |    *
                                |             |
                                |             |
                                |             |
                                |_____________|

                            Each time it does, we relax the scheme to bring
                            back the estimation of the solution in the face:
                                 _____________
                                |             | relaxed solution
                                |           * |
                                |             |
                                |             |
                                |             |
                                |_____________|

                            If the solution is outside the face more than three
                            times, we constraint the solution on the right
                            boundary of the face:
                                 _____________
                                |             |         real
                                |  constraint *     * solution
                                |   solution  |
                                |             |
                                |             |
                                |_____________|

                            */

                            /* Flag indicating if the correction brings us
                            outside of the cell.*/
                            bool check = false;
                            for (unsigned int k = 0; k < dim; ++k)
                              {
                                if (x_n_p1_ref[k] > 1.0 + 1e-12 ||
                                    x_n_p1_ref[k] < 0.0 - 1e-12)
                                  {
                                    check = true;

                                    /* Set the correction to put the solution on
                                    the face boundary. Select the minimum
                                    relaxation of the all directions to ensure
                                    the solution stays inside the face.*/
                                    if (correction[k] > 1e-12)
                                      {
                                        relaxation =
                                          std::min((1.0 - x_n_ref[k]) /
                                                     (correction[k] + 1e-12),
                                                   relaxation);
                                      }
                                    else if (correction[k] < -1e-12)
                                      {
                                        relaxation =
                                          std::min((0.0 - x_n_ref[k]) /
                                                     (correction[k] + 1e-12),
                                                   relaxation);
                                      }
                                  }
                              }

                            // Increment the outside_check if the correction
                            // brought us outside the face
                            if (check)
                              outside_check += 1;

                            // Re-compute the solution with the relaxation
                            x_n_p1_ref =
                              stencil_ref[0] + relaxation * correction;

                            // Transform the solution from reference to the real
                            // cell. This could be improved to not call
                            // fe_point_evaluation.
                            std::vector<Point<dim>> x_n_p1_ref_vec = {
                              x_n_p1_ref};
                            fe_point_evaluation.reinit(cell, x_n_p1_ref_vec);
                            Point<dim> x_n_p1_real =
                              fe_point_evaluation.quadrature_point(0);

                            // Update the solution.
                            x_n_ref  = x_n_p1_ref;
                            x_n_real = x_n_p1_real;

                            newton_it += 1;
                          } // End of the Newton solver.

                        // Compute the distance approximation: distance(x_I) =
                        // distance(x_n) + |x_n - x_I|
                        const Tensor<1, dim> x_n_to_x_I_real =
                          x_I_real - x_n_real;
                        fe_point_evaluation.evaluate(cell_dof_values,
                                                     EvaluationFlags::values);

                        double distance_value_at_x_n =
                          fe_point_evaluation.get_value(0);

                        double approx_distance =
                          compute_distance(x_n_to_x_I_real,
                                           distance_value_at_x_n);

                        // If the new distance is smaller than the previous,
                        // update the value and flag the change
                        if (distance(dof_indices[i]) > (approx_distance + 1e-8))
                          {
                            change                   = true;
                            distance(dof_indices[i]) = approx_distance;
                          }
                      } // End of the loop on the opposite faces
                  }     // End of the loop on the dofs
              }
          } // End of the loop on the cells

        exchange_distance();

        // Track the change flag across the processes
        change = Utilities::MPI::logical_or(change, mpi_communicator);

        count += 1;
      } // End of the iterative while loop

    // Update the hanging node values
    constraints.distribute(distance);
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::compute_signed_distance_from_distance()
  {
    for (auto p : this->locally_active_dofs)
      {
        signed_distance(p) = distance(p) * sgn(signed_distance_with_ghost(p));
      }

    // Update local ghost (signed_distance becomes read only)
    signed_distance.update_ghost_values();

    /* Copy distance to signed_distance_with_ghost to keep the knowledge of
    local ghost values and to have a read-only version of the vector. Ghost
    values of the signed distance are needed for volume computations.*/
    signed_distance_with_ghost = signed_distance;

    /* Zero out ghost DOFs to regain write functionalities in signed_distance
    (it becomes write-only, that is why we need signed_distance_with_ghost -
    to read the ghost values in it).*/
    signed_distance.zero_out_ghost_values();

    /* Copy the ghost values back in signed_distance (zero_out_ghost_values()
    puts zeros in ghost DOFs)*/
    for (auto p : this->locally_active_dofs)
      {
        signed_distance(p) = signed_distance_with_ghost(p);
      }
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::compute_cell_wise_volume_correction(
    const MPI_Comm &mpi_communicator)
  {
    FEPointEvaluation<1, dim> fe_point_evaluation(
      mapping, fe, update_jacobians | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    /* For the approximate L2 projection of the cell-wise correction (the
    projection for a given DOF corresponds to the average of the neighbor cell
    values)*/
    // double n_cells_per_dofs_inv = 1.0 / 4.0;
    // if constexpr (dim == 3)
    //   {
    //     n_cells_per_dofs_inv = 1.0 / 8.0;
    //   }

    // Re-initialize volume_correction vector.
    volume_correction = 0.0;
    // Store eta_n for all locally owned active cells that are intersected
    Vector<double> eta_cell(dof_handler.get_triangulation().n_active_cells());

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            const unsigned int cell_index = cell->global_active_cell_index();

            // The cell is not intersected, no need to correct the volume
            if (interface_reconstruction_vertices.find(cell_index) ==
                interface_reconstruction_vertices.end())
              {
                continue;
              }

            /* We want to find a cell wise correction to apply to the cell's dof
            values of the signed_distance so that the geometric cell wise volume
            encompassed by the level 0 of the signed_distance V_K and by the
            iso-contour 0.5 of the phase fraction V_K,VOF match. This is
            required because the computed distance doesn't belong to the Q1
            approximation space.

            We solve the non-linear problem: DeltaV_K(phi* + eta_K) = V_K,VOF -
            V_K(phi* + eta_K) = 0, where phi* is the redistanciated
            signed distance, eta_K is the correction on the signed_distance that
            we are looking for. We use the secant method to do so. See Ausas et
            al. (2010) for more details.*/

            // Get the level set values
            Vector<double> cell_level_set_dof_values(dofs_per_cell);

            cell->get_dof_values(level_set,
                                 cell_level_set_dof_values.begin(),
                                 cell_level_set_dof_values.end());

            // Compute the targeted volume to correct for
            double targeted_cell_volume =
              InterfaceTools::compute_cell_wise_volume(
                fe_point_evaluation,
                cell,
                cell_level_set_dof_values,
                0.0,
                fe.degree + 1);

            // Get the signed distance values to be corrected
            Vector<double> cell_dof_values(dofs_per_cell);
            cell->get_dof_values(signed_distance_with_ghost,
                                 cell_dof_values.begin(),
                                 cell_dof_values.end());

            /* Get cell size to initialize secant method. We use it to compute
            the first derivative value in the secant method */
            double cell_size;
            if (dim == 2)
              {
                cell_size = std::sqrt(4. * cell->measure() / M_PI);
              }
            else if (dim == 3)
              {
                cell_size = std::pow(6 * cell->measure() / M_PI, 1. / 3.);
              }

            /* Secant method. The subscript nm1 (or n minus 1) stands for the
            previous secant iteration (it = n-1), the subscript n stands for
            the current iteration and the subscript np1 stands for the next
            iteration (it = n+1).*/
            double inside_cell_volume_nm1 = 0.0;
            double inside_cell_volume_n   = 0.0;

            double delta_volume_nm1 = 0.0;
            double delta_volume_n   = 0.0;

            double delta_volume_prime = 0.0;

            double eta_nm1 = 0.0;
            double eta_n   = 1e-6 * cell_size;
            double eta_np1 = 0.0;

            // Compute the volume for the first initial value (eta_nm1)
            inside_cell_volume_nm1 =
              InterfaceTools::compute_cell_wise_volume(fe_point_evaluation,
                                                       cell,
                                                       cell_dof_values,
                                                       eta_nm1,
                                                       fe.degree + 1);
            delta_volume_nm1 = targeted_cell_volume - inside_cell_volume_nm1;

            /* Store the initial volume in the cell to limit the secant method
            in some case.*/
            const double initial_inside_cell_volume = inside_cell_volume_nm1;

            /* Check if there is enough volume to correct. If not, we don't
            correct.*/
            if (inside_cell_volume_nm1 < 1e-10 * cell_size ||
                inside_cell_volume_nm1 > (cell_size - 1e-10 * cell_size))
              {
                eta_n = 0.0;
                continue;
              }

            unsigned int secant_it     = 0;
            double       secant_update = 1.0;
            while (abs(secant_update) > 1e-10 &&
                   abs(delta_volume_nm1) > 1e-10 * initial_inside_cell_volume &&
                   secant_it < 20)
              {
                // If the cell is almost full or empty, we stop correcting the
                // volume.
                if (inside_cell_volume_nm1 < 1e-10 * cell_size ||
                    inside_cell_volume_nm1 > (cell_size - 1e-10 * cell_size))
                  {
                    eta_n = 0.0;
                    break;
                  }
                secant_it += 1;

                inside_cell_volume_n =
                  InterfaceTools::compute_cell_wise_volume(fe_point_evaluation,
                                                           cell,
                                                           cell_dof_values,
                                                           eta_n,
                                                           fe.degree + 1);

                delta_volume_n = targeted_cell_volume - inside_cell_volume_n;

                delta_volume_prime = (delta_volume_n - delta_volume_nm1) /
                                     (eta_n - eta_nm1 + 1e-16);

                secant_update = -delta_volume_n / (delta_volume_prime + 1e-16);
                eta_np1       = eta_n + secant_update;

                eta_nm1                = eta_n;
                eta_n                  = eta_np1;
                inside_cell_volume_nm1 = inside_cell_volume_n;
                delta_volume_nm1       = delta_volume_n;
              } // End secant method loop.

            if (secant_it >= 20)
              {
                eta_n = 0.0;
              }

            // Store eta_n
            eta_cell[cell_index] = eta_n;

            // Approximate L2 projection of the cell-wise (discontinuous)
            // correction to have a continuous correction at the dofs.
            // std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
            // cell->get_dof_indices(dof_indices);
            // for (unsigned int i = 0; i < dofs_per_cell; ++i)
            //   {
            //     volume_correction(dof_indices[i]) +=
            //       eta_n * n_cells_per_dofs_inv;
            //   }
          }
      } // End loop on cells.

    // Compute L2 projection of the cell-wise (discontinuous) correction to have
    // a continuous correction at the dofs.
    compute_volume_correction_L2_projection(eta_cell, mpi_communicator);

    volume_correction.compress(VectorOperation::add);
    volume_correction.update_ghost_values();
  }

  /**
   * @brief Compute the L2 projection of the cell-wise correction.
   *
   * @param[in] eta_cell Cell-wise correction of the signed distance
   *
   * @param[in] mpi_communicator MPI communicator
   */
  template <int dim>
  void
  SignedDistanceSolver<dim>::compute_volume_correction_L2_projection(
    Vector<double> &eta_cell,
    const MPI_Comm &mpi_communicator)
  {
    // Assemble System
    system_matrix_volume_correction = 0;
    system_rhs_volume_correction    = 0;

    FEValues<dim> fe_values(this->fe,
                            QGauss<dim>(fe.degree + 2),
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = fe_values.get_quadrature().size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            const unsigned int cell_index = cell->active_cell_index();

            // The cell is not intersected, no need to correct the volume
            if (interface_reconstruction_vertices.find(cell_index) ==
                interface_reconstruction_vertices.end())
              {
                continue;
              }

            fe_values.reinit(cell);

            cell_matrix = 0;
            cell_rhs    = 0;

            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
              {
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      cell_matrix(i, j) += (fe_values.shape_value(i, q_point) *
                                            fe_values.shape_value(j, q_point) *
                                            fe_values.JxW(q_point));
                    }

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    cell_rhs(i) +=
                      (fe_values.shape_value(i, q_point) *
                       eta_cell[cell_index] * fe_values.JxW(q_point));
                  }
              }
            cell->get_dof_indices(local_dof_indices);

            constraints.distribute_local_to_global(
              cell_matrix,
              cell_rhs,
              local_dof_indices,
              system_matrix_volume_correction,
              system_rhs_volume_correction);
          }
      }

    // Solve system
    const double       linear_solver_tolerance = 1e-12;
    const unsigned int max_iterations          = 100;

    VectorType completely_distributed_volume_correction(
      this->locally_owned_dofs, mpi_communicator);

    SolverControl solver_control(max_iterations,
                                 linear_solver_tolerance,
                                 true,
                                 true);

    TrilinosWrappers::SolverCG solver(solver_control);

    const double                                      ilu_fill = 0.0;
    const double                                      ilu_atol = 1e-12;
    const double                                      ilu_rtol = 1.0;
    TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
      ilu_fill, ilu_atol, ilu_rtol, 0);
    std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner =
      std::make_shared<TrilinosWrappers::PreconditionILU>();
    ilu_preconditioner->initialize(system_matrix_volume_correction,
                                   preconditionerOptions);

    solver.solve(system_matrix_volume_correction,
                 completely_distributed_volume_correction,
                 system_rhs_volume_correction,
                 *ilu_preconditioner);

    // Convert completely_distributed_volume_correction to dealii vector
    this->constraints.distribute(completely_distributed_volume_correction);
    LinearAlgebra::distributed::Vector<double> tmp_volume_correction(
      this->locally_owned_dofs, mpi_communicator);
    dealii::LinearAlgebra::ReadWriteVector<double> rwv(
      tmp_volume_correction.locally_owned_elements());
    rwv.reinit(completely_distributed_volume_correction);
    tmp_volume_correction.import_elements(rwv, dealii::VectorOperation::insert);

    volume_correction = tmp_volume_correction;
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::conserve_global_volume(
    const MPI_Comm &mpi_communicator)
  {
    /* We want to find a global correction function to apply to the dof value of
    the signed_distance so that the geometric global volume encompassed by the
    level 0 of the signed_distance V and by the iso-contour 0.5 of the phase
    fraction V_VOF match. This is required because the computed distance doesn't
    belong to the Q1 approximation space. We solve the non-linear problem:
    DeltaV(phi* + xi) = V_VOF - V(phi* + xi) = 0, where phi* is the
    redistanciated signed distance, xi = C*eta is the correction function on the
    signed_distance that we are looking for, with eta being the cell wise
    correction compute with compute_cell_wise_volume_correction() and C being a
    constant. We use the secant method to do so. See Ausas et al.
    (2010) for more details.*/

    /* Compute targeted global volume. It corresponds to the one enclosed by
    the level 0 of the level_set vector (same volume as the one enclosed
    by iso-contour 0.5 of the phase fraction).*/
    const double global_volume =
      compute_volume(mapping, dof_handler, fe, level_set, mpi_communicator);

    /* Initialization of values for the secant method. The subscript nm1 (or n
    minus 1) stands for the previous secant iteration (it = n-1), the
    subscript n stands for the current iteration and the subscript np1 stands
    for the next iteration (it = n+1).*/
    double global_volume_nm1 = 0.0;
    double global_volume_n   = 0.0;

    double global_delta_volume_nm1 = 0.0;
    double global_delta_volume_n   = 0.0;

    double global_delta_volume_prime = 0.0;

    /* Global constant C that we are solving for to obtain
        DeltaV(phi* + C*eta) = V_VOF - V(phi* + C*eta) = 0
    where eta is the cell wise correction computed with
    compute_cell_wise_volume_correction()*/
    double C_nm1 = 0.0;
    double C_n   = 0.0;
    double C_np1 = 0.0;

    /* Compute the volume and the difference with the targeted volume for 1st
    initial guess of the correction funtion (xi_nm1 = C_nm1*eta)*/
    C_nm1 = 1.0;
    LinearAlgebra::distributed::Vector<double> signed_distance_0(
      signed_distance_with_ghost);
    signed_distance_0.add(C_nm1, volume_correction);

    // Update_ghost_values is required for cell-wise volume computations
    signed_distance_0.update_ghost_values();

    global_volume_nm1 = compute_volume(
      mapping, dof_handler, fe, signed_distance_0, mpi_communicator);

    global_delta_volume_nm1 = global_volume - global_volume_nm1;

    // Initialize the 2nd initial guest
    C_n = 1e-6 * global_volume;

    // Store the initial volume for the stop criterion
    const double global_volume_0 = global_volume_nm1;

    // Initialize secant method it and update
    unsigned int secant_it     = 0;
    double       secant_update = 1.0;

    // Secant method
    while (abs(secant_update) > 1e-10 &&
           abs(global_delta_volume_nm1) > 1e-10 * global_volume_0 &&
           secant_it < 20)
      {
        secant_it += 1;

        LinearAlgebra::distributed::Vector<double> signed_distance_n(
          signed_distance_with_ghost);
        signed_distance_n.add(C_n, volume_correction);
        signed_distance_n.update_ghost_values();

        global_volume_n = compute_volume(
          mapping, dof_handler, fe, signed_distance_n, mpi_communicator);

        global_delta_volume_n = global_volume - global_volume_n;

        global_delta_volume_prime =
          (global_delta_volume_n - global_delta_volume_nm1) /
          (C_n - C_nm1 + 1e-16);

        secant_update =
          -global_delta_volume_n / (global_delta_volume_prime + 1e-16);

        C_np1 = C_n + secant_update;
        C_nm1 = C_n;
        C_n   = C_np1;

        global_volume_nm1       = global_volume_n;
        global_delta_volume_nm1 = global_delta_volume_n;
      }

    // If the secant method does not converge, do not correct.
    if (secant_it >= 20)
      C_n = 0.0;

    // Update signed_distance with the correction
    signed_distance.add(C_n, volume_correction);
    signed_distance.update_ghost_values();

    for (auto p : this->locally_active_dofs)
      {
        distance(p) = abs(signed_distance(p));
      }

    exchange_distance();
  }

  /**
   * @brief Set the map of local ids of the opposite faces to the given local dofs
   * (works for quad only).
   */
  template <int dim>
  inline void
  SignedDistanceSolver<dim>::set_dof_opposite_faces_map()
  {
    if constexpr (dim == 2)
      {
        dof_opposite_faces_map[0] = {1, 3};
        dof_opposite_faces_map[1] = {0, 3};
        dof_opposite_faces_map[2] = {1, 2};
        dof_opposite_faces_map[3] = {0, 2};

        if (fe.degree == 2)
          {
            dof_opposite_faces_map[4] = {1, 2, 3};
            dof_opposite_faces_map[5] = {0, 2, 3};
            dof_opposite_faces_map[6] = {0, 1, 3};
            dof_opposite_faces_map[7] = {0, 1, 2};
            dof_opposite_faces_map[8] = {0, 1, 2, 3};
          }
      }

    if constexpr (dim == 3)
      {
        dof_opposite_faces_map[0] = {1, 3, 5};
        dof_opposite_faces_map[1] = {0, 3, 5};
        dof_opposite_faces_map[2] = {1, 2, 5};
        dof_opposite_faces_map[3] = {0, 2, 5};
        dof_opposite_faces_map[4] = {1, 3, 4};
        dof_opposite_faces_map[5] = {0, 3, 4};
        dof_opposite_faces_map[6] = {1, 2, 4};
        dof_opposite_faces_map[7] = {0, 2, 4};

        if (fe.degree == 2)
          {
            dof_opposite_faces_map[8]  = {1, 2, 3, 5};
            dof_opposite_faces_map[9]  = {0, 2, 3, 5};
            dof_opposite_faces_map[10] = {0, 1, 3, 5};
            dof_opposite_faces_map[11] = {0, 1, 2, 5};
            dof_opposite_faces_map[12] = {1, 2, 3, 4};
            dof_opposite_faces_map[13] = {0, 2, 3, 4};
            dof_opposite_faces_map[14] = {0, 1, 3, 4};
            dof_opposite_faces_map[15] = {0, 1, 2, 4};
            dof_opposite_faces_map[16] = {1, 3, 4, 5};
            dof_opposite_faces_map[17] = {0, 3, 4, 5};
            dof_opposite_faces_map[18] = {1, 2, 4, 5};
            dof_opposite_faces_map[19] = {0, 2, 4, 5};
            dof_opposite_faces_map[20] = {1, 2, 3, 4, 5};
            dof_opposite_faces_map[21] = {0, 2, 3, 4, 5};
            dof_opposite_faces_map[22] = {0, 1, 3, 4, 5};
            dof_opposite_faces_map[23] = {0, 1, 2, 4, 5};
            dof_opposite_faces_map[24] = {0, 1, 2, 3, 5};
            dof_opposite_faces_map[25] = {0, 1, 2, 3, 4};
            dof_opposite_faces_map[26] = {0, 1, 2, 3, 4, 5};
          }
      }
  }

  /**
   * @brief Return the local ids of the opposite faces to the given dof
   * (works for quad only).
   *
   * @param[in] local_dof_id Local id of the dof in the cell
   *
   * @param[out] local_opposite_faces The vector containing the id of the
   * opposite faces
   */
  template <int dim>
  inline void
  SignedDistanceSolver<dim>::get_dof_opposite_faces(
    unsigned int               local_dof_id,
    std::vector<unsigned int> &local_opposite_faces)
  {
    local_opposite_faces = dof_opposite_faces_map.at(local_dof_id);
  }

  /**
   * @brief Return the number of opposite faces for the given dof
   * (works for quad only).
   *
   * @param[in] local_dof_id Local id of the dof in the cell
   *
   * @return n_opposite_faces_per_dof Number of opposite faces
   */
  template <int dim>
  inline unsigned int
  SignedDistanceSolver<dim>::get_n_opposite_faces_per_dof(
    unsigned int local_dof_id)
  {
    if (local_dof_id < 4)
      return dim;
    else if (local_dof_id < 8)
      return 3;
    else if (local_dof_id < 20)
      return 4;
    else if (local_dof_id < 26)
      return 5;
    else
      return 6;
  }

  /**
   * @brief
   * Return the face transformation jacobian (dim-1 x dim-1).
   *
   * @param[in] cell_transformation_jac Transformation jacobian of the cell (dim
   * x dim)
   *
   * @param[in] local_face_id Local id of the face
   *
   * @param[out] face_transformation_jac Face transformation jacobian (dim-1 x
   * dim-1) faces
   */
  template <int dim>
  inline void
  SignedDistanceSolver<dim>::get_face_transformation_jacobian(
    const DerivativeForm<1, dim, dim> &cell_transformation_jac,
    const unsigned int                 local_face_id,
    DerivativeForm<1, dim - 1, dim>   &face_transformation_jac)
  {
    for (unsigned int i = 0; i < dim; ++i)
      {
        unsigned int k = 0;
        for (unsigned int j = 0; j < dim; ++j)
          {
            if (local_face_id / 2 == j)
              continue;
            face_transformation_jac[i][k] = cell_transformation_jac[i][j];
            k += 1;
          }
      }
  }

  /**
   * @brief
   * Transform a point dim-1 in a reference face to a point dim in the
   * reference cell.
   *
   * @param[in] x_ref_face Point dim-1 in the reference face
   *
   * @param[in] local_face_id Local id of the face
   *
   * @return Point dim in the reference cell
   */
  template <int dim>
  inline Point<dim>
  SignedDistanceSolver<dim>::transform_ref_face_point_to_ref_cell(
    const Point<dim - 1> &x_ref_face,
    const unsigned int    local_face_id)
  {
    Point<dim> x_ref_cell;

    unsigned int j = 0;
    for (unsigned int i = 0; i < dim; ++i)
      {
        if (local_face_id / 2 == i)
          {
            x_ref_cell[i] = double(local_face_id % 2);
            continue;
          }
        x_ref_cell[i] = x_ref_face[j];
        j += 1;
      }

    return x_ref_cell;
  }

  template <int dim>
  inline void
  SignedDistanceSolver<dim>::compute_residual(
    const Tensor<1, dim>                  &x_n_to_x_I_real,
    const Tensor<1, dim>                  &distance_gradient,
    const DerivativeForm<1, dim - 1, dim> &transformation_jac,
    Tensor<1, dim - 1>                    &residual_ref)
  {
    Tensor<1, dim> residual_real =
      distance_gradient - (1.0 / x_n_to_x_I_real.norm()) * x_n_to_x_I_real;

    DerivativeForm<1, dim, dim - 1> transformation_jac_transpose =
      transformation_jac.transpose();

    for (unsigned int i = 0; i < dim - 1; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        residual_ref[i] +=
          transformation_jac_transpose[i][j] * residual_real[j];
  }

  template <int dim>
  inline std::vector<Point<dim>>
  SignedDistanceSolver<dim>::compute_numerical_jacobian_stencil(
    const Point<dim>  &x_ref,
    const unsigned int local_face_id,
    const double       perturbation)
  {
    std::vector<Point<dim>> stencil(2 * dim - 1);
    for (unsigned int i = 0; i < 2 * dim - 1; ++i)
      {
        stencil[i] = x_ref;
      }

    unsigned int              skip_index = local_face_id / 2;
    std::vector<unsigned int> j_index(dim - 1);

    // Set the coordinates (x,y,z) to be skipped (the coordinates not perturbed)
    unsigned int j = 0;
    for (unsigned int i = 0; i < dim; ++i)
      {
        if (i == skip_index)
          continue;
        j_index[j] = i;
        j += 1;
      }

    for (unsigned int i = 1; i < dim; ++i)
      {
        j = j_index[i - 1];
        stencil[2 * i - 1][j] -= perturbation;
        stencil[2 * i][j] += perturbation;
      }

    return stencil;
  }

  template <int dim>
  inline Tensor<1, dim>
  SignedDistanceSolver<dim>::transform_ref_face_correction_to_ref_cell(
    const Vector<double> &correction_ref_face,
    const unsigned int    local_face_id)
  {
    Tensor<1, dim> correction_ref_cell;

    unsigned int j = 0;
    for (unsigned int i = 0; i < dim; ++i)
      {
        if (local_face_id / 2 == i)
          {
            correction_ref_cell[i] = 0.0;
            continue;
          }
        correction_ref_cell[i] = correction_ref_face[j];
        j += 1;
      }

    return correction_ref_cell;
  }

  template <int dim>
  inline void
  SignedDistanceSolver<dim>::compute_numerical_jacobian(
    const std::vector<Point<dim>>     &stencil_real,
    const Point<dim>                  &x_I_real,
    const std::vector<Tensor<1, dim>> &distance_gradients,
    const std::vector<DerivativeForm<1, dim - 1, dim>>
                             &transformation_jacobians,
    const double              perturbation,
    LAPACKFullMatrix<double> &jacobian_matrix)
  {
    for (unsigned int i = 0; i < dim - 1; ++i)
      {
        const Tensor<1, dim> x_n_to_x_I_real_m1 =
          x_I_real - stencil_real[2 * i + 1];

        Tensor<1, dim - 1> residual_ref_m1;
        compute_residual(x_n_to_x_I_real_m1,
                         distance_gradients[2 * i + 1],
                         transformation_jacobians[2 * i + 1],
                         residual_ref_m1);

        const Tensor<1, dim> x_n_to_x_I_real_p1 =
          x_I_real - stencil_real[2 * i + 2];

        Tensor<1, dim - 1> residual_ref_p1;
        compute_residual(x_n_to_x_I_real_p1,
                         distance_gradients[2 * i + 2],
                         transformation_jacobians[2 * i + 2],
                         residual_ref_p1);

        for (unsigned int j = 0; j < dim - 1; ++j)
          {
            jacobian_matrix.set(j,
                                i,
                                (residual_ref_p1[j] - residual_ref_m1[j]) /
                                  (2.0 * perturbation));
          }
      }
  }

  template <int dim>
  inline double
  SignedDistanceSolver<dim>::compute_distance(
    const Tensor<1, dim> &x_n_to_x_I_real,
    const double          distance)
  {
    return distance + x_n_to_x_I_real.norm();
  }

  template <int dim>
  void
  SignedDistanceSolver<dim>::output_interface_reconstruction(
    const std::string  output_name,
    const std::string  output_path,
    const unsigned int it) const
  {
    Visualization<dim> intersection_data_out;

    intersection_data_out.build_patches(interface_reconstruction_vertices);

    intersection_data_out.write_vtu_with_pvtu_record(
      output_path, output_name, it, MPI_COMM_WORLD, 3);
  }
} // namespace InterfaceTools

template <int dim>
class AdvectionField : public TensorFunction<1, dim>
{
public:
  virtual Tensor<1, dim>
  value(const Point<dim> &p) const override;
};

template <int dim>
Tensor<1, dim>
AdvectionField<dim>::value(const Point<dim> &p) const
{
  const double   period = 2.0;
  Tensor<1, dim> value;
  value[0] =
    -(Utilities::pow(sin(numbers::PI * p[0]), 2) * sin(2 * numbers::PI * p[1]) *
      cos(numbers::PI * this->get_time() / period));
  value[1] = Utilities::pow(sin(numbers::PI * p[1]), 2) *
             sin(2 * numbers::PI * p[0]) *
             cos(numbers::PI * this->get_time() / period);

  return value;
}

template <int dim>
class InitialConditions : public Function<dim>
{
public:
  InitialConditions(const Point<dim> p_center_point,
                    const double     p_tanh_thickness);
  double
  value(const Point<dim> &p, const unsigned int component = 0) const override;

private:
  const Point<dim> center_point;
  const double     tanh_thickness;
};

template <int dim>
InitialConditions<dim>::InitialConditions(const Point<dim> p_center_point,
                                          const double     p_tanh_thickness)
  : center_point(p_center_point)
  , tanh_thickness(p_tanh_thickness)
{}

template <int dim>
double
InitialConditions<dim>::value(const Point<dim>  &p,
                              const unsigned int component) const
{
  (void)component;
  Assert(component == 0, ExcIndexRange(component, 0, 1));

  Tensor<1, dim> dist = center_point - p;

  return 0.5 - 0.5 * std::tanh((dist.norm() - 0.15) / tanh_thickness);
}

template <int dim>
class AdvectionProblem
{
public:
  AdvectionProblem(const Settings &parameters);
  void
  run();

private:
  using VectorType = TrilinosWrappers::MPI::Vector;
  using MatrixType = TrilinosWrappers::SparseMatrix;

  void
  make_grid();
  void
  setup_system();

  struct AssemblyScratchData
  {
    // Constructors
    AssemblyScratchData(const FiniteElement<dim> &fe);
    AssemblyScratchData(const AssemblyScratchData &scratch_data);

    // Set up FEValues to reuse them because their initialization is expensive.
    FEValues<dim> fe_values;

    // Previous phase values
    std::vector<double> previous_phase_values;

    std::vector<Tensor<1, dim>> advection_directions;

    // Velocity
    AdvectionField<dim> advection_field;
  };

  struct AssemblyCopyData
  {
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
  };

  void
  assemble_system();
  void
  local_assemble_system(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    AssemblyScratchData                                  &scratch,
    AssemblyCopyData                                     &copy_data);
  void
  copy_local_to_global(const AssemblyCopyData &copy_data);

  void
  set_initial_conditions();
  void
  solve();

  void
  compute_level_set_from_phase_fraction();
  void
  compute_phase_fraction_from_level_set();

  void
  reinitialize_phase_fraction_with_geometric_method();

  double
  monitor_volume(unsigned int time_iteration);

  void
  refine_grid();
  void
  output_results(const int time_iteration) const;


  parallel::distributed::Triangulation<dim> triangulation;
  const MappingQ<dim>                       mapping;

  const FE_Q<dim> fe;
  DoFHandler<dim> dof_handler;

  AffineConstraints<double> constraints;
  MatrixType                system_matrix;

  IndexSet locally_owned_dofs;

  IndexSet locally_relevant_dofs;
  IndexSet locally_active_dofs;

  VectorType solution;
  VectorType previous_solution;
  VectorType system_rhs;

  double dt;
  double time = 0.0;

  MPI_Comm           mpi_communicator;
  ConditionalOStream pcout;

  VectorType level_set;

  InterfaceTools::SignedDistanceSolver<dim> signed_distance_solver;

  TableHandler table_volume_monitoring;
  TableHandler table_error_monitoring;
  double       initial_volume;

  double   tanh_thickness;
  Settings parameters;
};

enum ActiveFEIndex
{
  lagrange = 0,
  nothing  = 1
};

template <int dim>
AdvectionProblem<dim>::AdvectionProblem(const Settings &p_parameters)
  : triangulation(MPI_COMM_WORLD,
                  typename Triangulation<dim>::MeshSmoothing(
                    Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening))
  , mapping(1)
  , fe(1)
  , dof_handler(triangulation)
  , dt(p_parameters.time_step)
  , mpi_communicator(MPI_COMM_WORLD)
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
  , signed_distance_solver(triangulation,
                           fe,
                           p_parameters.max_reinitialization_distance)
  , tanh_thickness(p_parameters.tanh_thickness)
  , parameters(p_parameters)
{}

template <int dim>
void
AdvectionProblem<dim>::make_grid()
{
  Point<dim> p_0 = Point<dim>();
  p_0[0]         = 0;
  for (unsigned int i = 1; i < dim; ++i)
    p_0[i] = 0;

  Point<dim> p_1 = Point<dim>();
  p_1[0]         = 1;
  for (unsigned int i = 1; i < dim; ++i)
    p_1[i] = 1;

  std::vector<unsigned int> repetitions(dim);
  repetitions[0] = 1;
  for (unsigned int i = 1; i < dim; ++i)
    repetitions[i] = 1;

  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            repetitions,
                                            p_0,
                                            p_1);

  triangulation.refine_global(parameters.initial_refinement);
}

template <int dim>
void
AdvectionProblem<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
  locally_active_dofs   = DoFTools::extract_locally_active_dofs(dof_handler);

  solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  previous_solution.reinit(locally_owned_dofs,
                           locally_relevant_dofs,
                           mpi_communicator);
  level_set.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs =*/false);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             dof_handler.locally_owned_dofs(),
                                             mpi_communicator,
                                             locally_relevant_dofs);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);
}

template <int dim>
void
AdvectionProblem<dim>::assemble_system()
{
  system_matrix = 0;
  system_rhs    = 0;

  WorkStream::run(dof_handler.begin_active(),
                  dof_handler.end(),
                  *this,
                  &AdvectionProblem::local_assemble_system,
                  &AdvectionProblem::copy_local_to_global,
                  AssemblyScratchData(fe),
                  AssemblyCopyData());

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

template <int dim>
AdvectionProblem<dim>::AssemblyScratchData::AssemblyScratchData(
  const FiniteElement<dim> &fe)
  : fe_values(fe,
              QGauss<dim>(fe.degree + 2),
              update_values | update_gradients | update_quadrature_points |
                update_JxW_values)
  , advection_directions(fe_values.get_quadrature().size())
{
  const unsigned int n_q_points = fe_values.get_quadrature().size();

  this->previous_phase_values = std::vector<double>(n_q_points);
}

template <int dim>
AdvectionProblem<dim>::AssemblyScratchData::AssemblyScratchData(
  const AssemblyScratchData &scratch_data)
  : fe_values(scratch_data.fe_values.get_fe(),
              scratch_data.fe_values.get_quadrature(),
              update_values | update_gradients | update_quadrature_points |
                update_JxW_values)
  , advection_directions(fe_values.get_quadrature().size())
{
  const unsigned int n_q_points = fe_values.get_quadrature().size();

  this->previous_phase_values = std::vector<double>(n_q_points);
}

template <int dim>
void
AdvectionProblem<dim>::local_assemble_system(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  AssemblyScratchData                                  &scratch_data,
  AssemblyCopyData                                     &copy_data)
{
  if (!cell->is_locally_owned())
    return;

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q_points =
    scratch_data.fe_values.get_quadrature().size();

  copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
  copy_data.cell_rhs.reinit(dofs_per_cell);

  copy_data.local_dof_indices.resize(dofs_per_cell);

  scratch_data.fe_values.reinit(cell);
  scratch_data.fe_values.get_function_values(
    previous_solution, scratch_data.previous_phase_values);

  scratch_data.advection_field.set_time(time);
  scratch_data.advection_field.value_list(
    scratch_data.fe_values.get_quadrature_points(),
    scratch_data.advection_directions);

  double cell_size;

  if (dim == 2)
    {
      cell_size = std::sqrt(4. * cell->measure() / M_PI);
    }
  else if (dim == 3)
    {
      cell_size = std::pow(6 * cell->measure() / M_PI, 1. / 3.);
    }

  double dt_inv = 1.0 / this->dt;

  const auto &sd = scratch_data;
  for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
    {
      const Tensor<1, dim> velocity = sd.advection_directions[q_point];

      const double u_mag = std::max(velocity.norm(), 1e-12);

      const double tau =
        1. / std::sqrt(
               Utilities::fixed_power<2>(dt_inv) +
               Utilities::fixed_power<2>(2. * u_mag / (cell_size / fe.degree)));
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              // LHS

              // Time integration
              copy_data.cell_matrix(i, j) +=
                dt_inv * sd.fe_values.shape_value(i, q_point) *
                sd.fe_values.shape_value(j, q_point) *
                sd.fe_values.JxW(q_point);

              // Advective term
              copy_data.cell_matrix(i, j) +=
                sd.fe_values.shape_value(i, q_point) *
                sd.advection_directions[q_point] *
                sd.fe_values.shape_grad(j, q_point) * sd.fe_values.JxW(q_point);

              // Stabilization term
              copy_data.cell_matrix(i, j) +=
                tau * sd.advection_directions[q_point] *
                sd.fe_values.shape_grad(i, q_point) *
                (sd.advection_directions[q_point] *
                   sd.fe_values.shape_grad(j, q_point) +
                 sd.fe_values.shape_value(j, q_point) * dt_inv) *
                sd.fe_values.JxW(q_point);
            }
          // RHS

          // Stabilization term
          copy_data.cell_rhs(i) += tau * sd.advection_directions[q_point] *
                                   sd.fe_values.shape_grad(i, q_point) *
                                   dt_inv *
                                   scratch_data.previous_phase_values[q_point] *
                                   sd.fe_values.JxW(q_point);

          copy_data.cell_rhs(i) += dt_inv *
                                   sd.fe_values.shape_value(i, q_point) *
                                   scratch_data.previous_phase_values[q_point] *
                                   sd.fe_values.JxW(q_point);
        }
    }
  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
AdvectionProblem<dim>::copy_local_to_global(const AssemblyCopyData &copy_data)
{
  constraints.distribute_local_to_global(copy_data.cell_matrix,
                                         copy_data.cell_rhs,
                                         copy_data.local_dof_indices,
                                         this->system_matrix,
                                         this->system_rhs);
}

template <int dim>
void
AdvectionProblem<dim>::set_initial_conditions()
{
  Point<dim> center = Point<dim>();

  center(0) = 0.5;
  center(1) = 0.75;
  if constexpr (dim == 3)
    center(2) = 0.5;

  VectorType completely_distributed_solution(this->locally_owned_dofs,
                                             mpi_communicator);
  VectorTools::interpolate(this->mapping,
                           this->dof_handler,
                           InitialConditions<dim>(center, tanh_thickness),
                           completely_distributed_solution);
  this->constraints.distribute(completely_distributed_solution);

  this->solution          = completely_distributed_solution;
  this->previous_solution = this->solution;
}

template <int dim>
void
AdvectionProblem<dim>::solve()
{
  SolverControl solver_control(1000, 1e-12 * system_rhs.l2_norm());
  TrilinosWrappers::SolverGMRES solver(solver_control);

  TrilinosWrappers::PreconditionILU                 preconditioner;
  TrilinosWrappers::PreconditionILU::AdditionalData data_ilu;

  preconditioner.initialize(system_matrix, data_ilu);

  VectorType completely_distributed_solution(this->locally_owned_dofs,
                                             mpi_communicator);

  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               preconditioner);

  VectorType residual(locally_owned_dofs, mpi_communicator);

  system_matrix.vmult(residual, completely_distributed_solution);
  pcout << "   Iterations required for convergence: "
        << solver_control.last_step() << '\n'
        << "   Max norm of residual:                " << residual.linfty_norm()
        << '\n';

  constraints.distribute(completely_distributed_solution);
  solution = completely_distributed_solution;
}

template <int dim>
void
AdvectionProblem<dim>::refine_grid()
{
  SolutionTransfer<dim, VectorType> solution_trans(dof_handler);

  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate(
    mapping,
    dof_handler,
    QGauss<dim - 1>(fe.degree + 1),
    typename std::map<types::boundary_id, const Function<dim, double> *>(),
    solution,
    estimated_error_per_cell,
    ComponentMask(),
    nullptr,
    0,
    triangulation.locally_owned_subdomain());

  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
    triangulation, estimated_error_per_cell, 0.9, 0.00005);

  if (triangulation.n_levels() > parameters.max_refinement_level)
    for (auto &cell : triangulation.active_cell_iterators_on_level(
           parameters.max_refinement_level))
      cell->clear_refine_flag();
  for (auto &cell : triangulation.active_cell_iterators_on_level(
         parameters.min_refinement_level))
    cell->clear_coarsen_flag();

  triangulation.prepare_coarsening_and_refinement();

  solution_trans.prepare_for_coarsening_and_refinement(solution);

  triangulation.execute_coarsening_and_refinement();

  setup_system();

  VectorType tmp_solution(this->locally_owned_dofs, mpi_communicator);

  solution_trans.interpolate(tmp_solution);

  constraints.distribute(tmp_solution);

  solution = tmp_solution;
}

template <int dim>
void
AdvectionProblem<dim>::compute_phase_fraction_from_level_set()
{
  VectorType solution_owned(this->locally_owned_dofs, mpi_communicator);
  std::map<types::global_dof_index, Point<dim>> dof_support_points =
    DoFTools::map_dofs_to_support_points(mapping, dof_handler);

  for (auto p : this->locally_owned_dofs)
    {
      const double signed_dist = level_set[p];
      solution_owned[p] = 0.5 - 0.5 * std::tanh(signed_dist / tanh_thickness);
    }
  constraints.distribute(solution_owned);

  solution = solution_owned;
}

template <int dim>
void
AdvectionProblem<dim>::compute_level_set_from_phase_fraction()
{
  VectorType level_set_owned(this->locally_owned_dofs, mpi_communicator);
  for (auto p : this->locally_owned_dofs)
    {
      const double phase      = solution[p];
      double       phase_sign = sgn(0.5 - phase);
      level_set_owned[p] =
        tanh_thickness *
        std::atanh(phase_sign * std::min(abs(0.5 - phase) / 0.5, 1.0 - 1e-12));
    }
  constraints.distribute(level_set_owned);

  level_set = level_set_owned;
}

template <int dim>
void
AdvectionProblem<dim>::reinitialize_phase_fraction_with_geometric_method()
{
  pcout << "In redistanciation..." << std::endl;

  compute_level_set_from_phase_fraction();

  signed_distance_solver.setup_dofs(mpi_communicator);

  signed_distance_solver.set_level_set_from_background_mesh(dof_handler,
                                                            level_set,
                                                            mpi_communicator);

  signed_distance_solver.solve(mpi_communicator);

  TrilinosWrappers::MPI::Vector tmp_local_level_set(this->locally_owned_dofs,
                                                    mpi_communicator);

  VectorTools::interpolate_to_different_mesh(
    signed_distance_solver.dof_handler,
    signed_distance_solver.get_level_set(mpi_communicator),
    dof_handler,
    tmp_local_level_set);

  level_set = tmp_local_level_set;

  compute_phase_fraction_from_level_set();
}

template <int dim>
double
AdvectionProblem<dim>::monitor_volume(unsigned int time_iteration)
{
  FEValues<dim> fe_values(fe,
                          QGauss<dim>(fe.degree + 1),
                          update_values | update_JxW_values);

  const unsigned int  n_q_points = fe_values.n_quadrature_points;
  std::vector<double> phase_values(n_q_points);

  double volume_sharp = 0.0;
  double volume_phase = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(solution, phase_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              volume_phase += fe_values.JxW(q) * phase_values[q];
            }
        }
    }
  volume_sharp = InterfaceTools::compute_volume(
    mapping, dof_handler, fe, level_set, mpi_communicator);
  volume_phase = Utilities::MPI::sum(volume_phase, mpi_communicator);

  table_volume_monitoring.add_value("time_iteration", time_iteration);

  table_volume_monitoring.add_value("time", time);
  table_volume_monitoring.set_scientific("time", true);

  table_volume_monitoring.add_value("volume_sharp", volume_sharp);
  table_volume_monitoring.set_scientific("volume_sharp", true);

  table_volume_monitoring.add_value("volume_phase", volume_phase);
  table_volume_monitoring.set_scientific("volume_phase", true);

  std::ofstream output(parameters.output_path + "/volume.dat");
  table_volume_monitoring.write_text(output);

  pcout << "Volume = " << volume_sharp << std::endl;

  return volume_sharp;
}

template <int dim>
void
AdvectionProblem<dim>::output_results(const int time_iteration) const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  data_out.add_data_vector(solution, "solution");
  data_out.add_data_vector(previous_solution, "previous_solution");
  data_out.add_data_vector(level_set, "level_set");

  data_out.build_patches(fe.degree);

  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::CompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);

  const std::string output_path = parameters.output_path;
  const std::string filename    = parameters.output_name;

  data_out.write_vtu_with_pvtu_record(
    output_path, filename, time_iteration, MPI_COMM_WORLD, 3);

  signed_distance_solver.output_interface_reconstruction(
    "interface_" + filename, output_path, time_iteration);
}

template <int dim>
void
AdvectionProblem<dim>::run()
{
  pcout << "Setup system..." << std::endl;
  make_grid();

  setup_system();
  set_initial_conditions();

  unsigned int it         = 0;
  double       final_time = parameters.time_end;

  for (unsigned int i = 0; i < parameters.initial_refinement_steps; ++i)
    refine_grid();

  previous_solution = solution;

  reinitialize_phase_fraction_with_geometric_method();

  output_results(it);
  monitor_volume(it);

  pcout << "Solve system..." << std::endl;
  while (time + dt < final_time)
    {
      it += 1;
      time += this->dt;

      assemble_system();

      pcout << "it = " << it << " time = " << time << std::endl;

      solve();

      reinitialize_phase_fraction_with_geometric_method();

      output_results(it);
      monitor_volume(it);

      refine_grid();

      previous_solution = solution;
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  Settings parameters;
  if (!parameters.parse((argc > 1) ? (argv[1]) : ""))
    return 0;

  try
    {
      switch (parameters.dimension)
        {
          case 2:
            {
              AdvectionProblem<2> advection_problem(parameters);
              advection_problem.run();
              break;
            }
          case 3:
            {
              AdvectionProblem<3> advection_problem(parameters);
              advection_problem.run();
              break;
            }

          default:
            Assert(false, ExcMessage("This program only works in 2d and 3d."));
        }
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
