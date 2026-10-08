// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_anderson_jackson_filter_h
#define lethe_anderson_jackson_filter_h

#include <core/boundary_conditions.h>
#include <core/immersed_solid_classifier.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>
#include <core/periodic_boundary.h>
#include <core/pvd_handler.h>
#include <core/vector.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/rtree.h>

#include <array>
#include <map>
#include <memory>
#include <utility>
#include <vector>

/**
 * @brief Anderson-Jackson filter of a resolved fluid-solid snapshot.
 *
 * The filter computes, at the support points (the filter centers) of a
 * continuous FE_Q space defined on the triangulation of the fluid, the
 * phase-averaged fields of Anderson and Jackson (1967):
 * \f[
 *   M(\mathbf{x}) = \int_\Omega g(|\mathbf{x}-\mathbf{y}|) \mathrm{d}V, \quad
 *   E(\mathbf{x}) = \int_\Omega I_f g \mathrm{d}V, \quad
 *   \bar{\mathbf{u}}_f(\mathbf{x}) = \frac{\int_\Omega I_f \mathbf{u} g
 *   \mathrm{d}V}{E(\mathbf{x})}, \quad
 *   \bar{p}_f(\mathbf{x}) = \frac{\int_\Omega I_f p g
 * \mathrm{d}V}{E(\mathbf{x})}, \f] where \f$g\f$ is a kernel of unit mass with
 * compact support and \f$I_f\f$ is the fluid indicator. The fluid volume
 * fraction is \f$\epsilon_f = E\f$ and the solid volume fraction is
 * \f$\epsilon_s = M - E\f$. Where the support of the kernel does not reach a
 * non-periodic boundary of the domain, the kernel mass is exactly one, and
 * both volume fractions are divided by its discrete value, which removes the
 * quadrature error of the kernel and ensures that
 * \f$\epsilon_f + \epsilon_s = 1\f$. The same division is applied where
 * the discrete mass exceeds one, since any excess is a quadrature error.
 * Otherwise, where the kernel is truncated by a boundary, the volume
 * fractions are optionally divided by the kernel mass \f$M\f$, which
 * renormalizes the kernel.
 *
 * Optionally, the part of the kernel outside of the domain, beyond the walls
 * where the velocity is imposed, is filled with fluid moving at the velocity
 * of the wall. The mass \f$M_e\f$ of this exterior part is \f$1 - M\f$,
 * attributed to the walls in proportion to the integral of the kernel over
 * the walls and over the other non-periodic boundaries, and its velocity is
 * the average of the wall velocity weighted by the kernel. The exterior fluid
 * is counted both in the fluid volume fraction,
 * \f$\epsilon_f = E + M_e\f$, and in the phase-averaged velocity, so that
 * their product remains the flux of fluid seen by the kernel. The integrals
 * over the boundaries are accumulated by the same sweep as the integrals over
 * the cells, from the boundary faces of the locally owned cells.
 *
 * The convolution is source-centric: every process integrates its locally
 * owned cells (the sources) once, and adds the contribution of each batch of
 * quadrature points to the filter centers within the support of the kernel.
 * The algorithm proceeds as follows:
 * 1. The filter centers are the locally owned and unconstrained degrees of
 *    freedom of the FE_Q space. The hanging nodes are recovered from the
 *    constraints afterwards.
 * 2. Every process describes its locally owned cells by a few bounding
 *    boxes, which are gathered on all the processes. Each filter center, and
 *    each of its periodic images, is sent to the processes whose boxes lie
 *    within the support radius of the kernel.
 * 3. Every process builds an R-tree of the filter centers it received and
 *    integrates its locally owned cells. The cells cut by an immersed solid
 *    are integrated with an iterated Gauss quadrature on which the fluid
 *    indicator is evaluated point by point.
 * 4. The partial integrals are returned to the owners of the filter centers,
 *    which sum and normalize them.
 *
 * The cost of the filter scales with the number of filter centers times the
 * number of quadrature points within the support of the kernel, i.e. as
 * \f$(R/h)^d\f$ per filter center, where \f$R\f$ is the support radius of the
 * kernel and \f$h\f$ is the size of the cells.
 *
 * The filter is independent of the solver that produced the solution: it only
 * requires the DoFHandler, the mapping and the solution of a velocity-pressure
 * finite element with dim+1 components, as well as a description of the
 * immersed solids through an ImmersedSolidClassifier.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class AndersonJacksonFilter
{
public:
  /// Type of the vectors that store the filtered fields.
  using FieldVectorType = dealii::LinearAlgebra::distributed::Vector<double>;

  /// Velocity imposed on the walls, keyed by boundary id. Each function has at
  /// least dim components, the first dim of which are the velocity.
  using WallVelocities = std::map<dealii::types::boundary_id,
                                  std::shared_ptr<const dealii::Function<dim>>>;

  /**
   * @brief Filtered fields, defined on the FE_Q space returned by
   * get_dof_handler(). The vectors contain the locally relevant values.
   */
  struct FilteredFields
  {
    /// Fluid volume fraction.
    FieldVectorType fluid_volume_fraction;

    /// Solid volume fraction.
    FieldVectorType solid_volume_fraction;

    /// Kernel mass, i.e. the integral of the kernel over the domain.
    FieldVectorType kernel_mass;

    /// Components of the phase-averaged fluid velocity.
    std::array<FieldVectorType, dim> velocity;

    /// Phase-averaged fluid pressure. It is zero if the pressure is not
    /// filtered.
    FieldVectorType pressure;

    /// Indicator of the filter centers where the phase-averaged velocity and
    /// pressure are defined (1) or not (0).
    FieldVectorType valid;
  };

  /**
   * @brief Global statistics of the last application of the filter. They are
   * identical on every process and independent of the partition of the mesh,
   * up to round-off errors.
   */
  struct FilterStatistics
  {
    /// Number of filter centers.
    dealii::types::global_dof_index n_filter_centers = 0;

    /// Number of source cells lying in the fluid.
    dealii::types::global_cell_index n_fluid_cells = 0;

    /// Number of source cells lying in an immersed solid.
    dealii::types::global_cell_index n_solid_cells = 0;

    /// Number of source cells cut by an immersed solid.
    dealii::types::global_cell_index n_cut_cells = 0;

    /// Number of filter centers where the phase-averaged velocity and
    /// pressure are not defined.
    dealii::types::global_dof_index n_undefined_averages = 0;

    /// Minimum of the kernel mass over the filter centers.
    double minimum_kernel_mass = 0.;

    /// Maximum of the kernel mass over the filter centers.
    double maximum_kernel_mass = 0.;

    /// Minimum of the fluid volume fraction over the filter centers.
    double minimum_fluid_volume_fraction = 0.;

    /// Maximum of the fluid volume fraction over the filter centers.
    double maximum_fluid_volume_fraction = 0.;

    /// Volume of the source cells integrated with the source quadratures.
    double source_volume = 0.;

    /// Volume of the immersed solids integrated with the source quadratures.
    double source_solid_volume = 0.;

    /// Integral of the solid volume fraction over the domain. For a periodic
    /// domain, it equals the volume of the solids up to discretization
    /// errors.
    double solid_volume_fraction_integral = 0.;
  };

  /**
   * @brief Constructor.
   *
   * @param[in] filter_parameters Parameters of the filter.
   *
   * @param[in] timer_parameters Parameters controlling the timer summary.
   *
   * @param[in] mpi_communicator MPI communicator of the triangulation.
   */
  AndersonJacksonFilter(
    const Parameters::AndersonJacksonFilter &filter_parameters,
    const Parameters::Timer                 &timer_parameters,
    const MPI_Comm                           mpi_communicator);

  /**
   * @brief Filter a velocity-pressure solution.
   *
   * @param[in] fluid_dof_handler DoFHandler of the velocity-pressure solution.
   * Its finite element has dim velocity components followed by the pressure.
   * The filtered fields are defined on the same triangulation.
   *
   * @param[in] mapping Mapping of the triangulation.
   *
   * @param[in] fluid_solution Velocity-pressure solution with its ghost
   * values.
   *
   * @param[in] classifier Description of the immersed solids.
   *
   * @param[in] periodic_boundaries Pairs of periodic boundaries of the
   * triangulation. The support radius of the kernel may exceed the period:
   * the kernel is then periodized, i.e. a source point contributes to a
   * filter center through several of its periodic images.
   *
   * @param[in] wall_velocities Velocity imposed on the walls, used to fill the
   * part of the kernel beyond the walls with fluid when the parameters request
   * it. The boundaries that are not listed truncate the kernel.
   */
  void
  apply(const dealii::DoFHandler<dim>        &fluid_dof_handler,
        const dealii::Mapping<dim>           &mapping,
        const GlobalVectorType               &fluid_solution,
        const ImmersedSolidClassifier<dim>   &classifier,
        const Parameters::PeriodicBoundaries &periodic_boundaries,
        const WallVelocities                 &wall_velocities = {});

  /**
   * @brief Write the filtered fields in the output folder, as vtu files and a
   * pvd record.
   *
   * @param[in] mapping Mapping of the triangulation.
   *
   * @param[in] time Time associated with the filtered fields.
   *
   * @param[in] group_files Number of vtu files written.
   */
  void
  write_output(const dealii::Mapping<dim> &mapping,
               const double                time,
               const unsigned int          group_files);

  /**
   * @brief Print the summary of the timer, unless the timer is disabled.
   */
  void
  print_timer_summary() const;

  /**
   * @brief Return the timer of the filter, in which the caller can time its
   * own sections, e.g. the loading of the snapshot.
   *
   * @return Timer of the filter.
   */
  dealii::TimerOutput &
  get_computing_timer()
  {
    return computing_timer;
  }

  /**
   * @brief Return the DoFHandler of the FE_Q space of the filtered fields.
   *
   * @return DoFHandler of the filtered fields.
   */
  const dealii::DoFHandler<dim> &
  get_dof_handler() const
  {
    return dof_handler;
  }

  /**
   * @brief Return the filtered fields of the last application of the filter.
   *
   * @return Filtered fields.
   */
  const FilteredFields &
  get_filtered_fields() const
  {
    return filtered_fields;
  }

  /**
   * @brief Return the statistics of the last application of the filter.
   *
   * @return Global statistics of the filter.
   */
  const FilterStatistics &
  get_statistics() const
  {
    return statistics;
  }

private:
  /**
   * @brief Moments of the kernel accumulated at a filter center.
   */
  struct FilterMoments
  {
    /// Integral of the kernel.
    double kernel_mass = 0.;

    /// Integral of the kernel times the fluid indicator.
    double fluid_volume = 0.;

    /// Integral of the kernel times the fluid indicator times the velocity.
    dealii::Tensor<1, dim> velocity_moment;

    /// Integral of the kernel times the fluid indicator times the pressure.
    double pressure_moment = 0.;

    /// Integral of the kernel over the non-periodic boundaries. It is zero
    /// where the support of the kernel does not reach these boundaries.
    double boundary_weight = 0.;

    /// Integral of the kernel over the walls where the velocity is imposed.
    double wall_weight = 0.;

    /// Integral of the kernel times the imposed velocity over the walls.
    dealii::Tensor<1, dim> wall_velocity_moment;

    /**
     * @brief Add the moments of another contribution.
     *
     * @param[in] other Moments to add.
     *
     * @return Reference to the updated moments.
     */
    FilterMoments &
    operator+=(const FilterMoments &other);
  };

  /**
   * @brief Batch of source quadrature points, stored as a structure of arrays
   * so that the loop over the points of the batch can be vectorized. The
   * points of consecutive cells are gathered in the same batch so that a
   * single search of the filter centers serves several cells.
   */
  struct SourceBatch
  {
    /// Coordinates of the points.
    std::array<std::vector<double>, dim> coordinates;

    /// Quadrature weights times the Jacobian determinant.
    std::vector<double> weights;

    /// Weights times the fluid indicator.
    std::vector<double> fluid_weights;

    /// Fluid weights times the components of the velocity.
    std::array<std::vector<double>, dim> velocity_moments;

    /// Fluid weights times the pressure.
    std::vector<double> pressure_moments;

    /// Lower corner of the bounding box of the points.
    dealii::Point<dim> lower_corner;

    /// Upper corner of the bounding box of the points.
    dealii::Point<dim> upper_corner;

    /**
     * @brief Return the number of points of the batch.
     *
     * @return Number of points.
     */
    unsigned int
    size() const
    {
      return weights.size();
    }

    /**
     * @brief Remove all the points of the batch, keeping the allocated memory.
     */
    void
    clear();
  };

  /**
   * @brief Contiguous range of filter centers received from the same
   * process.
   */
  struct OwnerBlock
  {
    /// Process that owns the filter centers.
    unsigned int process;

    /// Index of the first filter center of the range.
    unsigned int begin;

    /// Index past the last filter center of the range.
    unsigned int end;
  };

  /// R-tree of the filter centers integrated by this process, paired with
  /// their index.
  using TargetTree = dealii::RTree<std::pair<dealii::Point<dim>, unsigned int>>;

  /**
   * @brief Apply the filter with a given kernel.
   *
   * @tparam KernelType Type of the kernel.
   *
   * @param[in] kernel Kernel of the filter.
   *
   * @param[in] fluid_dof_handler DoFHandler of the velocity-pressure solution.
   *
   * @param[in] mapping Mapping of the triangulation.
   *
   * @param[in] fluid_solution Velocity-pressure solution with its ghost
   * values.
   *
   * @param[in] classifier Description of the immersed solids.
   *
   * @param[in] periodic_boundaries Pairs of periodic boundaries.
   *
   * @param[in] wall_velocities Velocity imposed on the walls.
   */
  template <typename KernelType>
  void
  apply_kernel(const KernelType                     &kernel,
               const dealii::DoFHandler<dim>        &fluid_dof_handler,
               const dealii::Mapping<dim>           &mapping,
               const GlobalVectorType               &fluid_solution,
               const ImmersedSolidClassifier<dim>   &classifier,
               const Parameters::PeriodicBoundaries &periodic_boundaries,
               const WallVelocities                 &wall_velocities);

  /**
   * @brief Set up the FE_Q space of the filtered fields, its hanging node
   * constraints and the vectors of the filtered fields.
   *
   * @param[in] fluid_dof_handler DoFHandler of the velocity-pressure solution.
   */
  void
  setup_filter_space(const dealii::DoFHandler<dim> &fluid_dof_handler);

  /**
   * @brief Generate the filter centers owned by this process, i.e. the
   * locally owned and unconstrained degrees of freedom of the FE_Q space.
   *
   * @param[in] mapping Mapping of the triangulation.
   */
  void
  generate_owned_targets(const dealii::Mapping<dim> &mapping);

  /**
   * @brief Send every owned filter center, and its periodic images, to the
   * processes whose locally owned cells lie within the support of the kernel.
   *
   * @param[in] mapping Mapping of the triangulation.
   *
   * @param[in] periodic_boundaries Pairs of periodic boundaries.
   */
  void
  distribute_targets(const dealii::Mapping<dim>           &mapping,
                     const Parameters::PeriodicBoundaries &periodic_boundaries);

  /**
   * @brief Build the R-tree of the filter centers integrated by this process.
   *
   * @return R-tree of the filter centers.
   */
  TargetTree
  build_target_tree();

  /**
   * @brief Integrate the locally owned cells and accumulate their
   * contribution to the filter centers within the support of the kernel.
   *
   * @tparam KernelType Type of the kernel.
   *
   * @param[in] kernel Kernel of the filter.
   *
   * @param[in] target_tree R-tree of the filter centers.
   *
   * @param[in] fluid_dof_handler DoFHandler of the velocity-pressure solution.
   *
   * @param[in] mapping Mapping of the triangulation.
   *
   * @param[in] fluid_solution Velocity-pressure solution with its ghost
   * values.
   *
   * @param[in] classifier Description of the immersed solids.
   *
   * @param[in] wall_velocities Velocity imposed on the walls. The kernel is
   * always integrated over the non-periodic boundary faces, to identify the
   * filter centers whose kernel is truncated by a boundary. The velocity of
   * the walls is only integrated if the parameters request the extension of
   * the velocity beyond the walls.
   */
  template <typename KernelType>
  void
  integrate_source_cells(const KernelType                   &kernel,
                         const TargetTree                   &target_tree,
                         const dealii::DoFHandler<dim>      &fluid_dof_handler,
                         const dealii::Mapping<dim>         &mapping,
                         const GlobalVectorType             &fluid_solution,
                         const ImmersedSolidClassifier<dim> &classifier,
                         const WallVelocities               &wall_velocities);

  /**
   * @brief Accumulate the integrals of the kernel over a boundary face at the
   * filter centers within the support of the kernel.
   *
   * @tparam KernelType Type of the kernel.
   *
   * @param[in] kernel Kernel of the filter.
   *
   * @param[in] target_tree R-tree of the filter centers.
   *
   * @param[in] points Quadrature points of the face.
   *
   * @param[in] weights Quadrature weights of the face times the surface
   * element.
   *
   * @param[in] wall_velocity_values Velocity imposed at the quadrature points
   * if the face belongs to a wall, empty otherwise.
   */
  template <typename KernelType>
  void
  accumulate_boundary_face(
    const KernelType                          &kernel,
    const TargetTree                          &target_tree,
    const std::vector<dealii::Point<dim>>     &points,
    const std::vector<double>                 &weights,
    const std::vector<dealii::Tensor<1, dim>> &wall_velocity_values);

  /**
   * @brief Accumulate the contribution of a batch of source points to the
   * filter centers within the support of the kernel, and clear the batch.
   *
   * @tparam KernelType Type of the kernel.
   *
   * @param[in] kernel Kernel of the filter.
   *
   * @param[in] target_tree R-tree of the filter centers.
   *
   * @param[in,out] batch Batch of source points, cleared on exit.
   */
  template <typename KernelType>
  void
  accumulate_batch(const KernelType &kernel,
                   const TargetTree &target_tree,
                   SourceBatch      &batch);

  /**
   * @brief Return the partial integrals to the owners of the filter centers,
   * which sum them.
   */
  void
  exchange_partial_integrals();

  /**
   * @brief Compute the filtered fields from the moments of the owned filter
   * centers, and fill the constrained and ghost values.
   */
  void
  normalize_filtered_fields();

  /**
   * @brief Reduce the statistics of the filter over all the processes.
   *
   * @param[in] mapping Mapping of the triangulation.
   */
  void
  compute_statistics(const dealii::Mapping<dim> &mapping);

  /**
   * @brief Print the statistics of the filter according to the verbosity.
   */
  void
  print_statistics() const;

  /// Parameters of the filter.
  const Parameters::AndersonJacksonFilter filter_parameters;

  /// Parameters of the timer.
  const Parameters::Timer timer_parameters;

  /// MPI communicator of the triangulation.
  const MPI_Comm mpi_communicator;

  /// Output stream of the first process.
  dealii::ConditionalOStream pcout;

  /// Timer of the stages of the filter.
  dealii::TimerOutput computing_timer;

  /// Finite element of the filtered fields.
  std::unique_ptr<dealii::FiniteElement<dim>> fe;

  /// DoFHandler of the filtered fields.
  dealii::DoFHandler<dim> dof_handler;

  /// Locally owned degrees of freedom of the filtered fields.
  dealii::IndexSet locally_owned_dofs;

  /// Locally relevant degrees of freedom of the filtered fields.
  dealii::IndexSet locally_relevant_dofs;

  /// Hanging node constraints of the filtered fields.
  dealii::AffineConstraints<double> constraints;

  /// Filtered fields.
  FilteredFields filtered_fields;

  /// Statistics of the last application of the filter.
  FilterStatistics statistics;

  /// Polynomial degree of the filtered fields.
  unsigned int output_degree = 0;

  /// Number of Gauss points per direction of the source quadratures.
  unsigned int n_quadrature_points = 0;

  /// Radius of the support of the kernel.
  double support_radius = 0.;

  /// Mass of the untruncated gaussian retained by the kernel (1 for the
  /// top-hat and Wendland kernels).
  double retained_mass_fraction = 1.;

  /// Locations of the filter centers owned by this process.
  std::vector<dealii::Point<dim>> owned_target_locations;

  /// Degrees of freedom of the filter centers owned by this process.
  std::vector<dealii::types::global_dof_index> owned_target_dofs;

  /// Moments of the filter centers owned by this process, summed over all
  /// the processes.
  std::vector<FilterMoments> owned_target_moments;

  /// Locations of the filter centers integrated by this process: the owned
  /// filter centers and their periodic images that lie within the support of
  /// the kernel of the owned cells, followed by the filter centers received
  /// from the other processes. It is released once the R-tree is built.
  std::vector<dealii::Point<dim>> integration_target_locations;

  /// Index of each integrated filter center in the list of owned filter
  /// centers of its owner.
  std::vector<unsigned int> integration_target_owner_indices;

  /// Partial moments of the integrated filter centers.
  std::vector<FilterMoments> integration_target_moments;

  /// Ranges of integrated filter centers sharing the same owner. The first
  /// range contains the filter centers owned by this process.
  std::vector<OwnerBlock> owner_blocks;

  /// Number of locally owned cells lying in the fluid.
  dealii::types::global_cell_index n_local_fluid_cells = 0;

  /// Number of locally owned cells lying in an immersed solid.
  dealii::types::global_cell_index n_local_solid_cells = 0;

  /// Number of locally owned cells cut by an immersed solid.
  dealii::types::global_cell_index n_local_cut_cells = 0;

  /// Volume of the locally owned cells integrated with the source
  /// quadratures.
  double local_source_volume = 0.;

  /// Volume of the immersed solids in the locally owned cells integrated with
  /// the source quadratures.
  double local_source_solid_volume = 0.;

  /// Minimum of the kernel mass over the owned filter centers.
  double local_minimum_kernel_mass = 0.;

  /// Maximum of the kernel mass over the owned filter centers.
  double local_maximum_kernel_mass = 0.;

  /// Minimum of the fluid volume fraction over the owned filter centers.
  double local_minimum_fluid_volume_fraction = 0.;

  /// Maximum of the fluid volume fraction over the owned filter centers.
  double local_maximum_fluid_volume_fraction = 0.;

  /// Number of owned filter centers where the phase averages are not
  /// defined.
  dealii::types::global_dof_index n_local_undefined_averages = 0;

  /// Record of the output files for the pvd file.
  PVDHandler pvd_handler;

  /// Number of outputs written.
  unsigned int n_outputs = 0;
};

/**
 * @brief Build the velocity imposed on the walls of a fluid dynamics
 * simulation, for the extension of the filtered velocity beyond the walls. The
 * walls are the boundaries with a noslip, function or function weak boundary
 * condition. The other boundaries (e.g. slip, outlets and periodic
 * boundaries) are not walls.
 *
 * @tparam dim Number of spatial dimensions.
 *
 * @param[in,out] boundary_conditions Boundary conditions of the fluid
 * dynamics. The time of their functions is set, and the returned functions
 * refer to them, so they must outlive the returned functions.
 *
 * @param[in] time Time at which the imposed velocity is evaluated.
 *
 * @return Velocity imposed on each wall, keyed by boundary id.
 */
template <int dim>
typename AndersonJacksonFilter<dim>::WallVelocities
make_wall_velocities(
  BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const double                                   time);

#endif
