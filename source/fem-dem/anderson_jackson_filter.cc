// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/filter_kernels.h>
#include <core/lethe_grid_tools.h>
#include <core/solutions_output.h>
#include <core/utilities.h>

#include <solvers/postprocessors.h>

#include <fem-dem/anderson_jackson_filter.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/data_out.h>

#include <boost/geometry/index/predicates.hpp>
#include <boost/iterator/function_output_iterator.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace dealii;

namespace
{
  /// Maximum number of points of a batch of source points.
  constexpr unsigned int maximum_batch_size = 1024;

  /// Maximum extent of a batch of source points, relative to the support
  /// radius of the kernel. The filter centers of a batch are searched within
  /// the support radius of the bounding box of the batch, so the number of
  /// kernel evaluations that return zero grows with the extent of the batch.
  constexpr double maximum_batch_extent_fraction = 0.1;

  /// Number of bins per direction used to group the locally owned cells of a
  /// process into a few bounding boxes.
  constexpr unsigned int n_box_bins_per_direction = 4;

  /// Relative enlargement of the bounding box of a cell given by the mapping.
  /// This box is built from the support points of the mapping, which can
  /// slightly underestimate the extent of a curved cell.
  constexpr double cell_box_relative_margin = 0.01;

  /// Number of digits after the decimal point of the printed statistics.
  constexpr unsigned int statistics_precision = 6;

  /// Width of the labels of the printed statistics, so that the values are
  /// aligned.
  constexpr std::size_t statistics_label_width = 42;

  /**
   * @brief Indent a label of the printed statistics and pad it with spaces so
   * that the values are aligned.
   *
   * @param[in] label Label of a statistic.
   *
   * @return Indented and padded label.
   */
  std::string
  statistics_label(const std::string &label)
  {
    const std::size_t n_spaces = (label.size() < statistics_label_width) ?
                                   statistics_label_width - label.size() :
                                   1;
    return "   " + label + std::string(n_spaces, ' ');
  }

  /**
   * @brief Request sent to a process to integrate a filter center.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  struct FilterTargetRequest
  {
    /// Location of the filter center, or of one of its periodic images.
    std::array<double, dim> location;

    /// Index of the filter center among the filter centers of its owner.
    unsigned int owner_index;

    /**
     * @brief Serialize the request. The requests are trivially copyable and
     * are sent bit for bit by Utilities::MPI::some_to_some(), but
     * Utilities::pack() also compiles its serialization path.
     *
     * @tparam Archive Type of the archive.
     *
     * @param[in,out] archive Archive.
     *
     * @param[in] version Version of the archive (unused).
     */
    template <class Archive>
    void
    serialize(Archive &archive, const unsigned int /*version*/)
    {
      for (double &coordinate : location)
        archive &coordinate;
      archive &owner_index;
    }
  };

  /**
   * @brief Partial moments of a filter center returned to its owner.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  struct FilterTargetContribution
  {
    /// Index of the filter center among the filter centers of its owner.
    unsigned int owner_index;

    /// Partial integral of the kernel.
    double kernel_mass;

    /// Partial integral of the kernel times the fluid indicator.
    double fluid_volume;

    /// Partial integral of the kernel times the fluid indicator times the
    /// velocity.
    std::array<double, dim> velocity_moment;

    /// Partial integral of the kernel times the fluid indicator times the
    /// pressure.
    double pressure_moment;

    /// Partial integral of the kernel over the non-periodic boundaries.
    double boundary_weight;

    /// Partial integral of the kernel over the walls.
    double wall_weight;

    /// Partial integral of the kernel times the velocity imposed on the walls.
    std::array<double, dim> wall_velocity_moment;

    /**
     * @brief Serialize the contribution. The contributions are trivially
     * copyable and are sent bit for bit by Utilities::MPI::some_to_some(), but
     * Utilities::pack() also compiles its serialization path.
     *
     * @tparam Archive Type of the archive.
     *
     * @param[in,out] archive Archive.
     *
     * @param[in] version Version of the archive (unused).
     */
    template <class Archive>
    void
    serialize(Archive &archive, const unsigned int /*version*/)
    {
      archive &owner_index;
      archive &kernel_mass;
      archive &fluid_volume;
      for (double &component : velocity_moment)
        archive &component;
      archive &pressure_moment;
      archive &boundary_weight;
      archive &wall_weight;
      for (double &component : wall_velocity_moment)
        archive &component;
    }
  };

  /**
   * @brief Pack the partial moments of a filter center to return them to its
   * owner.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @tparam MomentsType Type of the moments of the filter.
   *
   * @param[in] owner_index Index of the filter center among the filter
   * centers of its owner.
   *
   * @param[in] moments Partial moments of the filter center.
   *
   * @return Contribution sent to the owner.
   */
  template <int dim, typename MomentsType>
  FilterTargetContribution<dim>
  make_contribution(const unsigned int owner_index, const MomentsType &moments)
  {
    FilterTargetContribution<dim> contribution;
    contribution.owner_index     = owner_index;
    contribution.kernel_mass     = moments.kernel_mass;
    contribution.fluid_volume    = moments.fluid_volume;
    contribution.pressure_moment = moments.pressure_moment;
    contribution.boundary_weight = moments.boundary_weight;
    contribution.wall_weight     = moments.wall_weight;
    for (unsigned int d = 0; d < dim; ++d)
      {
        contribution.velocity_moment[d]      = moments.velocity_moment[d];
        contribution.wall_velocity_moment[d] = moments.wall_velocity_moment[d];
      }
    return contribution;
  }

  /**
   * @brief Unpack the partial moments of a filter center received from
   * another process.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @tparam MomentsType Type of the moments of the filter.
   *
   * @param[in] contribution Contribution received from another process.
   *
   * @return Partial moments of the filter center.
   */
  template <typename MomentsType, int dim>
  MomentsType
  unpack_contribution(const FilterTargetContribution<dim> &contribution)
  {
    MomentsType moments;
    moments.kernel_mass     = contribution.kernel_mass;
    moments.fluid_volume    = contribution.fluid_volume;
    moments.pressure_moment = contribution.pressure_moment;
    moments.boundary_weight = contribution.boundary_weight;
    moments.wall_weight     = contribution.wall_weight;
    for (unsigned int d = 0; d < dim; ++d)
      {
        moments.velocity_moment[d]      = contribution.velocity_moment[d];
        moments.wall_velocity_moment[d] = contribution.wall_velocity_moment[d];
      }
    return moments;
  }

  /**
   * @brief Return the squared distance between a point and an axis-aligned
   * box.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] point Point.
   *
   * @param[in] lower_corner Lower corner of the box.
   *
   * @param[in] upper_corner Upper corner of the box.
   *
   * @return Squared distance, which is zero if the point lies in the box.
   */
  template <int dim>
  double
  squared_distance_to_box(const Point<dim> &point,
                          const Point<dim> &lower_corner,
                          const Point<dim> &upper_corner)
  {
    double squared_distance = 0.;
    for (unsigned int d = 0; d < dim; ++d)
      {
        const double distance = std::max(
          {0., lower_corner[d] - point[d], point[d] - upper_corner[d]});
        squared_distance += distance * distance;
      }
    return squared_distance;
  }

  /**
   * @brief Check that a message can be sent by
   * Utilities::MPI::some_to_some(), which sends the number of objects followed
   * by the objects as a single buffer whose size is an int.
   *
   * @tparam T Type of the objects of the message.
   *
   * @param[in] message Objects to send.
   */
  template <typename T>
  void
  check_message_size(const std::vector<T> &message)
  {
    static_assert(std::is_trivially_copyable_v<T>,
                  "The messages of the filter are copied bit for bit.");
    const std::size_t message_size =
      sizeof(std::size_t) + message.size() * sizeof(T);
    AssertThrow(message_size <=
                  static_cast<std::size_t>(std::numeric_limits<int>::max()),
                ExcMessage(
                  "A message exchanged by the Anderson-Jackson filter exceeds "
                  "the 2 GB limit of MPI messages. Use more processes or a "
                  "smaller filter width."));
  }

  /**
   * @brief Describe the region covered by the locally owned cells with a few
   * bounding boxes. The cells are grouped by the position of their center in
   * a regular grid of bins covering the locally owned cells, and each bin is
   * described by the bounding box of its cells.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] triangulation Triangulation.
   *
   * @param[in] mapping Mapping of the triangulation.
   *
   * @return Bounding boxes of the locally owned cells. It is empty if the
   * process does not own any cell.
   */
  template <int dim>
  std::vector<BoundingBox<dim>>
  compute_local_source_boxes(const Triangulation<dim> &triangulation,
                             const Mapping<dim>       &mapping)
  {
    std::vector<BoundingBox<dim>> cell_boxes;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        if (!cell->is_locally_owned())
          continue;

        BoundingBox<dim> cell_box = mapping.get_bounding_box(cell);
        const double     margin   = cell_box_relative_margin *
                              cell_box.get_boundary_points().first.distance(
                                cell_box.get_boundary_points().second);
        cell_box.extend(margin);
        cell_boxes.push_back(cell_box);
      }

    if (cell_boxes.empty())
      return {};

    BoundingBox<dim> subdomain_box = cell_boxes.front();
    for (const BoundingBox<dim> &cell_box : cell_boxes)
      subdomain_box.merge_with(cell_box);
    const auto &[subdomain_lower_corner, subdomain_upper_corner] =
      subdomain_box.get_boundary_points();

    const unsigned int n_bins =
      Utilities::fixed_power<dim>(n_box_bins_per_direction);
    std::vector<BoundingBox<dim>> bin_boxes(n_bins);
    std::vector<bool>             bin_is_used(n_bins, false);
    for (const BoundingBox<dim> &cell_box : cell_boxes)
      {
        const Point<dim> center = cell_box.center();
        unsigned int     bin    = 0;
        unsigned int     stride = 1;
        for (unsigned int d = 0; d < dim; ++d)
          {
            const double relative_position =
              (center[d] - subdomain_lower_corner[d]) /
              (subdomain_upper_corner[d] - subdomain_lower_corner[d]);
            const unsigned int bin_in_direction =
              std::min(n_box_bins_per_direction - 1,
                       static_cast<unsigned int>(n_box_bins_per_direction *
                                                 relative_position));
            bin += bin_in_direction * stride;
            stride *= n_box_bins_per_direction;
          }

        if (bin_is_used[bin])
          bin_boxes[bin].merge_with(cell_box);
        else
          {
            bin_boxes[bin]   = cell_box;
            bin_is_used[bin] = true;
          }
      }

    std::vector<BoundingBox<dim>> source_boxes;
    for (unsigned int bin = 0; bin < n_bins; ++bin)
      if (bin_is_used[bin])
        source_boxes.push_back(bin_boxes[bin]);
    return source_boxes;
  }
} // namespace

template <int dim>
typename AndersonJacksonFilter<dim>::FilterMoments &
AndersonJacksonFilter<dim>::FilterMoments::operator+=(
  const FilterMoments &other)
{
  kernel_mass += other.kernel_mass;
  fluid_volume += other.fluid_volume;
  velocity_moment += other.velocity_moment;
  pressure_moment += other.pressure_moment;
  boundary_weight += other.boundary_weight;
  wall_weight += other.wall_weight;
  wall_velocity_moment += other.wall_velocity_moment;
  return *this;
}

template <int dim>
void
AndersonJacksonFilter<dim>::SourceBatch::clear()
{
  for (unsigned int d = 0; d < dim; ++d)
    {
      coordinates[d].clear();
      velocity_moments[d].clear();
    }
  weights.clear();
  fluid_weights.clear();
  pressure_moments.clear();
}

template <int dim>
AndersonJacksonFilter<dim>::AndersonJacksonFilter(
  const Parameters::AndersonJacksonFilter &filter_parameters,
  const Parameters::Timer                 &timer_parameters,
  const MPI_Comm                           mpi_communicator)
  : filter_parameters(filter_parameters)
  , timer_parameters(timer_parameters)
  , mpi_communicator(mpi_communicator)
  , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  , computing_timer(mpi_communicator,
                    pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
{}

template <int dim>
void
AndersonJacksonFilter<dim>::apply(
  const DoFHandler<dim>                &fluid_dof_handler,
  const Mapping<dim>                   &mapping,
  const GlobalVectorType               &fluid_solution,
  const ImmersedSolidClassifier<dim>   &classifier,
  const Parameters::PeriodicBoundaries &periodic_boundaries,
  const WallVelocities                 &wall_velocities)
{
  // The kernel is a template argument of the stages that evaluate it, so that
  // the kernel evaluation in the innermost loop is not a runtime dispatch.
  if (filter_parameters.kernel_type == Parameters::FilterKernelType::gaussian)
    {
      const GaussianFilterKernel<dim> kernel(filter_parameters.filter_width,
                                             filter_parameters.gaussian_cutoff);
      retained_mass_fraction = kernel.retained_mass_fraction();
      apply_kernel(kernel,
                   fluid_dof_handler,
                   mapping,
                   fluid_solution,
                   classifier,
                   periodic_boundaries,
                   wall_velocities);
    }
  else
    {
      const TopHatFilterKernel<dim> kernel(filter_parameters.filter_width);
      retained_mass_fraction = 1.;
      apply_kernel(kernel,
                   fluid_dof_handler,
                   mapping,
                   fluid_solution,
                   classifier,
                   periodic_boundaries,
                   wall_velocities);
    }
}

template <int dim>
template <typename KernelType>
void
AndersonJacksonFilter<dim>::apply_kernel(
  const KernelType                     &kernel,
  const DoFHandler<dim>                &fluid_dof_handler,
  const Mapping<dim>                   &mapping,
  const GlobalVectorType               &fluid_solution,
  const ImmersedSolidClassifier<dim>   &classifier,
  const Parameters::PeriodicBoundaries &periodic_boundaries,
  const WallVelocities                 &wall_velocities)
{
  support_radius = kernel.support_radius();

  {
    TimerOutput::Scope t(computing_timer, "Setup filter space");
    setup_filter_space(fluid_dof_handler);
  }

  {
    TimerOutput::Scope t(computing_timer, "Generate filter centers");
    generate_owned_targets(mapping);
  }

  {
    TimerOutput::Scope t(computing_timer, "Distribute filter centers");
    distribute_targets(mapping, periodic_boundaries);
  }

  TargetTree target_tree;
  {
    TimerOutput::Scope t(computing_timer, "Build spatial search");
    target_tree = build_target_tree();
  }

  {
    TimerOutput::Scope t(computing_timer, "Integrate source cells");
    integrate_source_cells(kernel,
                           target_tree,
                           fluid_dof_handler,
                           mapping,
                           fluid_solution,
                           classifier,
                           wall_velocities);
  }

  {
    TimerOutput::Scope t(computing_timer, "Return filter contributions");
    exchange_partial_integrals();
  }

  {
    TimerOutput::Scope t(computing_timer, "Normalize filter fields");
    normalize_filtered_fields();
    compute_statistics(mapping);
  }

  print_statistics();
}

template <int dim>
void
AndersonJacksonFilter<dim>::setup_filter_space(
  const DoFHandler<dim> &fluid_dof_handler)
{
  const FiniteElement<dim> &fluid_fe = fluid_dof_handler.get_fe();
  AssertThrow(fluid_fe.n_components() == dim + 1,
              ExcMessage("The Anderson-Jackson filter requires a "
                         "velocity-pressure finite element with dim velocity "
                         "components followed by the pressure."));

  const Triangulation<dim> &triangulation =
    fluid_dof_handler.get_triangulation();
  AssertThrow(triangulation.all_reference_cells_are_hyper_cube(),
              ExcMessage("The Anderson-Jackson filter only supports "
                         "quadrilateral and hexahedral meshes."));

  const unsigned int velocity_degree = fluid_fe.get_sub_fe(0, 1).degree;
  output_degree                      = (filter_parameters.output_degree > 0) ?
                                         filter_parameters.output_degree :
                                         velocity_degree;
  n_quadrature_points = (filter_parameters.n_quadrature_points > 0) ?
                          filter_parameters.n_quadrature_points :
                          velocity_degree + 1;

  fe = std::make_unique<FE_Q<dim>>(output_degree);
  dof_handler.reinit(triangulation);
  dof_handler.distribute_dofs(*fe);

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  // Only the hanging nodes are constrained. The degrees of freedom on the
  // periodic boundaries remain independent filter centers, which evaluate the
  // same filtered values.
  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  for (FieldVectorType *field : {&filtered_fields.fluid_volume_fraction,
                                 &filtered_fields.solid_volume_fraction,
                                 &filtered_fields.kernel_mass,
                                 &filtered_fields.pressure,
                                 &filtered_fields.valid})
    field->reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  for (FieldVectorType &velocity_component : filtered_fields.velocity)
    velocity_component.reinit(locally_owned_dofs,
                              locally_relevant_dofs,
                              mpi_communicator);
}

template <int dim>
void
AndersonJacksonFilter<dim>::generate_owned_targets(const Mapping<dim> &mapping)
{
  owned_target_locations.clear();
  owned_target_dofs.clear();

  const Quadrature<dim> support_point_quadrature(fe->get_unit_support_points());
  FEValues<dim>         fe_values(mapping,
                          *fe,
                          support_point_quadrature,
                          update_quadrature_points);

  std::vector<types::global_dof_index> local_dof_indices(fe->n_dofs_per_cell());
  std::vector<bool> is_visited(locally_owned_dofs.n_elements(), false);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
        {
          // The constrained degrees of freedom (hanging nodes) are not filter
          // centers: their values are interpolated from the constraints.
          const types::global_dof_index dof = local_dof_indices[i];
          if (!locally_owned_dofs.is_element(dof) ||
              constraints.is_constrained(dof))
            continue;

          const auto index = locally_owned_dofs.index_within_set(dof);
          if (is_visited[index])
            continue;
          is_visited[index] = true;

          owned_target_locations.push_back(fe_values.quadrature_point(i));
          owned_target_dofs.push_back(dof);
        }
    }
}

template <int dim>
void
AndersonJacksonFilter<dim>::distribute_targets(
  const Mapping<dim>                   &mapping,
  const Parameters::PeriodicBoundaries &periodic_boundaries)
{
  const Triangulation<dim> &triangulation = dof_handler.get_triangulation();
  const unsigned int        this_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);
  const double support_radius_squared = support_radius * support_radius;

  // Every process describes the region covered by its source points with a
  // few bounding boxes, which are gathered on every process. Processes that
  // do not own any cell contribute no box, but still take part in the
  // collective communications.
  const std::vector<BoundingBox<dim>> local_boxes =
    compute_local_source_boxes(triangulation, mapping);
  const RTree<std::pair<BoundingBox<dim>, unsigned int>> global_boxes =
    GridTools::build_global_description_tree(local_boxes, mpi_communicator);

  // The periodic translations are computed from the coarse mesh, so that they
  // are known on every process, including the processes that do not own
  // cells at the periodic boundaries.
  const std::vector<LetheGridTools::PeriodicTranslation<dim>> translations =
    LetheGridTools::compute_periodic_translations(triangulation,
                                                  periodic_boundaries);
  for (const auto &translation : translations)
    AssertThrow(
      2. * support_radius <=
        std::abs(translation.offset[translation.direction]),
      ExcMessage(
        "The support radius of the kernel of the Anderson-Jackson filter (" +
        std::to_string(support_radius) +
        ") exceeds half of the period of the domain in direction " +
        std::to_string(translation.direction) +
        ". A source point could then contribute to a filter center through "
        "several periodic images."));

  std::vector<std::pair<Point<dim>, unsigned int>>              local_targets;
  std::map<unsigned int, std::vector<FilterTargetRequest<dim>>> requests;
  std::vector<Point<dim>>                                       images;
  std::vector<unsigned int>                                     processes;
  for (unsigned int t = 0; t < owned_target_locations.size(); ++t)
    {
      const Point<dim> &location = owned_target_locations[t];

      // A filter center within the support radius of a periodic boundary also
      // gathers the sources located across that boundary, through its image
      // translated to the other side of the domain. Since the support radius
      // does not exceed half of the period, a filter center has at most one
      // image per periodic direction, and the images along several directions
      // are combined to reach the edges and corners of the domain.
      images.assign(1, location);
      for (const auto &translation : translations)
        {
          const unsigned int n_images   = images.size();
          const double       coordinate = location[translation.direction];
          if (std::abs(coordinate - translation.principal_coordinate) <
              support_radius)
            for (unsigned int i = 0; i < n_images; ++i)
              images.push_back(images[i] + translation.offset);
          else if (std::abs(coordinate - translation.neighbor_coordinate) <
                   support_radius)
            for (unsigned int i = 0; i < n_images; ++i)
              images.push_back(images[i] - translation.offset);
        }

      for (const Point<dim> &image : images)
        {
          // Processes owning cells within the support radius of the image.
          processes.clear();
          BoundingBox<dim> search_box(std::make_pair(image, image));
          search_box.extend(support_radius);
          global_boxes.query(
            boost::geometry::index::intersects(search_box),
            boost::make_function_output_iterator(
              [&](const std::pair<BoundingBox<dim>, unsigned int> &box) {
                const auto &[lower_corner, upper_corner] =
                  box.first.get_boundary_points();
                if (squared_distance_to_box(image, lower_corner, upper_corner) <
                    support_radius_squared)
                  processes.push_back(box.second);
              }));
          std::sort(processes.begin(), processes.end());
          processes.erase(std::unique(processes.begin(), processes.end()),
                          processes.end());

          for (const unsigned int process : processes)
            {
              if (process == this_process)
                local_targets.emplace_back(image, t);
              else
                {
                  FilterTargetRequest<dim> request;
                  for (unsigned int d = 0; d < dim; ++d)
                    request.location[d] = image[d];
                  request.owner_index = t;
                  requests[process].push_back(request);
                }
            }
        }
    }

  for (const auto &process_requests : requests)
    check_message_size(process_requests.second);
  const std::map<unsigned int, std::vector<FilterTargetRequest<dim>>>
    received_requests =
      Utilities::MPI::some_to_some(mpi_communicator, requests);

  // The filter centers integrated by this process are stored in blocks sharing
  // the same owner, starting with the filter centers owned by this process.
  integration_target_locations.clear();
  integration_target_owner_indices.clear();
  owner_blocks.clear();

  for (const auto &[location, owner_index] : local_targets)
    {
      integration_target_locations.push_back(location);
      integration_target_owner_indices.push_back(owner_index);
    }
  owner_blocks.push_back(
    OwnerBlock{this_process,
               0,
               static_cast<unsigned int>(integration_target_locations.size())});

  for (const auto &[process, process_requests] : received_requests)
    {
      const unsigned int begin = integration_target_locations.size();
      for (const FilterTargetRequest<dim> &request : process_requests)
        {
          Point<dim> location;
          for (unsigned int d = 0; d < dim; ++d)
            location[d] = request.location[d];
          integration_target_locations.push_back(location);
          integration_target_owner_indices.push_back(request.owner_index);
        }
      owner_blocks.push_back(OwnerBlock{
        process,
        begin,
        static_cast<unsigned int>(integration_target_locations.size())});
    }

  integration_target_moments.assign(integration_target_locations.size(),
                                    FilterMoments());
}

template <int dim>
typename AndersonJacksonFilter<dim>::TargetTree
AndersonJacksonFilter<dim>::build_target_tree()
{
  std::vector<std::pair<Point<dim>, unsigned int>> leaves;
  leaves.reserve(integration_target_locations.size());
  for (unsigned int i = 0; i < integration_target_locations.size(); ++i)
    leaves.emplace_back(integration_target_locations[i], i);

  // The R-tree stores a copy of the locations, which are released to limit
  // the peak memory when the halo of filter centers is large.
  std::vector<Point<dim>>().swap(integration_target_locations);

  return pack_rtree(leaves);
}

template <int dim>
template <typename KernelType>
void
AndersonJacksonFilter<dim>::integrate_source_cells(
  const KernelType                   &kernel,
  const TargetTree                   &target_tree,
  const DoFHandler<dim>              &fluid_dof_handler,
  const Mapping<dim>                 &mapping,
  const GlobalVectorType             &fluid_solution,
  const ImmersedSolidClassifier<dim> &classifier,
  const WallVelocities               &wall_velocities)
{
  const FiniteElement<dim> &fluid_fe = fluid_dof_handler.get_fe();

  // The cells cut by an immersed solid are integrated with an iterated Gauss
  // quadrature, on which the discontinuous fluid indicator is evaluated point
  // by point.
  const QGauss<dim>    regular_quadrature(n_quadrature_points);
  const QIterated<dim> cut_quadrature(QGauss<1>(n_quadrature_points),
                                      filter_parameters.cut_cell_subdivisions);
  const UpdateFlags    update_flags =
    update_values | update_quadrature_points | update_JxW_values;
  FEValues<dim> regular_fe_values(mapping,
                                  fluid_fe,
                                  regular_quadrature,
                                  update_flags);
  FEValues<dim> cut_fe_values(mapping, fluid_fe, cut_quadrature, update_flags);

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  // The boundary faces are integrated to extend the velocity beyond the walls.
  const bool extend_velocity_beyond_walls =
    filter_parameters.extend_velocity_beyond_walls && !wall_velocities.empty();
  FEFaceValues<dim>           face_fe_values(mapping,
                                   fluid_fe,
                                   QGauss<dim - 1>(n_quadrature_points),
                                   update_quadrature_points |
                                     update_JxW_values);
  std::vector<Tensor<1, dim>> wall_velocity_values;

  // The values are zero-initialized, so that the moments of the solid cells,
  // whose values are not evaluated, are zero.
  std::vector<Tensor<1, dim>> velocity_values;
  std::vector<double>         pressure_values;
  std::vector<double>         fluid_indicator;
  std::vector<unsigned int>   cutting_solids;

  SourceBatch        batch;
  const unsigned int batch_capacity =
    maximum_batch_size +
    std::max(regular_quadrature.size(), cut_quadrature.size());
  for (unsigned int d = 0; d < dim; ++d)
    {
      batch.coordinates[d].reserve(batch_capacity);
      batch.velocity_moments[d].reserve(batch_capacity);
    }
  batch.weights.reserve(batch_capacity);
  batch.fluid_weights.reserve(batch_capacity);
  batch.pressure_moments.reserve(batch_capacity);

  const double maximum_batch_extent =
    maximum_batch_extent_fraction * kernel.support_radius();

  n_local_fluid_cells       = 0;
  n_local_solid_cells       = 0;
  n_local_cut_cells         = 0;
  local_source_volume       = 0.;
  local_source_solid_volume = 0.;

  // The locally owned cells are swept in the order of the triangulation,
  // which follows a space-filling curve, so that consecutive cells, and thus
  // the cells of a batch, are close to each other.
  for (const auto &cell : fluid_dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      const CellPhase phase =
        classifier.classify_cell(cell, mapping, cutting_solids);
      FEValues<dim> &fe_values =
        (phase == CellPhase::cut) ? cut_fe_values : regular_fe_values;
      fe_values.reinit(cell);

      const unsigned int             n_q_points = fe_values.n_quadrature_points;
      const std::vector<Point<dim>> &quadrature_points =
        fe_values.get_quadrature_points();
      const std::vector<double> &JxW = fe_values.get_JxW_values();

      velocity_values.resize(n_q_points);
      pressure_values.resize(n_q_points);
      if (phase == CellPhase::fluid)
        {
          ++n_local_fluid_cells;
          fluid_indicator.assign(n_q_points, 1.);
        }
      else if (phase == CellPhase::solid)
        {
          ++n_local_solid_cells;
          fluid_indicator.assign(n_q_points, 0.);
        }
      else
        {
          ++n_local_cut_cells;
          classifier.fill_fluid_indicator(cutting_solids,
                                          quadrature_points,
                                          fluid_indicator);
        }

      // The solid cells only contribute to the kernel mass.
      if (phase != CellPhase::solid)
        {
          fe_values[velocities].get_function_values(fluid_solution,
                                                    velocity_values);
          if (filter_parameters.filter_pressure)
            fe_values[pressure].get_function_values(fluid_solution,
                                                    pressure_values);
        }

      // Bounding box of the quadrature points of the cell.
      Point<dim> cell_lower_corner = quadrature_points[0];
      Point<dim> cell_upper_corner = quadrature_points[0];
      for (const Point<dim> &point : quadrature_points)
        for (unsigned int d = 0; d < dim; ++d)
          {
            cell_lower_corner[d] = std::min(cell_lower_corner[d], point[d]);
            cell_upper_corner[d] = std::max(cell_upper_corner[d], point[d]);
          }

      // The batch is processed before adding the cell if the cell would make
      // it too large or too wide.
      if (batch.size() > 0)
        {
          double merged_extent = 0.;
          for (unsigned int d = 0; d < dim; ++d)
            merged_extent =
              std::max(merged_extent,
                       std::max(batch.upper_corner[d], cell_upper_corner[d]) -
                         std::min(batch.lower_corner[d], cell_lower_corner[d]));
          if (batch.size() + n_q_points > maximum_batch_size ||
              merged_extent > maximum_batch_extent)
            accumulate_batch(kernel, target_tree, batch);
        }

      if (batch.size() == 0)
        {
          batch.lower_corner = cell_lower_corner;
          batch.upper_corner = cell_upper_corner;
        }
      else
        for (unsigned int d = 0; d < dim; ++d)
          {
            batch.lower_corner[d] =
              std::min(batch.lower_corner[d], cell_lower_corner[d]);
            batch.upper_corner[d] =
              std::max(batch.upper_corner[d], cell_upper_corner[d]);
          }

      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          const double fluid_weight = fluid_indicator[q] * JxW[q];
          for (unsigned int d = 0; d < dim; ++d)
            {
              batch.coordinates[d].push_back(quadrature_points[q][d]);
              batch.velocity_moments[d].push_back(fluid_weight *
                                                  velocity_values[q][d]);
            }
          batch.weights.push_back(JxW[q]);
          batch.fluid_weights.push_back(fluid_weight);
          batch.pressure_moments.push_back(fluid_weight * pressure_values[q]);

          local_source_volume += JxW[q];
          local_source_solid_volume += JxW[q] - fluid_weight;
        }

      if (!extend_velocity_beyond_walls || !cell->at_boundary())
        continue;

      // The faces on the periodic boundaries are interior faces of the
      // periodic domain and do not truncate the kernel.
      for (const unsigned int f : cell->face_indices())
        {
          if (!cell->at_boundary(f) || cell->has_periodic_neighbor(f))
            continue;

          face_fe_values.reinit(cell, f);
          const std::vector<Point<dim>> &face_points =
            face_fe_values.get_quadrature_points();

          wall_velocity_values.clear();
          const auto wall = wall_velocities.find(cell->face(f)->boundary_id());
          if (wall != wall_velocities.end())
            {
              wall_velocity_values.resize(face_points.size());
              for (unsigned int q = 0; q < face_points.size(); ++q)
                for (unsigned int d = 0; d < dim; ++d)
                  wall_velocity_values[q][d] =
                    wall->second->value(face_points[q], d);
            }

          accumulate_boundary_face(kernel,
                                   target_tree,
                                   face_points,
                                   face_fe_values.get_JxW_values(),
                                   wall_velocity_values);
        }
    }

  accumulate_batch(kernel, target_tree, batch);
}

template <int dim>
template <typename KernelType>
void
AndersonJacksonFilter<dim>::accumulate_boundary_face(
  const KernelType                  &kernel,
  const TargetTree                  &target_tree,
  const std::vector<Point<dim>>     &points,
  const std::vector<double>         &weights,
  const std::vector<Tensor<1, dim>> &wall_velocity_values)
{
  const bool   is_wall                = !wall_velocity_values.empty();
  const double support_radius_squared = kernel.support_radius_squared();

  Point<dim> lower_corner = points[0];
  Point<dim> upper_corner = points[0];
  for (const Point<dim> &point : points)
    for (unsigned int d = 0; d < dim; ++d)
      {
        lower_corner[d] = std::min(lower_corner[d], point[d]);
        upper_corner[d] = std::max(upper_corner[d], point[d]);
      }

  // Accumulate the integrals of the kernel over the face at a filter center.
  const auto accumulate_target =
    [&](const std::pair<Point<dim>, unsigned int> &target) {
      const Point<dim> &center = target.first;
      if (squared_distance_to_box(center, lower_corner, upper_corner) >=
          support_radius_squared)
        return;

      double         boundary_weight = 0.;
      Tensor<1, dim> wall_velocity_moment;
      for (unsigned int q = 0; q < points.size(); ++q)
        {
          const double weight = kernel.value_from_squared_distance(
                                  center.distance_square(points[q])) *
                                weights[q];
          boundary_weight += weight;
          if (is_wall)
            wall_velocity_moment += weight * wall_velocity_values[q];
        }

      FilterMoments &moments = integration_target_moments[target.second];
      moments.boundary_weight += boundary_weight;
      if (is_wall)
        {
          moments.wall_weight += boundary_weight;
          moments.wall_velocity_moment += wall_velocity_moment;
        }
    };

  BoundingBox<dim> search_box(std::make_pair(lower_corner, upper_corner));
  search_box.extend(kernel.support_radius());
  target_tree.query(boost::geometry::index::intersects(search_box),
                    boost::make_function_output_iterator(accumulate_target));
}

template <int dim>
template <typename KernelType>
void
AndersonJacksonFilter<dim>::accumulate_batch(const KernelType &kernel,
                                             const TargetTree &target_tree,
                                             SourceBatch      &batch)
{
  const unsigned int n_points = batch.size();
  if (n_points == 0)
    return;

  const double support_radius_squared = kernel.support_radius_squared();

  std::array<const double *, dim> coordinates;
  std::array<const double *, dim> velocity_moments;
  for (unsigned int d = 0; d < dim; ++d)
    {
      coordinates[d]      = batch.coordinates[d].data();
      velocity_moments[d] = batch.velocity_moments[d].data();
    }
  const double *weights          = batch.weights.data();
  const double *fluid_weights    = batch.fluid_weights.data();
  const double *pressure_moments = batch.pressure_moments.data();

  // Accumulate the contribution of every point of the batch to a filter
  // center.
  const auto accumulate_target =
    [&](const std::pair<Point<dim>, unsigned int> &target) {
      const Point<dim> &center = target.first;

      // The search box also contains filter centers near its corners, which
      // lie beyond the support radius of every point of the batch.
      if (squared_distance_to_box(center,
                                  batch.lower_corner,
                                  batch.upper_corner) >= support_radius_squared)
        return;

      double         kernel_mass     = 0.;
      double         fluid_volume    = 0.;
      double         pressure_moment = 0.;
      Tensor<1, dim> velocity_moment;
      for (unsigned int q = 0; q < n_points; ++q)
        {
          double squared_distance = 0.;
          for (unsigned int d = 0; d < dim; ++d)
            {
              const double difference = center[d] - coordinates[d][q];
              squared_distance += difference * difference;
            }
          const double g = kernel.value_from_squared_distance(squared_distance);
          kernel_mass += g * weights[q];
          fluid_volume += g * fluid_weights[q];
          for (unsigned int d = 0; d < dim; ++d)
            velocity_moment[d] += g * velocity_moments[d][q];
          pressure_moment += g * pressure_moments[q];
        }

      FilterMoments &moments = integration_target_moments[target.second];
      moments.kernel_mass += kernel_mass;
      moments.fluid_volume += fluid_volume;
      moments.velocity_moment += velocity_moment;
      moments.pressure_moment += pressure_moment;
    };

  BoundingBox<dim> search_box(
    std::make_pair(batch.lower_corner, batch.upper_corner));
  search_box.extend(kernel.support_radius());
  target_tree.query(boost::geometry::index::intersects(search_box),
                    boost::make_function_output_iterator(accumulate_target));

  batch.clear();
}

template <int dim>
void
AndersonJacksonFilter<dim>::exchange_partial_integrals()
{
  const unsigned int this_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);

  owned_target_moments.assign(owned_target_locations.size(), FilterMoments());

  std::map<unsigned int, std::vector<FilterTargetContribution<dim>>>
    contributions;
  for (const OwnerBlock &block : owner_blocks)
    for (unsigned int i = block.begin; i < block.end; ++i)
      {
        // A filter center that received no contribution lies beyond the
        // support of the kernel of every source point of this process.
        const FilterMoments &moments = integration_target_moments[i];
        if (moments.kernel_mass <= 0. && moments.boundary_weight <= 0.)
          continue;

        const unsigned int owner_index = integration_target_owner_indices[i];
        if (block.process == this_process)
          owned_target_moments[owner_index] += moments;
        else
          contributions[block.process].push_back(
            make_contribution<dim>(owner_index, moments));
      }

  for (const auto &process_contributions : contributions)
    check_message_size(process_contributions.second);
  const std::map<unsigned int, std::vector<FilterTargetContribution<dim>>>
    received_contributions =
      Utilities::MPI::some_to_some(mpi_communicator, contributions);

  // The contributions are summed in the order of the ranks of the processes,
  // so that the result does not depend on the order of arrival of the
  // messages.
  for (const auto &process_contributions : received_contributions)
    for (const FilterTargetContribution<dim> &contribution :
         process_contributions.second)
      {
        AssertIndexRange(contribution.owner_index, owned_target_moments.size());
        owned_target_moments[contribution.owner_index] +=
          unpack_contribution<FilterMoments>(contribution);
      }
}

template <int dim>
void
AndersonJacksonFilter<dim>::normalize_filtered_fields()
{
  // Setting the vectors to zero also zeroes their ghost values, which keeps
  // them consistent while the owned values are written.
  filtered_fields.fluid_volume_fraction = 0.;
  filtered_fields.solid_volume_fraction = 0.;
  filtered_fields.kernel_mass           = 0.;
  filtered_fields.pressure              = 0.;
  filtered_fields.valid                 = 0.;
  for (FieldVectorType &velocity_component : filtered_fields.velocity)
    velocity_component = 0.;

  local_minimum_kernel_mass           = std::numeric_limits<double>::max();
  local_maximum_kernel_mass           = std::numeric_limits<double>::lowest();
  local_minimum_fluid_volume_fraction = std::numeric_limits<double>::max();
  local_maximum_fluid_volume_fraction = std::numeric_limits<double>::lowest();
  n_local_undefined_averages          = 0;

  for (unsigned int t = 0; t < owned_target_dofs.size(); ++t)
    {
      const FilterMoments          &moments      = owned_target_moments[t];
      const types::global_dof_index dof          = owned_target_dofs[t];
      const double                  kernel_mass  = moments.kernel_mass;
      const double                  fluid_volume = moments.fluid_volume;

      // The fluid volume never exceeds the kernel mass, since both integrate
      // the same non-negative terms and the fluid volume omits some of them.
      const bool has_mass              = kernel_mass > 0.;
      double     fluid_volume_fraction = fluid_volume;
      double     solid_volume_fraction = kernel_mass - fluid_volume;
      if (filter_parameters.normalize_at_domain_boundaries)
        {
          fluid_volume_fraction = has_mass ? fluid_volume / kernel_mass : 0.;
          solid_volume_fraction =
            has_mass ? (kernel_mass - fluid_volume) / kernel_mass : 0.;
        }

      filtered_fields.kernel_mass(dof)           = kernel_mass;
      filtered_fields.fluid_volume_fraction(dof) = fluid_volume_fraction;
      filtered_fields.solid_volume_fraction(dof) = solid_volume_fraction;

      // The phase averages are only defined where the kernel contains enough
      // fluid, relative to its mass.
      const bool is_defined =
        has_mass && fluid_volume > 0. &&
        fluid_volume >= filter_parameters.minimum_fluid_fraction * kernel_mass;
      if (is_defined)
        {
          // The part of the kernel outside of the domain, of mass 1 - M, is
          // attributed to the walls in proportion to the integral of the
          // kernel over the walls and over all the non-periodic boundaries,
          // and it is filled with fluid moving at the kernel-weighted average
          // of the wall velocity. The wall integrals are only accumulated if
          // the velocity is extended beyond the walls.
          double         velocity_weight = fluid_volume;
          Tensor<1, dim> velocity_moment = moments.velocity_moment;
          if (moments.wall_weight > 0.)
            {
              const double exterior_mass = std::max(0., 1. - kernel_mass) *
                                           moments.wall_weight /
                                           moments.boundary_weight;
              velocity_moment += (exterior_mass / moments.wall_weight) *
                                 moments.wall_velocity_moment;
              velocity_weight += exterior_mass;
            }

          for (unsigned int d = 0; d < dim; ++d)
            filtered_fields.velocity[d](dof) =
              velocity_moment[d] / velocity_weight;
          if (filter_parameters.filter_pressure)
            filtered_fields.pressure(dof) =
              moments.pressure_moment / fluid_volume;
          filtered_fields.valid(dof) = 1.;
        }
      else
        ++n_local_undefined_averages;

      local_minimum_kernel_mass =
        std::min(local_minimum_kernel_mass, kernel_mass);
      local_maximum_kernel_mass =
        std::max(local_maximum_kernel_mass, kernel_mass);
      local_minimum_fluid_volume_fraction =
        std::min(local_minimum_fluid_volume_fraction, fluid_volume_fraction);
      local_maximum_fluid_volume_fraction =
        std::max(local_maximum_fluid_volume_fraction, fluid_volume_fraction);
    }

  // The values of the hanging nodes are interpolated from the constraints.
  for (FieldVectorType *field : {&filtered_fields.fluid_volume_fraction,
                                 &filtered_fields.solid_volume_fraction,
                                 &filtered_fields.kernel_mass,
                                 &filtered_fields.pressure,
                                 &filtered_fields.valid})
    {
      constraints.distribute(*field);
      field->update_ghost_values();
    }
  for (FieldVectorType &velocity_component : filtered_fields.velocity)
    {
      constraints.distribute(velocity_component);
      velocity_component.update_ghost_values();
    }
}

template <int dim>
void
AndersonJacksonFilter<dim>::compute_statistics(const Mapping<dim> &mapping)
{
  statistics.n_filter_centers =
    Utilities::MPI::sum(static_cast<types::global_dof_index>(
                          owned_target_dofs.size()),
                        mpi_communicator);
  statistics.n_fluid_cells =
    Utilities::MPI::sum(n_local_fluid_cells, mpi_communicator);
  statistics.n_solid_cells =
    Utilities::MPI::sum(n_local_solid_cells, mpi_communicator);
  statistics.n_cut_cells =
    Utilities::MPI::sum(n_local_cut_cells, mpi_communicator);
  statistics.n_undefined_averages =
    Utilities::MPI::sum(n_local_undefined_averages, mpi_communicator);
  statistics.minimum_kernel_mass =
    Utilities::MPI::min(local_minimum_kernel_mass, mpi_communicator);
  statistics.maximum_kernel_mass =
    Utilities::MPI::max(local_maximum_kernel_mass, mpi_communicator);
  statistics.minimum_fluid_volume_fraction =
    Utilities::MPI::min(local_minimum_fluid_volume_fraction, mpi_communicator);
  statistics.maximum_fluid_volume_fraction =
    Utilities::MPI::max(local_maximum_fluid_volume_fraction, mpi_communicator);
  statistics.source_volume =
    Utilities::MPI::sum(local_source_volume, mpi_communicator);
  statistics.source_solid_volume =
    Utilities::MPI::sum(local_source_solid_volume, mpi_communicator);

  // Integral of the solid volume fraction. For a periodic domain without
  // renormalization, it equals the volume of the solids, since the kernel has
  // unit mass, which verifies that every source point reached every filter
  // center within its support.
  const QGauss<dim>   quadrature(output_degree + 1);
  FEValues<dim>       fe_values(mapping,
                          *fe,
                          quadrature,
                          update_values | update_JxW_values);
  std::vector<double> solid_volume_fraction_values(quadrature.size());
  double              local_integral = 0.;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      fe_values.get_function_values(filtered_fields.solid_volume_fraction,
                                    solid_volume_fraction_values);
      for (const unsigned int q : fe_values.quadrature_point_indices())
        local_integral += solid_volume_fraction_values[q] * fe_values.JxW(q);
    }
  statistics.solid_volume_fraction_integral =
    Utilities::MPI::sum(local_integral, mpi_communicator);
}

template <int dim>
void
AndersonJacksonFilter<dim>::print_statistics() const
{
  if (filter_parameters.verbosity == Parameters::Verbosity::quiet)
    return;

  announce_string(pcout, "Anderson-Jackson filter");
  pcout << std::scientific << std::setprecision(statistics_precision);
  if (filter_parameters.kernel_type == Parameters::FilterKernelType::gaussian)
    {
      pcout << statistics_label("Kernel:") << "gaussian" << std::endl;
      pcout << statistics_label("Standard deviation:")
            << filter_parameters.filter_width << std::endl;
      pcout << statistics_label("Support radius:") << support_radius
            << std::endl;
      pcout << statistics_label("Retained mass of the gaussian:")
            << retained_mass_fraction << std::endl;
    }
  else
    {
      pcout << statistics_label("Kernel:") << "top-hat" << std::endl;
      pcout << statistics_label("Support radius:") << support_radius
            << std::endl;
    }
  pcout << statistics_label("Polynomial degree of the filtered fields:")
        << output_degree << std::endl;
  pcout << statistics_label("Quadrature points per direction:")
        << n_quadrature_points << std::endl;
  pcout << statistics_label("Cut cell subdivisions:")
        << filter_parameters.cut_cell_subdivisions << std::endl;
  pcout << statistics_label("Number of filter centers:")
        << statistics.n_filter_centers << std::endl;
  pcout << statistics_label("Number of fluid, solid and cut cells:")
        << statistics.n_fluid_cells << " " << statistics.n_solid_cells << " "
        << statistics.n_cut_cells << std::endl;
  pcout << statistics_label("Kernel mass (min, max):")
        << statistics.minimum_kernel_mass << " "
        << statistics.maximum_kernel_mass << std::endl;
  pcout << statistics_label("Fluid volume fraction (min, max):")
        << statistics.minimum_fluid_volume_fraction << " "
        << statistics.maximum_fluid_volume_fraction << std::endl;
  pcout << statistics_label("Number of undefined phase averages:")
        << statistics.n_undefined_averages << std::endl;
  pcout << statistics_label("Volume of the source cells:")
        << statistics.source_volume << std::endl;
  pcout << statistics_label("Solid volume of the source cells:")
        << statistics.source_solid_volume << std::endl;
  pcout << statistics_label("Integral of the solid volume fraction:")
        << statistics.solid_volume_fraction_integral << std::endl;

  if (filter_parameters.verbosity == Parameters::Verbosity::extra_verbose)
    {
      // Distribution of the work among the processes: the owned filter
      // centers, and the integrated filter centers, which include the
      // periodic images and the halo of filter centers of the other
      // processes.
      const Utilities::MPI::MinMaxAvg owned_statistics =
        Utilities::MPI::min_max_avg(
          static_cast<double>(owned_target_dofs.size()), mpi_communicator);
      const Utilities::MPI::MinMaxAvg integrated_statistics =
        Utilities::MPI::min_max_avg(static_cast<double>(
                                      integration_target_moments.size()),
                                    mpi_communicator);
      pcout << std::defaultfloat;
      pcout << statistics_label("Owned filter centers per process:") << "min "
            << owned_statistics.min << ", max " << owned_statistics.max
            << ", average " << owned_statistics.avg << std::endl;
      pcout << statistics_label("Integrated filter centers per process:")
            << "min " << integrated_statistics.min << ", max "
            << integrated_statistics.max << ", average "
            << integrated_statistics.avg << std::endl;
      pcout << std::scientific;
    }
}

template <int dim>
void
AndersonJacksonFilter<dim>::write_output(const Mapping<dim> &mapping,
                                         const double        time,
                                         const unsigned int  group_files)
{
  TimerOutput::Scope t(computing_timer, "Write output");

  // The velocity is copied to a vector-valued space so that it is visualized
  // as a vector.
  const FESystem<dim> velocity_fe(*fe, dim);
  DoFHandler<dim>     velocity_dof_handler(dof_handler.get_triangulation());
  velocity_dof_handler.distribute_dofs(velocity_fe);
  const IndexSet velocity_owned_dofs =
    velocity_dof_handler.locally_owned_dofs();
  FieldVectorType velocity(velocity_owned_dofs,
                           DoFTools::extract_locally_relevant_dofs(
                             velocity_dof_handler),
                           mpi_communicator);

  std::vector<types::global_dof_index> scalar_dof_indices(
    fe->n_dofs_per_cell());
  std::vector<types::global_dof_index> vector_dof_indices(
    velocity_fe.n_dofs_per_cell());
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      cell->get_dof_indices(scalar_dof_indices);
      cell->as_dof_handler_iterator(velocity_dof_handler)
        ->get_dof_indices(vector_dof_indices);
      for (unsigned int i = 0; i < scalar_dof_indices.size(); ++i)
        for (unsigned int d = 0; d < dim; ++d)
          {
            const types::global_dof_index vector_dof =
              vector_dof_indices[velocity_fe.component_to_system_index(d, i)];
            if (velocity_owned_dofs.is_element(vector_dof))
              velocity(vector_dof) =
                filtered_fields.velocity[d](scalar_dof_indices[i]);
          }
    }
  velocity.update_ghost_values();

  // The gradient of the filtered velocity is the gradient of its finite
  // element interpolant. The postprocessor is declared before the DataOut
  // object, which refers to it until it is destroyed.
  const GradientPostprocessor<dim> velocity_gradient(
    "filtered_velocity_gradient");

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(filtered_fields.fluid_volume_fraction,
                           "fluid_volume_fraction");
  data_out.add_data_vector(filtered_fields.solid_volume_fraction,
                           "solid_volume_fraction");
  data_out.add_data_vector(filtered_fields.kernel_mass, "kernel_mass");
  if (filter_parameters.filter_pressure)
    data_out.add_data_vector(filtered_fields.pressure, "filtered_pressure");
  data_out.add_data_vector(filtered_fields.valid, "valid");
  data_out.add_data_vector(
    velocity_dof_handler,
    velocity,
    std::vector<std::string>(dim, "filtered_velocity"),
    std::vector<DataComponentInterpretation::DataComponentInterpretation>(
      dim, DataComponentInterpretation::component_is_part_of_vector));
  data_out.add_data_vector(velocity_dof_handler, velocity, velocity_gradient);
  data_out.build_patches(mapping,
                         output_degree,
                         DataOut<dim>::curved_inner_cells);

  // The output folder is created by the first process before any process
  // writes in it.
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    create_output_folder(filter_parameters.output_folder);
  const int ierr = MPI_Barrier(mpi_communicator);
  AssertThrowMPI(ierr);

  write_vtu_and_pvd<dim>(pvd_handler,
                         data_out,
                         filter_parameters.output_folder,
                         filter_parameters.output_name,
                         time,
                         n_outputs,
                         group_files,
                         mpi_communicator);
  ++n_outputs;
}

template <int dim>
void
AndersonJacksonFilter<dim>::print_timer_summary() const
{
  if (timer_parameters.type == Parameters::Timer::Type::none)
    return;

  announce_string(pcout, "Anderson-Jackson filter");
  pcout << std::defaultfloat;
  computing_timer.print_summary();
  pcout << std::scientific;
}

template <int dim>
typename AndersonJacksonFilter<dim>::WallVelocities
make_wall_velocities(
  BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const double                                   time)
{
  typename AndersonJacksonFilter<dim>::WallVelocities wall_velocities;
  for (const auto &[id, type] : boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::noslip)
        wall_velocities[id] =
          std::make_shared<Functions::ZeroFunction<dim>>(dim);
      else if (type == BoundaryConditions::BoundaryType::function ||
               type == BoundaryConditions::BoundaryType::function_weak)
        {
          BoundaryConditions::NSBoundaryFunctions<dim> &functions =
            *boundary_conditions.navier_stokes_functions.at(id);
          functions.u.set_time(time);
          functions.v.set_time(time);
          functions.w.set_time(time);
          wall_velocities[id] =
            std::make_shared<NavierStokesFunctionDefined<dim>>(&functions.u,
                                                               &functions.v,
                                                               &functions.w);
        }
    }
  return wall_velocities;
}

template class AndersonJacksonFilter<2>;
template class AndersonJacksonFilter<3>;

template AndersonJacksonFilter<2>::WallVelocities
make_wall_velocities(
  BoundaryConditions::NSBoundaryConditions<2> &boundary_conditions,
  const double                                 time);
template AndersonJacksonFilter<3>::WallVelocities
make_wall_velocities(
  BoundaryConditions::NSBoundaryConditions<3> &boundary_conditions,
  const double                                 time);
