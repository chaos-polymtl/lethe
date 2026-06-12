#include <core/manifolds.h>
#include <core/pvd_handler.h>
#include <core/simulation_control.h>
#include <core/solutions_output.h>

#include <solvers/fluid_dynamics_matrix_based.h>
#include <solvers/solid_phase.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>


// using namespace dealii;



class SolidBlockDiagonalPreconditioner : public Subscriptor
{
public:
  SolidBlockDiagonalPreconditioner(
    const std::shared_ptr<TrilinosWrappers::PreconditionAMG> &velocity_amg,
    const std::shared_ptr<TrilinosWrappers::PreconditionILU> &velocity_ilu,
    const std::shared_ptr<TrilinosWrappers::PreconditionILU> &alpha_ilu,
    const bool                                                use_amg)
    : velocity_amg(velocity_amg)
    , velocity_ilu(velocity_ilu)
    , alpha_ilu(alpha_ilu)
    , use_amg(use_amg)
  {}

  void
  vmult(TrilinosWrappers::MPI::BlockVector       &dst,
        const TrilinosWrappers::MPI::BlockVector &src) const
  {
    if (dst.block(0).size() > 0)
      {
        if (use_amg && velocity_amg)
          velocity_amg->vmult(dst.block(0), src.block(0));
        else if (velocity_ilu)
          velocity_ilu->vmult(dst.block(0), src.block(0));
        else
          dst.block(0) = src.block(0);
      }

    if (dst.block(1).size() > 0)
      {
        if (alpha_ilu)
          alpha_ilu->vmult(dst.block(1), src.block(1));
        else
          dst.block(1) = src.block(1);
      }
  }

private:
  std::shared_ptr<TrilinosWrappers::PreconditionAMG> velocity_amg;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> velocity_ilu;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> alpha_ilu;
  bool                                               use_amg;
};

template <int dim>
void
SolidPhaseSolver<dim>::setup_linear_preconditioners()
{
  velocity_amg.reset();
  velocity_ilu.reset();
  alpha_ilu.reset();

  const bool use_amg = (parameters.amg_sweeps > 0);

  if (use_amg)
    {
      velocity_amg = std::make_shared<TrilinosWrappers::PreconditionAMG>();

      TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
      amg_data.elliptic              = parameters.amg_elliptic;
      amg_data.higher_order_elements = (order > 1);
      amg_data.smoother_type         = parameters.amg_smoother_type.c_str();
      amg_data.smoother_sweeps       = parameters.amg_sweeps;
      amg_data.aggregation_threshold = parameters.amg_agg_threshold;
      amg_data.output_details        = false;

      velocity_amg->initialize(system_matrix.block(0, 0), amg_data);
    }
  else
    {
      velocity_ilu = std::make_shared<TrilinosWrappers::PreconditionILU>();

      TrilinosWrappers::PreconditionILU::AdditionalData ilu_data;
      ilu_data.overlap = parameters.ilu_overlap;

      velocity_ilu->initialize(system_matrix.block(0, 0), ilu_data);
    }

  alpha_ilu = std::make_shared<TrilinosWrappers::PreconditionILU>();

  TrilinosWrappers::PreconditionILU::AdditionalData ilu_data;
  ilu_data.overlap = parameters.ilu_overlap;
  alpha_ilu->initialize(system_matrix.block(1, 1), ilu_data);
}


// template <int dim>
// class SolidInitialValues : public Function<dim>
// {
// public:
//   SolidInitialValues(const double &alpha0, const double &u0)
//     : Function<dim>(dim + 1)
//     , alpha_0(alpha0)
//     , u_0(u0)
//   {}


//   virtual double
//   value(const Point<dim> & /*p*/, const unsigned int component) const
//   override
//   {
//     if (component < dim)
//       {
//         return component == 0 ? u_0 : 0.0;
//       }
//     else
//       {
//         return alpha_0;
//       }
//   }


// private:
//   const double alpha_0;
//   const double u_0;
// };



// template <int dim>
// class SolidBoundaryValues : public Function<dim>
// {
// public:
//   SolidBoundaryValues(const double alpha_inlet, const Tensor<1, dim>
//   &u_inlet)
//     : Function<dim>(dim + 1)
//     , alpha_in(alpha_inlet)
//     , u_in(u_inlet)
//   {}

//   double
//   value(const Point<dim> &, const unsigned int component) const override
//   {
//     if (component < dim)
//       {
//         return u_in[component];
//       }
//     else
//       {
//         return alpha_in;
//       }
//   }

// private:
//   const double         alpha_in;
//   const Tensor<1, dim> u_in;
// };

template <int dim>
class SolidVelocityFunctionDefined : public Function<dim>
{
public:
  SolidVelocityFunctionDefined(Function<dim> *u,
                               Function<dim> *v,
                               Function<dim> *w = nullptr)
    : Function<dim>(dim + 1)
    , u(u)
    , v(v)
    , w(w)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const override
  {
    if (component == 0)
      return u->value(p);
    if (component == 1 && dim > 1)
      return v->value(p);
    if (component == 2 && dim > 2)
      return w->value(p);
    return 0.0;
  }

private:
  Function<dim> *u;
  Function<dim> *v;
  Function<dim> *w;
};

template <int dim>
class SolidAlphaFunctionDefined : public Function<dim>
{
public:
  explicit SolidAlphaFunctionDefined(Function<dim> *alpha_function)
    : Function<dim>(dim + 1)
    , alpha_function(alpha_function)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const override
  {
    if (component == dim)
      return alpha_function->value(p);
    return 0.0;
  }

private:
  Function<dim> *alpha_function;
};



template <int dim>
SolidPhaseSolver<dim>::SolidPhaseSolver(
  const SolidPhaseParameters<dim>           &p,
  const std::shared_ptr<SimulationControl>  &simulation_control_in,
  parallel::distributed::Triangulation<dim> &tria,
  MPI_Comm                                   comm)
  : parameters(p)
  , mpi_communicator(comm)
  , simulation_control(simulation_control_in)
  , order(parameters.order)
  , fe(FE_Q<dim>(order), dim, FE_Q<dim>(order), 1)
  , triangulation(tria)
  , dof_handler(triangulation)
  , rho_s(parameters.rho_s)
  , beta(parameters.beta)
  , timestep_number(0)
  , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  , computing_timer(mpi_communicator,
                    pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
  , output_every(parameters.output_every)
  , digits(parameters.digits)
  , output_folder(parameters.output_folder)
  , output_prefix(parameters.output_prefix)
{
  AssertThrow(simulation_control != nullptr,
              ExcMessage(
                "simulation_control is null in SolidPhaseSolver constructor."));
}



template <int dim>
void
SolidPhaseSolver<dim>::setup_dofs()
{
  if (parameters.verbose_assembly)
    pcout << "setting up dofs...\n" << std::endl;

  TimerOutput::Scope t(computing_timer, "setup");

  dof_handler.distribute_dofs(fe);

  std::vector<unsigned int> block_component(fe.n_components(), 0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise(dof_handler, block_component);

  const auto dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, block_component);

  AssertDimension(dofs_per_block.size(), 2);

  const types::global_dof_index n_u = dofs_per_block[0];
  const types::global_dof_index n_a = dofs_per_block[1];

  pcout << "   Number of orders of freedom: " << dof_handler.n_dofs() << " ("
        << n_u << "+" << n_a << ")" << std::endl;

  locally_owned    = dof_handler.locally_owned_dofs();
  locally_relevant = DoFTools::extract_locally_relevant_dofs(dof_handler);

  owned_partitioning.clear();
  relevant_partitioning.clear();

  owned_partitioning.emplace_back(locally_owned.get_view(0, n_u));
  owned_partitioning.emplace_back(locally_owned.get_view(n_u, n_u + n_a));

  relevant_partitioning.emplace_back(locally_relevant.get_view(0, n_u));
  relevant_partitioning.emplace_back(locally_relevant.get_view(n_u, n_u + n_a));

  constraints.clear();
  constraints.reinit(locally_owned, locally_relevant);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  AssertThrow(simulation_control != nullptr,
              ExcMessage("simulation_control is null in setup_dofs()."));

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar alpha(dim);

  const ComponentMask vel_mask   = fe.component_mask(velocities);
  const ComponentMask alpha_mask = fe.component_mask(alpha);

  const double time = simulation_control->get_current_time();

  for (const auto &[id, type] : parameters.boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::noslip)
        {
          VectorTools::interpolate_boundary_values(dof_handler,
                                                   id,
                                                   Functions::ZeroFunction<dim>(
                                                     dim + 1),
                                                   constraints,
                                                   vel_mask);
        }
      else if (type == BoundaryConditions::BoundaryType::slip)
        {
          std::set<types::boundary_id> slip_boundaries;
          slip_boundaries.insert(id);

          VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                          0,
                                                          slip_boundaries,
                                                          constraints);
        }
      else if (type == BoundaryConditions::BoundaryType::function)
        {
          auto &velocity_functions =
            *parameters.boundary_conditions.velocity_functions.at(id);

          velocity_functions.u.set_time(time);
          velocity_functions.v.set_time(time);
          if constexpr (dim == 3)
            velocity_functions.w.set_time(time);

          auto &alpha_function =
            *parameters.boundary_conditions.alpha_functions.at(id);
          alpha_function.set_time(time);

          SolidVelocityFunctionDefined<dim> velocity_bc(
            &velocity_functions.u,
            &velocity_functions.v,
            (dim == 3 ? &velocity_functions.w : nullptr));

          SolidAlphaFunctionDefined<dim> alpha_bc(&alpha_function);

          VectorTools::interpolate_boundary_values(
            dof_handler, id, velocity_bc, constraints, vel_mask);

          VectorTools::interpolate_boundary_values(
            dof_handler, id, alpha_bc, constraints, alpha_mask);
        }
      else if (type == BoundaryConditions::BoundaryType::outlet)
        {
          // Natural boundary condition.

          // Outlet flux terms are assembled in assemble_boundary_face().
        }


      else if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            dof_handler,
            id,
            parameters.boundary_conditions.periodic_neighbor_id.at(id),
            parameters.boundary_conditions.periodic_direction.at(id),
            constraints);
        }
      else
        {
          AssertThrow(
            false,
            ExcMessage(
              "Unsupported solid boundary type in setup_dofs() for boundary ID " +
              std::to_string(static_cast<unsigned int>(id))));
        }
    }

  constraints.close();

  system_matrix.clear();

  TrilinosWrappers::BlockSparsityPattern sp(owned_partitioning,
                                            owned_partitioning,
                                            relevant_partitioning,
                                            mpi_communicator);

  DoFTools::make_sparsity_pattern(dof_handler,
                                  sp,
                                  constraints,
                                  false,
                                  Utilities::MPI::this_mpi_process(
                                    mpi_communicator));

  sp.compress();
  system_matrix.reinit(sp);

  solution.reinit(owned_partitioning, mpi_communicator);
  old_solution.reinit(owned_partitioning, mpi_communicator);
  older_solution.reinit(owned_partitioning, mpi_communicator);
  system_rhs.reinit(owned_partitioning, mpi_communicator);

  picard_solution.reinit(owned_partitioning, mpi_communicator);

  locally_relevant_solution.reinit(owned_partitioning,
                                   relevant_partitioning,
                                   mpi_communicator);

  locally_relevant_old_solution.reinit(owned_partitioning,
                                       relevant_partitioning,
                                       mpi_communicator);

  locally_relevant_older_solution.reinit(owned_partitioning,
                                         relevant_partitioning,
                                         mpi_communicator);

  locally_relevant_picard_solution.reinit(owned_partitioning,
                                          relevant_partitioning,
                                          mpi_communicator);

  picard_solution_ptr = nullptr;
}


template <int dim>
void
SolidPhaseSolver<dim>::set_fluid_velocity_field(
  const DoFHandler<dim>               &fluid_dh,
  const Mapping<dim>                  &fluid_mapping,
  const TrilinosWrappers::MPI::Vector &fluid_solution)
{
  fluid_dof_handler_ptr    = &fluid_dh;
  fluid_mapping_ptr        = &fluid_mapping;
  fluid_solution_ptr       = &fluid_solution;
  has_fluid_velocity_field = true;
}



template <int dim>
struct SolidCopyData
{
  FullMatrix<double>                   local_matrix;
  Vector<double>                       local_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  bool                                 cell_is_local = false;

  SolidCopyData(const unsigned int dofs_per_cell = 0)
    : local_matrix(dofs_per_cell, dofs_per_cell)
    , local_rhs(dofs_per_cell)
    , local_dof_indices(dofs_per_cell)
  {}

  void
  reset()
  {
    local_matrix = 0.0;
    local_rhs    = 0.0;
  }
};


template <int dim>
class SolidScratchData
{
public:
  SolidScratchData(const FiniteElement<dim>  &solid_fe_in,
                   const Quadrature<dim>     &cell_quadrature_in,
                   const Quadrature<dim - 1> &face_quadrature_in,
                   const Mapping<dim>        &fluid_mapping_in,
                   const FiniteElement<dim>  &fluid_fe_in)
    : solid_fe(&solid_fe_in)
    , cell_quadrature(&cell_quadrature_in)
    , face_quadrature(&face_quadrature_in)
    , fluid_mapping(&fluid_mapping_in)
    , fluid_fe(&fluid_fe_in)
    , fe_values(solid_fe_in,
                cell_quadrature_in,
                update_values | update_gradients | update_quadrature_points |
                  update_JxW_values)
    , fe_face_values(solid_fe_in,
                     face_quadrature_in,
                     update_values | update_normal_vectors |
                       update_quadrature_points | update_JxW_values)
    , fluid_fe_values(fluid_mapping_in,
                      fluid_fe_in,
                      cell_quadrature_in,
                      update_values | update_gradients)
    , n_q_points(cell_quadrature_in.size())
    , n_face_q_points(face_quadrature_in.size())
    , n_dofs(solid_fe_in.n_dofs_per_cell())
    , cell_diameter(0.0)

    // Cell data
    , JxW(n_q_points)
    , quadrature_points(n_q_points)
    , components(n_dofs)

    , phi_u(n_q_points, std::vector<Tensor<1, dim>>(n_dofs))
    , grad_phi_u(n_q_points, std::vector<Tensor<2, dim>>(n_dofs))

    , phi_a(n_q_points, std::vector<double>(n_dofs))
    , grad_phi_a(n_q_points, std::vector<Tensor<1, dim>>(n_dofs))

    , old_velocity_values(n_q_points)
    , old_velocity_gradients(n_q_points)
    , old_alpha_values(n_q_points)
    , old_alpha_gradients(n_q_points)

    , older_velocity_values(n_q_points)
    , older_alpha_values(n_q_points)

    , picard_velocity_values(n_q_points)
    , picard_velocity_gradients(n_q_points)
    , picard_velocity_divergences(n_q_points, 0.0)
    , picard_alpha_values(n_q_points)
    , picard_alpha_gradients(n_q_points)

    , fluid_velocity_values(n_q_points)
    , fluid_velocity_gradients(n_q_points)
    , fluid_velocity_divergences(n_q_points, 0.0)

    // Face data
    , face_JxW(n_face_q_points)
    , face_quadrature_points(n_face_q_points)
    , face_normals(n_face_q_points)

    , face_picard_velocity_values(n_face_q_points)
    , face_picard_alpha_values(n_face_q_points)

    , face_phi_u(n_face_q_points, std::vector<Tensor<1, dim>>(n_dofs))
    , face_phi_a(n_face_q_points, std::vector<double>(n_dofs))
  {
    for (unsigned int i = 0; i < n_dofs; ++i)
      components[i] = solid_fe_in.system_to_component_index(i).first;
  }

  SolidScratchData(const SolidScratchData<dim> &other)
    : solid_fe(other.solid_fe)
    , cell_quadrature(other.cell_quadrature)
    , face_quadrature(other.face_quadrature)
    , fluid_mapping(other.fluid_mapping)
    , fluid_fe(other.fluid_fe)
    , fe_values(*solid_fe,
                *cell_quadrature,
                update_values | update_gradients | update_quadrature_points |
                  update_JxW_values)
    , fe_face_values(*solid_fe,
                     *face_quadrature,
                     update_values | update_normal_vectors |
                       update_quadrature_points | update_JxW_values)
    , fluid_fe_values(*fluid_mapping,
                      *fluid_fe,
                      *cell_quadrature,
                      update_values | update_gradients)
    , n_q_points(other.n_q_points)
    , n_face_q_points(other.n_face_q_points)
    , n_dofs(other.n_dofs)
    , cell_diameter(0.0)

    // Cell data
    , JxW(n_q_points)
    , quadrature_points(n_q_points)
    , components(other.components)

    , phi_u(n_q_points, std::vector<Tensor<1, dim>>(n_dofs))
    , grad_phi_u(n_q_points, std::vector<Tensor<2, dim>>(n_dofs))

    , phi_a(n_q_points, std::vector<double>(n_dofs))
    , grad_phi_a(n_q_points, std::vector<Tensor<1, dim>>(n_dofs))

    , old_velocity_values(n_q_points)
    , old_velocity_gradients(n_q_points)
    , old_alpha_values(n_q_points)
    , old_alpha_gradients(n_q_points)

    , older_velocity_values(n_q_points)
    , older_alpha_values(n_q_points)

    , picard_velocity_values(n_q_points)
    , picard_velocity_gradients(n_q_points)
    , picard_velocity_divergences(n_q_points, 0.0)
    , picard_alpha_values(n_q_points)
    , picard_alpha_gradients(n_q_points)

    , fluid_velocity_values(n_q_points)
    , fluid_velocity_gradients(n_q_points)
    , fluid_velocity_divergences(n_q_points, 0.0)

    // Face data
    , face_JxW(n_face_q_points)
    , face_quadrature_points(n_face_q_points)
    , face_normals(n_face_q_points)

    , face_picard_velocity_values(n_face_q_points)
    , face_picard_alpha_values(n_face_q_points)

    , face_phi_u(n_face_q_points, std::vector<Tensor<1, dim>>(n_dofs))
    , face_phi_a(n_face_q_points, std::vector<double>(n_dofs))
  {}

  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const typename DoFHandler<dim>::active_cell_iterator &fluid_cell,
         const TrilinosWrappers::MPI::BlockVector             &old_solution,
         const TrilinosWrappers::MPI::BlockVector             &older_solution,
         const TrilinosWrappers::MPI::BlockVector             &picard_solution,
         const TrilinosWrappers::MPI::Vector                  &fluid_solution,
         const FEValuesExtractors::Vector                     &velocities,
         const FEValuesExtractors::Scalar                     &alpha,
         const FEValuesExtractors::Vector                     &fluid_velocities)
  {
    fe_values.reinit(cell);
    fluid_fe_values.reinit(fluid_cell);
    cell_diameter = cell->diameter();

    // Old solution: U^n
    fe_values[velocities].get_function_values(old_solution,
                                              old_velocity_values);
    fe_values[velocities].get_function_gradients(old_solution,
                                                 old_velocity_gradients);

    fe_values[alpha].get_function_values(old_solution, old_alpha_values);
    fe_values[alpha].get_function_gradients(old_solution, old_alpha_gradients);

    // Older solution: U^{n-1}
    fe_values[velocities].get_function_values(older_solution,
                                              older_velocity_values);
    fe_values[alpha].get_function_values(older_solution, older_alpha_values);

    // Picard solution: U^k
    fe_values[velocities].get_function_values(picard_solution,
                                              picard_velocity_values);
    fe_values[velocities].get_function_gradients(picard_solution,
                                                 picard_velocity_gradients);

    fe_values[alpha].get_function_values(picard_solution, picard_alpha_values);

    fe_values[alpha].get_function_gradients(picard_solution,
                                            picard_alpha_gradients);

    // Fluid velocity field
    fluid_fe_values[fluid_velocities].get_function_values(
      fluid_solution, fluid_velocity_values);

    fluid_fe_values[fluid_velocities].get_function_gradients(
      fluid_solution, fluid_velocity_gradients);

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        JxW[q]               = fe_values.JxW(q);
        quadrature_points[q] = fe_values.quadrature_point(q);

        picard_velocity_divergences[q] = trace(picard_velocity_gradients[q]);

        fluid_velocity_divergences[q] = trace(fluid_velocity_gradients[q]);

        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            phi_u[q][i]      = fe_values[velocities].value(i, q);
            grad_phi_u[q][i] = fe_values[velocities].gradient(i, q);

            phi_a[q][i]      = fe_values[alpha].value(i, q);
            grad_phi_a[q][i] = fe_values[alpha].gradient(i, q);
          }
      }
  }

  void
  reinit_face(const typename DoFHandler<dim>::active_cell_iterator &cell,
              const unsigned int                                    face,
              const TrilinosWrappers::MPI::BlockVector &picard_solution,
              const FEValuesExtractors::Vector         &velocities,
              const FEValuesExtractors::Scalar         &alpha)
  {
    fe_face_values.reinit(cell, face);

    fe_face_values[velocities].get_function_values(picard_solution,
                                                   face_picard_velocity_values);

    fe_face_values[alpha].get_function_values(picard_solution,
                                              face_picard_alpha_values);

    for (unsigned int q = 0; q < n_face_q_points; ++q)
      {
        face_JxW[q]               = fe_face_values.JxW(q);
        face_quadrature_points[q] = fe_face_values.quadrature_point(q);
        face_normals[q]           = fe_face_values.normal_vector(q);

        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            face_phi_u[q][i] = fe_face_values[velocities].value(i, q);

            face_phi_a[q][i] = fe_face_values[alpha].value(i, q);
          }
      }
  }

  const FiniteElement<dim>  *solid_fe;
  const Quadrature<dim>     *cell_quadrature;
  const Quadrature<dim - 1> *face_quadrature;
  const Mapping<dim>        *fluid_mapping;
  const FiniteElement<dim>  *fluid_fe;

  FEValues<dim>     fe_values;
  FEFaceValues<dim> fe_face_values;
  FEValues<dim>     fluid_fe_values;

  unsigned int n_q_points;
  unsigned int n_face_q_points;
  unsigned int n_dofs;
  double       cell_diameter;

  // Cell data
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;



  std::vector<unsigned int> components;

  std::vector<std::vector<Tensor<1, dim>>> phi_u;
  std::vector<std::vector<Tensor<2, dim>>> grad_phi_u;

  std::vector<std::vector<double>>         phi_a;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi_a;

  std::vector<Tensor<1, dim>> old_velocity_values;
  std::vector<Tensor<2, dim>> old_velocity_gradients;

  std::vector<double>         old_alpha_values;
  std::vector<Tensor<1, dim>> old_alpha_gradients;

  std::vector<Tensor<1, dim>> older_velocity_values;
  std::vector<double>         older_alpha_values;

  std::vector<Tensor<1, dim>> picard_velocity_values;
  std::vector<Tensor<2, dim>> picard_velocity_gradients;
  std::vector<double>         picard_velocity_divergences;
  std::vector<double>         picard_alpha_values;
  std::vector<Tensor<1, dim>> picard_alpha_gradients;

  std::vector<Tensor<1, dim>> fluid_velocity_values;
  std::vector<Tensor<2, dim>> fluid_velocity_gradients;
  std::vector<double>         fluid_velocity_divergences;

  // Face data
  std::vector<double>         face_JxW;
  std::vector<Point<dim>>     face_quadrature_points;
  std::vector<Tensor<1, dim>> face_normals;

  std::vector<Tensor<1, dim>> face_picard_velocity_values;
  std::vector<double>         face_picard_alpha_values;

  std::vector<std::vector<Tensor<1, dim>>> face_phi_u;
  std::vector<std::vector<double>>         face_phi_a;
};



template <int dim>
class SolidAssemblerBDF : public SolidTimeAssemblerBase<dim>
{
public:
  SolidAssemblerBDF(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const SolidPhaseParameters<dim>          &parameters)
    : simulation_control(simulation_control)
    , rho_s(parameters.rho_s)
    , dt(simulation_control->get_time_step())
    , assembly_method(Parameters::SimulationControl::TimeSteppingMethod::bdf1)
  {}

  void
  update_coefficients(const double /*unused*/) override
  {
    AssertThrow(simulation_control != nullptr,
                ExcMessage("simulation_control is null in SolidAssemblerBDF."));

    dt              = simulation_control->get_time_step();
    assembly_method = simulation_control->get_assembly_method();
  }

  void
  assemble(const SolidScratchData<dim> &scratch,
           SolidCopyData<dim>          &copy) const override
  {
    switch (assembly_method)
      {
        case Parameters::SimulationControl::TimeSteppingMethod::steady:
          return;

        case Parameters::SimulationControl::TimeSteppingMethod::steady_bdf:
        case Parameters::SimulationControl::TimeSteppingMethod::bdf1:
          assemble_bdf1(scratch, copy);
          return;

        case Parameters::SimulationControl::TimeSteppingMethod::bdf2:
          assemble_bdf2(scratch, copy);
          return;

        default:
          AssertThrow(
            false,
            ExcMessage(
              "Unsupported time stepping method in SolidAssemblerBDF."));
      }
  }

private:
  void
  assemble_bdf1(const SolidScratchData<dim> &scratch,
                SolidCopyData<dim>          &copy) const
  {
    const unsigned int n_q   = scratch.n_q_points;
    const unsigned int n_dof = scratch.n_dofs;

    for (unsigned int q = 0; q < n_q; ++q)
      {
        const double          JxW   = scratch.JxW[q];
        const Tensor<1, dim> &u_old = scratch.old_velocity_values[q];
        const double          a_old = scratch.old_alpha_values[q];

        /*
         * Picard coefficient for the unknown conservative momentum:
         *
         *   (alpha u)^{k+1} ≈ alpha^k u^{k+1}
         *
         * This appears only on the matrix side.
         */
        const double a_pic = scratch.picard_alpha_values[q];


        for (unsigned int i = 0; i < n_dof; ++i)
          {
            const unsigned int comp_i = scratch.components[i];

            // --------------------------------------------------
            // Alpha BDF1 time term:
            //
            // (1/dt) ∫ alpha^{k+1} q_i dΩ
            // =
            // (1/dt) ∫ alpha^n q_i dΩ
            // --------------------------------------------------
            if (comp_i == dim)
              {
                const double q_i = scratch.phi_a[q][i];

                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] != dim)
                      continue;

                    const double alpha_j = scratch.phi_a[q][j];

                    copy.local_matrix(i, j) += (q_i * alpha_j / dt) * JxW;
                  }

                copy.local_rhs(i) += (q_i * a_old / dt) * JxW;
              }

            // --------------------------------------------------
            // Momentum BDF1 time term:
            //
            // (rho_s/dt) ∫ alpha^k u^{k+1} · v_i dΩ
            // =
            // (rho_s/dt) ∫ alpha^n u^n · v_i dΩ
            //
            // Matrix uses alpha^k.
            // RHS uses alpha^n u^n.
            // --------------------------------------------------
            else if (comp_i < dim)
              {
                const Tensor<1, dim> &v_i = scratch.phi_u[q][i];

                const Tensor<1, dim> &u_pic = scratch.picard_velocity_values[q];

                // --------------------------------------------------
                // A_uu time block:
                //
                // (rho_s/dt) ∫ alpha^k u^{k+1} · v_i dΩ
                // --------------------------------------------------
                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] != comp_i)
                      continue;

                    const Tensor<1, dim> &u_j = scratch.phi_u[q][j];

                    copy.local_matrix(i, j) +=
                      (rho_s * a_pic / dt) * (v_i * u_j) * JxW;
                  }

                // --------------------------------------------------
                // A_u_alpha time coupling:
                //
                // (rho_s/dt) ∫ alpha^{k+1} u^k · v_i dΩ
                //
                // This is the new off-diagonal block A_{u alpha}.
                // --------------------------------------------------
                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] != dim)
                      continue;

                    const double alpha_j = scratch.phi_a[q][j];

                    copy.local_matrix(i, j) +=
                      (rho_s / dt) * alpha_j * (v_i * u_pic) * JxW;
                  }

                // --------------------------------------------------
                // RHS old conservative momentum:
                //
                // (rho_s/dt) ∫ alpha^n u^n · v_i dΩ
                // --------------------------------------------------
                copy.local_rhs(i) += (rho_s * a_old / dt) * (v_i * u_old) * JxW;

                // --------------------------------------------------
                // RHS correction from product linearization:
                //
                // (rho_s/dt) ∫ alpha^k u^k · v_i dΩ
                // --------------------------------------------------
                copy.local_rhs(i) += (rho_s * a_pic / dt) * (v_i * u_pic) * JxW;
              }
          }
      }
  }


  void
  assemble_bdf2(const SolidScratchData<dim> &scratch,
                SolidCopyData<dim>          &copy) const
  {
    const unsigned int n_q   = scratch.n_q_points;
    const unsigned int n_dof = scratch.n_dofs;

    for (unsigned int q = 0; q < n_q; ++q)
      {
        const double          JxW     = scratch.JxW[q];
        const double          a_old   = scratch.old_alpha_values[q];
        const double          a_older = scratch.older_alpha_values[q];
        const Tensor<1, dim> &u_old   = scratch.old_velocity_values[q];
        const Tensor<1, dim> &u_older = scratch.older_velocity_values[q];

        for (unsigned int i = 0; i < n_dof; ++i)
          {
            const unsigned int comp_i = scratch.components[i];

            // ----------------
            // Alpha equation
            // ----------------
            if (comp_i == dim)
              {
                const double w_i = scratch.phi_a[q][i];

                // implicit BDF2 mass: (3/2)/dt
                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] != dim)
                      continue;

                    const double a_j = scratch.phi_a[q][j];
                    copy.local_matrix(i, j) += (1.5 / dt) * (w_i * a_j) * JxW;
                  }

                // RHS: (2 a^n - 0.5 a^{n-1}) / dt
                const double rhs_alpha = (2.0 * a_old - 0.5 * a_older) / dt;
                copy.local_rhs(i) += w_i * rhs_alpha * JxW;
              }

            // ----------------
            // Momentum equation
            // ----------------
            else if (comp_i < dim)
              {
                const Tensor<1, dim> &phi_u_i = scratch.phi_u[q][i];

                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] != comp_i)
                      continue;

                    const Tensor<1, dim> &phi_u_j = scratch.phi_u[q][j];

                    // implicit mass: (3/2)(rho_s * a^{n+1})/dt
                    copy.local_matrix(i, j) +=
                      (1.5 * rho_s * scratch.picard_alpha_values[q] / dt) *
                      (phi_u_i * phi_u_j) * JxW;
                  }

                // RHS: (rho_s/dt)( 2 a^n u^n  - 0.5 a^{n-1} u^{n-1} )
                Tensor<1, dim> momentum_bdf_rhs =
                  (2.0 * a_old * u_old - 0.5 * a_older * u_older) / dt;

                copy.local_rhs(i) += (phi_u_i * momentum_bdf_rhs) * JxW;
              }
          }
      }
  }

  std::shared_ptr<SimulationControl>                simulation_control;
  const double                                      rho_s;
  double                                            dt;
  Parameters::SimulationControl::TimeSteppingMethod assembly_method;
};



template <int dim>
class SolidCoreAssembler
{
public:
  explicit SolidCoreAssembler(
    const SolidPhaseParameters<dim>          &parameters,
    const std::shared_ptr<SimulationControl> &simulation_control)
    : parameters(parameters)
    , simulation_control(simulation_control)
    , rho_s(parameters.rho_s)
    , beta(parameters.beta)
  {}

  void
  assemble_cell(const SolidScratchData<dim> &scratch,
                SolidCopyData<dim>          &copy) const
  {
    const unsigned int n_q   = scratch.n_q_points;
    const unsigned int n_dof = scratch.n_dofs;

    auto &A = copy.local_matrix;
    auto &F = copy.local_rhs;

    const double dt = simulation_control->get_time_step();

    // SUPG terms
    const bool   use_alpha_supg       = true;
    const double alpha_supg_factor    = 5.0;
    const double tiny_velocity        = 1e-12;
    const double tiny_cell_diameter   = 1e-12;
    const bool   use_momentum_supg    = true;
    const double momentum_supg_factor = 5.0;

    // grad-div stabilization terms
    const bool   use_solid_grad_div = true;
    const double grad_div_factor    = 0;

    for (unsigned int q = 0; q < n_q; ++q)
      {
        const double JxW = scratch.JxW[q];

        const Tensor<1, dim> &u_pic = scratch.picard_velocity_values[q];
        const double          a_pic = scratch.picard_alpha_values[q];

        const Tensor<1, dim> &u_f = scratch.fluid_velocity_values[q];

        Tensor<1, dim> u_const;
        u_const[0] = 1.0;
        if constexpr (dim >= 2)
          u_const[1] = 0.0;
        if constexpr (dim == 3)
          u_const[2] = 0.0;



        for (unsigned int i = 0; i < n_dof; ++i)
          {
            const unsigned int comp_i = scratch.components[i];

            // --------------------------------------------------
            // Continuity equation:
            // --------------------------------------------------
            if (comp_i == dim)
              {
                const Tensor<1, dim> &grad_q_i = scratch.grad_phi_a[q][i];

                const double h =
                  std::max(scratch.cell_diameter, tiny_cell_diameter);

                const double u_norm = std::max(u_pic.norm(), tiny_velocity);

                const double div_u_pic = scratch.picard_velocity_divergences[q];

                const double tau_alpha =
                  alpha_supg_factor /
                  std::sqrt((1.0 / dt) * (1.0 / dt) +
                            (2.0 * u_norm / h) * (2.0 * u_norm / h));

                const double u_dot_grad_q_i = u_pic * grad_q_i;

                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] != dim)
                      continue;

                    const double alpha_j = scratch.phi_a[q][j];

                    const Tensor<1, dim> &grad_alpha_j =
                      scratch.grad_phi_a[q][j];

                    double val = 0.0;

                    // --------------------------------------------------
                    // Standard conservative Galerkin volume term:
                    //
                    // - ∫ alpha^{k+1} u^k · grad(q_i) dΩ
                    // --------------------------------------------------
                    val += -alpha_j * (u_pic * grad_q_i);

                    // --------------------------------------------------
                    // Alpha SUPG stabilization:
                    // --------------------------------------------------
                    if (use_alpha_supg)
                      {
                        val += tau_alpha * u_dot_grad_q_i * (alpha_j / dt);

                        val +=
                          tau_alpha * u_dot_grad_q_i * (u_pic * grad_alpha_j);

                        val +=
                          tau_alpha * u_dot_grad_q_i * (alpha_j * div_u_pic);
                      }

                    A(i, j) += val * JxW;
                  }



                // SUPG RHS part:

                if (use_alpha_supg)
                  {
                    F(i) += tau_alpha * u_dot_grad_q_i *
                            (scratch.old_alpha_values[q] / dt) * JxW;
                  }


                // -------------------------------
                // A_alpha_u
                // -------------------------------
                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] >= dim)
                      continue;

                    const Tensor<1, dim> &u_trial_j = scratch.phi_u[q][j];

                    A(i, j) += -a_pic * (u_trial_j * grad_q_i) * JxW;
                  }

                // -------------------------------
                // RHS correction
                // -------------------------------
                F(i) += -a_pic * (u_pic * grad_q_i) * JxW;
              }



            // --------------------------------------------------
            // Momentum equation:
            // --------------------------------------------------
            else if (comp_i < dim)
              {
                const Tensor<1, dim> &v_i = scratch.phi_u[q][i];

                const Tensor<2, dim> &grad_v_i = scratch.grad_phi_u[q][i];

                const double h =
                  std::max(scratch.cell_diameter, tiny_cell_diameter);

                const double u_norm = std::max(u_pic.norm(), tiny_velocity);

                const double tau_u =
                  momentum_supg_factor /
                  std::sqrt((1.0 / dt) * (1.0 / dt) +
                            (2.0 * u_norm / h) * (2.0 * u_norm / h));

                const Tensor<1, dim> supg_test_u = grad_v_i * u_pic;

                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] >= dim)
                      continue;

                    const Tensor<1, dim> &u_trial_j = scratch.phi_u[q][j];

                    const Tensor<2, dim> &grad_u_j = scratch.grad_phi_u[q][j];


                    double val = 0.0;

                    val += -rho_s * a_pic * (u_trial_j * (grad_v_i * u_pic));

                    // + beta alpha^k v_i · u_trial_j
                    val += beta * a_pic * (v_i * u_trial_j);

                    if (use_momentum_supg)
                      {
                        Tensor<1, dim> residual_matrix_part;

                        // rho_s alpha^k u^{k+1} / dt
                        residual_matrix_part +=
                          (rho_s * a_pic / dt) * u_trial_j;

                        // rho_s alpha^k (u^k · grad) u^{k+1}
                        residual_matrix_part +=
                          rho_s * a_pic * (grad_u_j * u_pic);

                        // beta alpha^k u^{k+1}
                        residual_matrix_part += beta * a_pic * u_trial_j;

                        val += tau_u * (supg_test_u * residual_matrix_part);

                        val += tau_u * rho_s * a_pic *
                               (supg_test_u * (grad_u_j * u_pic));
                      }

                    A(i, j) += val * JxW;
                  }

                // Alpha columns: A_u_alpha
                // -------------------------------
                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] != dim)
                      continue;

                    const double alpha_j = scratch.phi_a[q][j];

                    // linearised convection

                    A(i, j) +=
                      -rho_s * alpha_j * (u_pic * (grad_v_i * u_pic)) * JxW;

                    // Linearized drag alpha-coupling:
                    // beta alpha^{k+1} (u_s^k - u_f) · v_i
                    A(i, j) += beta * alpha_j * (v_i * (u_pic - u_f)) * JxW;
                  }



                // --------------------------------------------------
                // RHS correction from convection product linearization:
                //
                // - rho_s ∫ alpha^k
                //          (u^k ⊗ u^k) : grad(v_i) dΩ
                // --------------------------------------------------
                F(i) += -rho_s * a_pic * (u_pic * (grad_v_i * u_pic)) * JxW;


                // RHS drag:
                F(i) += beta * a_pic * (v_i * u_pic) * JxW;


                if (use_momentum_supg)
                  {
                    Tensor<1, dim> residual_rhs_part;

                    // rho_s alpha^n u^n / dt
                    residual_rhs_part +=
                      (rho_s * scratch.old_alpha_values[q] / dt) *
                      scratch.old_velocity_values[q];

                    // beta alpha^k u_f
                    residual_rhs_part += beta * a_pic * u_f;

                    F(i) += tau_u * (supg_test_u * residual_rhs_part) * JxW;
                  }

                // --------------------------------------------------
                // grad-div stabilization
                // --------------------------------------------------
                if (use_solid_grad_div)
                  {
                    const double div_v_i = trace(grad_v_i);

                    const double h =
                      std::max(scratch.cell_diameter, tiny_cell_diameter);

                    const double u_norm = std::max(u_pic.norm(), tiny_velocity);

                    const double gamma_gd = grad_div_factor * u_norm * h;

                    const Tensor<1, dim> &grad_a_pic =
                      scratch.picard_alpha_gradients[q];

                    const double div_u_pic =
                      scratch.picard_velocity_divergences[q];

                    const double div_alpha_u_pic =
                      (u_pic * grad_a_pic) + a_pic * div_u_pic;

                    // --------------------------------------------------
                    // Velocity columns: contribution from
                    // div(alpha^k u^{k+1})
                    //
                    // div(alpha^k u_j)
                    // = grad(alpha^k) · u_j + alpha^k div(u_j)
                    // --------------------------------------------------
                    for (unsigned int j = 0; j < n_dof; ++j)
                      {
                        if (scratch.components[j] >= dim)
                          continue;

                        const Tensor<1, dim> &u_trial_j = scratch.phi_u[q][j];

                        const Tensor<2, dim> &grad_u_j =
                          scratch.grad_phi_u[q][j];

                        const double div_u_j = trace(grad_u_j);

                        const double residual_u_part =
                          (grad_a_pic * u_trial_j) + a_pic * div_u_j;

                        A(i, j) +=
                          rho_s * gamma_gd * residual_u_part * div_v_i * JxW;
                      }

                    // --------------------------------------------------
                    // Alpha columns: contribution from
                    // alpha^{k+1}/dt + div(alpha^{k+1} u^k)
                    //
                    // div(alpha_j u^k)
                    // = u^k · grad(alpha_j) + alpha_j div(u^k)
                    // --------------------------------------------------
                    for (unsigned int j = 0; j < n_dof; ++j)
                      {
                        if (scratch.components[j] != dim)
                          continue;

                        const double alpha_j = scratch.phi_a[q][j];

                        const Tensor<1, dim> &grad_alpha_j =
                          scratch.grad_phi_a[q][j];

                        const double residual_alpha_part =
                          alpha_j / dt + (u_pic * grad_alpha_j) +
                          alpha_j * div_u_pic;

                        A(i, j) += rho_s * gamma_gd * residual_alpha_part *
                                   div_v_i * JxW;
                      }

                    // --------------------------------------------------
                    // RHS contribution:
                    //
                    // alpha^n/dt + div(alpha^k u^k)
                    // --------------------------------------------------
                    F(i) +=
                      rho_s * gamma_gd *
                      ((scratch.old_alpha_values[q] / dt) + div_alpha_u_pic) *
                      div_v_i * JxW;
                  }
              }
          }
      }
  }


  void
  assemble_boundary_face(const SolidScratchData<dim> &scratch,
                         const types::boundary_id     boundary_id,
                         SolidCopyData<dim>          &copy) const
  {
    auto &A = copy.local_matrix;

    const unsigned int n_face_q = scratch.n_face_q_points;
    const unsigned int n_dof    = scratch.n_dofs;


    if (!is_outlet_boundary(boundary_id))
      return;

    for (unsigned int q = 0; q < n_face_q; ++q)
      {
        const double JxW = scratch.face_JxW[q];

        const Tensor<1, dim> &normal = scratch.face_normals[q];

        const Tensor<1, dim> &u_pic = scratch.face_picard_velocity_values[q];

        const double a_pic = scratch.face_picard_alpha_values[q];

        const double un = u_pic * normal;

        Tensor<1, dim> u_const;

        u_const[0] = 1.0;

        if constexpr (dim >= 2)
          u_const[1] = 0.0;

        if constexpr (dim == 3)
          u_const[2] = 0.0;

        const double un_const = u_const * normal;



        for (unsigned int i = 0; i < n_dof; ++i)
          {
            const unsigned int comp_i = scratch.components[i];

            // --------------------------------------------------
            // Continuity boundary flux:
            //
            // + ∫ alpha^{k+1} (u^k · n) q_i dΓ
            // --------------------------------------------------
            if (comp_i == dim)
              {
                if (un <= 0.0)
                  continue;

                const double q_i = scratch.face_phi_a[q][i];

                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] != dim)
                      continue;

                    const double alpha_j = scratch.face_phi_a[q][j];

                    A(i, j) += q_i * alpha_j * un * JxW;
                  }

                // --------------------------------------------------
                // A_alpha_u coupling from continuity outlet flux:
                //
                // + ∫ q_i alpha^k (u^{k+1} · n) dΓ
                // --------------------------------------------------
                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] >= dim)
                      continue;

                    const Tensor<1, dim> &u_trial_j = scratch.face_phi_u[q][j];

                    A(i, j) += q_i * a_pic * (u_trial_j * normal) * JxW;
                  }

                // --------------------------------------------------
                // RHS correction from continuity outlet flux linearization:
                //
                // + ∫ q_i alpha^k (u^k · n) dΓ
                // --------------------------------------------------
                copy.local_rhs(i) += q_i * a_pic * un * JxW;
              }

            // --------------------------------------------------
            // Momentum boundary flux:
            //
            // + rho_s ∫ alpha^k (u^k · n)
            //          v_i · u^{k+1} dΓ
            // --------------------------------------------------
            else if (comp_i < dim)
              {
                if (un <= 0.0)
                  continue;

                const Tensor<1, dim> &v_i = scratch.face_phi_u[q][i];

                // Velocity columns: A_uu
                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] >= dim)
                      continue;

                    const Tensor<1, dim> &u_trial_j = scratch.face_phi_u[q][j];

                    A(i, j) += rho_s * a_pic * un * (v_i * u_trial_j) * JxW;
                  }

                // Alpha columns: A_u_alpha
                for (unsigned int j = 0; j < n_dof; ++j)
                  {
                    if (scratch.components[j] != dim)
                      continue;

                    const double alpha_j = scratch.face_phi_a[q][j];

                    A(i, j) += rho_s * alpha_j * un * (v_i * u_pic) * JxW;
                  }

                // RHS correction
                copy.local_rhs(i) += rho_s * a_pic * un * (v_i * u_pic) * JxW;
              }
          }
      }
  }

private:
  // bool
  // is_strong_or_wall_boundary(const types::boundary_id boundary_id) const
  // {
  //   const auto it = parameters.boundary_conditions.type.find(boundary_id);

  //   if (it == parameters.boundary_conditions.type.end())
  //     {
  //       // Boundary not listed in the solid BC section:
  //       // treat as natural outlet/do-nothing.
  //       return false;
  //     }

  //   const auto type = it->second;

  //   if (type == BoundaryConditions::BoundaryType::function ||
  //       type == BoundaryConditions::BoundaryType::noslip ||
  //       type == BoundaryConditions::BoundaryType::slip ||
  //       type == BoundaryConditions::BoundaryType::periodic)
  //     {
  //       return true;
  //     }

  //   return false;
  // }

  BoundaryConditions::BoundaryType
  get_boundary_type(const types::boundary_id boundary_id) const
  {
    const auto it = parameters.boundary_conditions.type.find(boundary_id);

    AssertThrow(it != parameters.boundary_conditions.type.end(),
                ExcMessage(
                  "No solid boundary condition specified for boundary ID " +
                  std::to_string(static_cast<unsigned int>(boundary_id))));

    return it->second;
  }

  bool
  is_inlet_boundary(const types::boundary_id boundary_id) const
  {
    return get_boundary_type(boundary_id) ==
           BoundaryConditions::BoundaryType::function;
  }

  bool
  is_outlet_boundary(const types::boundary_id boundary_id) const
  {
    return get_boundary_type(boundary_id) ==
           BoundaryConditions::BoundaryType::outlet;
  }

  bool
  is_wall_boundary(const types::boundary_id boundary_id) const
  {
    const auto type = get_boundary_type(boundary_id);

    return type == BoundaryConditions::BoundaryType::noslip ||
           type == BoundaryConditions::BoundaryType::slip;
  }

  bool
  is_periodic_boundary(const types::boundary_id boundary_id) const
  {
    return get_boundary_type(boundary_id) ==
           BoundaryConditions::BoundaryType::periodic;
  }

  bool
  is_strong_or_wall_boundary(const types::boundary_id boundary_id) const
  {
    const auto type = get_boundary_type(boundary_id);

    return type == BoundaryConditions::BoundaryType::function ||
           type == BoundaryConditions::BoundaryType::noslip ||
           type == BoundaryConditions::BoundaryType::slip ||
           type == BoundaryConditions::BoundaryType::periodic;
  }

  const SolidPhaseParameters<dim>         &parameters;
  const std::shared_ptr<SimulationControl> simulation_control;
  const double                             rho_s;
  const double                             beta;
};


template <int dim>
void
SolidPhaseSolver<dim>::check_solid_mass_and_divergence() const
{
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar alpha(dim);

  const QGauss<dim>     cell_quadrature(order + 1);
  const QGauss<dim - 1> face_quadrature(order + 1);

  FEValues<dim> fe_values(fe,
                          cell_quadrature,
                          update_values | update_gradients | update_JxW_values);

  FEFaceValues<dim> fe_face_values(fe,
                                   face_quadrature,
                                   update_values | update_normal_vectors |
                                     update_JxW_values);

  const unsigned int n_q      = cell_quadrature.size();
  const unsigned int n_face_q = face_quadrature.size();

  std::vector<Tensor<1, dim>> u_values(n_q);
  std::vector<Tensor<2, dim>> grad_u_values(n_q);

  std::vector<double>         alpha_values(n_q);
  std::vector<double>         alpha_old_values(n_q);
  std::vector<Tensor<1, dim>> grad_alpha_values(n_q);

  std::vector<Tensor<1, dim>> face_u_values(n_face_q);
  std::vector<double>         face_alpha_values(n_face_q);

  double local_mass_new = 0.0;
  double local_mass_old = 0.0;
  double local_volume   = 0.0;

  double local_div_min = std::numeric_limits<double>::max();
  double local_div_max = -std::numeric_limits<double>::max();
  double local_div_l2  = 0.0;

  double local_alpha_res_l2  = 0.0;
  double local_alpha_res_min = std::numeric_limits<double>::max();
  double local_alpha_res_max = -std::numeric_limits<double>::max();

  double local_flux_total  = 0.0;
  double local_flux_wall   = 0.0;
  double local_flux_inlet  = 0.0;
  double local_flux_outlet = 0.0;

  double local_outlet_backflow_flux = 0.0;
  double local_outlet_backflow_area = 0.0;

  const double dt = simulation_control->get_time_step();

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      fe_values[velocities].get_function_values(locally_relevant_solution,
                                                u_values);
      fe_values[velocities].get_function_gradients(locally_relevant_solution,
                                                   grad_u_values);

      fe_values[alpha].get_function_values(locally_relevant_solution,
                                           alpha_values);
      fe_values[alpha].get_function_gradients(locally_relevant_solution,
                                              grad_alpha_values);

      fe_values[alpha].get_function_values(locally_relevant_old_solution,
                                           alpha_old_values);

      for (unsigned int q = 0; q < n_q; ++q)
        {
          const double JxW = fe_values.JxW(q);

          const double div_u = trace(grad_u_values[q]);

          const double alpha_res =
            (alpha_values[q] - alpha_old_values[q]) / dt +
            (u_values[q] * grad_alpha_values[q]) + alpha_values[q] * div_u;

          local_mass_new += alpha_values[q] * JxW;
          local_mass_old += alpha_old_values[q] * JxW;
          local_volume += JxW;

          local_div_min = std::min(local_div_min, div_u);
          local_div_max = std::max(local_div_max, div_u);
          local_div_l2 += div_u * div_u * JxW;

          local_alpha_res_min = std::min(local_alpha_res_min, alpha_res);
          local_alpha_res_max = std::max(local_alpha_res_max, alpha_res);
          local_alpha_res_l2 += alpha_res * alpha_res * JxW;
        }

      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
           ++face)
        {
          if (!cell->face(face)->at_boundary())
            continue;

          fe_face_values.reinit(cell, face);

          fe_face_values[velocities].get_function_values(
            locally_relevant_solution, face_u_values);

          fe_face_values[alpha].get_function_values(locally_relevant_solution,
                                                    face_alpha_values);

          const types::boundary_id boundary_id =
            cell->face(face)->boundary_id();

          const auto boundary_it =
            parameters.boundary_conditions.type.find(boundary_id);

          AssertThrow(
            boundary_it != parameters.boundary_conditions.type.end(),
            ExcMessage(
              "No solid boundary condition specified for mesh boundary ID " +
              std::to_string(static_cast<unsigned int>(boundary_id))));

          const auto boundary_type = boundary_it->second;

          for (unsigned int q = 0; q < n_face_q; ++q)
            {
              const double JxW = fe_face_values.JxW(q);

              const Tensor<1, dim> normal = fe_face_values.normal_vector(q);

              const double un = face_u_values[q] * normal;

              const double flux = face_alpha_values[q] * un * JxW;

              // Function boundaries are currently treated as inlets.
              if (boundary_type == BoundaryConditions::BoundaryType::function)
                {
                  local_flux_inlet += flux;
                  local_flux_total += flux;
                }
              else if (boundary_type ==
                       BoundaryConditions::BoundaryType::outlet)
                {
                  local_flux_outlet += flux;
                  local_flux_total += flux;

                  if (un < 0.0)
                    {
                      local_outlet_backflow_flux += flux;
                      local_outlet_backflow_area += JxW;
                    }
                }
              else if (boundary_type ==
                         BoundaryConditions::BoundaryType::noslip ||
                       boundary_type == BoundaryConditions::BoundaryType::slip)
                {
                  local_flux_wall += flux;
                  local_flux_total += flux;
                }
              else if (boundary_type ==
                       BoundaryConditions::BoundaryType::periodic)
                {
                  // Periodic boundaries are paired internally.
                  // Do not include them in the physical boundary-flux total.
                }
              else
                {
                  AssertThrow(
                    false,
                    ExcMessage(
                      "Unsupported solid boundary type in mass diagnostics for "
                      "boundary ID " +
                      std::to_string(static_cast<unsigned int>(boundary_id))));
                }
            }
        }
    }

  const double mass_new = Utilities::MPI::sum(local_mass_new, mpi_communicator);

  const double mass_old = Utilities::MPI::sum(local_mass_old, mpi_communicator);

  const double volume = Utilities::MPI::sum(local_volume, mpi_communicator);

  const double div_min = Utilities::MPI::min(local_div_min, mpi_communicator);

  const double div_max = Utilities::MPI::max(local_div_max, mpi_communicator);

  const double div_l2 =
    std::sqrt(Utilities::MPI::sum(local_div_l2, mpi_communicator));

  const double div_rms = div_l2 / std::sqrt(std::max(volume, 1e-30));

  const double alpha_res_min =
    Utilities::MPI::min(local_alpha_res_min, mpi_communicator);

  const double alpha_res_max =
    Utilities::MPI::max(local_alpha_res_max, mpi_communicator);

  const double alpha_res_l2 =
    std::sqrt(Utilities::MPI::sum(local_alpha_res_l2, mpi_communicator));

  const double alpha_res_rms =
    alpha_res_l2 / std::sqrt(std::max(volume, 1e-30));

  const double flux_total =
    Utilities::MPI::sum(local_flux_total, mpi_communicator);

  const double flux_wall =
    Utilities::MPI::sum(local_flux_wall, mpi_communicator);

  const double flux_inlet =
    Utilities::MPI::sum(local_flux_inlet, mpi_communicator);

  const double flux_outlet =
    Utilities::MPI::sum(local_flux_outlet, mpi_communicator);

  const double outlet_backflow_flux =
    Utilities::MPI::sum(local_outlet_backflow_flux, mpi_communicator);

  const double outlet_backflow_area =
    Utilities::MPI::sum(local_outlet_backflow_area, mpi_communicator);

  const double dMdt = (mass_new - mass_old) / dt;

  const double mass_balance_error = dMdt + flux_total;

  const double mass_balance_relative =
    std::abs(mass_balance_error) /
    (std::abs(dMdt) + std::abs(flux_total) + 1e-30);

  pcout << "solid diagnostics:\n"
        << "  alpha mass old      = " << mass_old << "\n"
        << "  alpha mass new      = " << mass_new << "\n"
        << "  dM/dt               = " << dMdt << "\n"
        << "  boundary flux total = " << flux_total << "\n"
        << "    wall flux         = " << flux_wall << "\n"
        << "    inlet flux        = " << flux_inlet << "\n"
        << "    outlet flux       = " << flux_outlet << "\n"
        << "  mass balance error  = " << mass_balance_error << "\n"
        << "  relative error      = " << mass_balance_relative << "\n"
        << "  div(u_s) min/max    = " << div_min << " , " << div_max << "\n"
        << "  div(u_s) L2/RMS     = " << div_l2 << " , " << div_rms << "\n"
        << "  alpha residual min/max = " << alpha_res_min << " , "
        << alpha_res_max << "\n"
        << "  alpha residual L2/RMS  = " << alpha_res_l2 << " , "
        << alpha_res_rms << "\n"
        << "  outlet backflow flux   = " << outlet_backflow_flux << "\n"
        << "  outlet backflow area   = " << outlet_backflow_area << std::endl;
}



template <int dim>
void
SolidPhaseSolver<dim>::assemble_local_system(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  SolidScratchData<dim>                                &scratch,
  SolidCopyData<dim>                                   &copy)
{
  copy.cell_is_local = cell->is_locally_owned();
  if (!copy.cell_is_local)
    return;

  copy.reset();

  AssertThrow(time_assembler != nullptr,
              ExcMessage("time_assembler is null in assemble_local_system."));
  AssertThrow(core_assembler != nullptr,
              ExcMessage("core_assembler is null in assemble_local_system."));

  AssertThrow(has_fluid_velocity_field,
              ExcMessage("Fluid velocity field must be set before assembly."));
  AssertThrow(picard_solution_ptr != nullptr,
              ExcMessage(
                "picard_solution_ptr is null in assemble_local_system."));
  AssertThrow(fluid_dof_handler_ptr != nullptr,
              ExcMessage("fluid_dof_handler_ptr is null."));
  AssertThrow(fluid_mapping_ptr != nullptr,
              ExcMessage("fluid_mapping_ptr is null."));
  AssertThrow(fluid_solution_ptr != nullptr,
              ExcMessage("fluid_solution_ptr is null."));

  const FEValuesExtractors::Vector solid_velocities(0);
  const FEValuesExtractors::Scalar solid_alpha(dim);
  const FEValuesExtractors::Vector fluid_velocities(0);

  typename DoFHandler<dim>::active_cell_iterator fluid_cell(
    &(fluid_dof_handler_ptr->get_triangulation()),
    cell->level(),
    cell->index(),
    fluid_dof_handler_ptr);

  scratch.reinit(cell,
                 fluid_cell,
                 locally_relevant_old_solution,
                 locally_relevant_older_solution,
                 *picard_solution_ptr,
                 *fluid_solution_ptr,
                 solid_velocities,
                 solid_alpha,
                 fluid_velocities);

  time_assembler->assemble(scratch, copy);
  core_assembler->assemble_cell(scratch, copy);



  for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
    {
      if (!cell->face(face)->at_boundary())
        continue;

      scratch.reinit_face(
        cell, face, *picard_solution_ptr, solid_velocities, solid_alpha);

      core_assembler->assemble_boundary_face(scratch,
                                             cell->face(face)->boundary_id(),
                                             copy);
    }



  cell->get_dof_indices(copy.local_dof_indices);
}


template <int dim>
void
SolidPhaseSolver<dim>::assemble_system()
{
  TimerOutput::Scope t(computing_timer, "assembly");

  AssertThrow(has_fluid_velocity_field,
              ExcMessage("Fluid velocity field must be set before assembly."));
  AssertThrow(picard_solution_ptr != nullptr,
              ExcMessage("picard_solution_ptr is null in assemble_system."));
  AssertThrow(fluid_dof_handler_ptr != nullptr,
              ExcMessage("fluid_dof_handler_ptr is null."));
  AssertThrow(fluid_mapping_ptr != nullptr,
              ExcMessage("fluid_mapping_ptr is null."));
  AssertThrow(fluid_solution_ptr != nullptr,
              ExcMessage("fluid_solution_ptr is null."));

  system_matrix = 0.0;
  system_rhs    = 0.0;

  locally_relevant_old_solution = old_solution;
  locally_relevant_old_solution.update_ghost_values();

  locally_relevant_older_solution = older_solution;
  locally_relevant_older_solution.update_ghost_values();

  const QGauss<dim>     cell_quadrature(order + 1);
  const QGauss<dim - 1> face_quadrature(order + 1);

  SolidScratchData<dim> scratch_data(fe,
                                     cell_quadrature,
                                     face_quadrature,
                                     *fluid_mapping_ptr,
                                     fluid_dof_handler_ptr->get_fe());

  AssertThrow(time_assembler != nullptr, ExcMessage("time_assembler is null."));
  AssertThrow(core_assembler != nullptr, ExcMessage("core_assembler is null."));

  WorkStream::run(dof_handler.begin_active(),
                  dof_handler.end(),
                  *this,
                  &SolidPhaseSolver<dim>::assemble_local_system,
                  &SolidPhaseSolver<dim>::copy_local_system_to_global,
                  scratch_data,
                  SolidCopyData<dim>(fe.n_dofs_per_cell()));

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}


template <int dim>
void
SolidPhaseSolver<dim>::setup_time_assembler()
{
  time_assembler =
    std::make_shared<SolidAssemblerBDF<dim>>(simulation_control, parameters);
}

template <int dim>
void
SolidPhaseSolver<dim>::setup_core_assembler()
{
  core_assembler =
    std::make_shared<SolidCoreAssembler<dim>>(parameters, simulation_control);
}

template <int dim>
void
SolidPhaseSolver<dim>::copy_local_system_to_global(
  const SolidCopyData<dim> &copy)
{
  if (!copy.cell_is_local)
    return;

  constraints.distribute_local_to_global(copy.local_matrix,
                                         copy.local_rhs,
                                         copy.local_dof_indices,
                                         system_matrix,
                                         system_rhs);
}



template <int dim>
void
SolidPhaseSolver<dim>::solve()
{
  TimerOutput::Scope t(computing_timer, "solve");


  if (parameters.solver_verbose)
    pcout << "Solving solid system... " << std::flush;

  TrilinosWrappers::MPI::BlockVector distributed_solution(owned_partitioning,
                                                          mpi_communicator);
  distributed_solution = solution;

  constraints.distribute(distributed_solution);
  distributed_solution.compress(VectorOperation::insert);

  PrimitiveVectorMemory<TrilinosWrappers::MPI::BlockVector> mem;

  const double rhs_norm = system_rhs.l2_norm();
  const double tol =
    std::max(parameters.solver_abs_tol, parameters.solver_rel_tol * rhs_norm);

  SolverControl solver_control(parameters.solver_max_it, tol);



  // if (rebuild_preconditioner || !linear_preconditioners_are_initialized)
  //   {
  //     TimerOutput::Scope t_prec(computing_timer,
  //                               "setup linear preconditioners");

  setup_linear_preconditioners();

  //     linear_preconditioners_are_initialized = true;
  //   }

  const bool use_amg = (parameters.amg_sweeps > 0);

  SolidBlockDiagonalPreconditioner preconditioner(velocity_amg,
                                                  velocity_ilu,
                                                  alpha_ilu,
                                                  use_amg);

  SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(
    solver_control,
    mem,
    SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData(
      parameters.solver_restart));

  try
    {
      solver.solve(system_matrix,
                   distributed_solution,
                   system_rhs,
                   preconditioner);
    }
  catch (SolverControl::NoConvergence &)
    {
      if (parameters.solver_verbose)
        pcout << "\nNoConvergence: retrying with relaxed settings...\n";

      SolverControl retry_control(parameters.solver_max_it, tol);

      SolverFGMRES<TrilinosWrappers::MPI::BlockVector> retry_solver(
        retry_control,
        mem,
        SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData(
          std::max(parameters.solver_restart, 200u)));

      retry_solver.solve(system_matrix,
                         distributed_solution,
                         system_rhs,
                         preconditioner);

      solver_control = retry_control;
    }

  constraints.distribute(distributed_solution);
  distributed_solution.compress(VectorOperation::insert);

  solution = distributed_solution;
  solution.compress(VectorOperation::insert);

  locally_relevant_solution = solution;
  locally_relevant_solution.update_ghost_values();

  if (parameters.solver_verbose)
    pcout << solver_control.last_step()
          << " iterations. residual=" << solver_control.last_value()
          << std::endl;
}


template <int dim>
const TrilinosWrappers::MPI::Vector &
SolidPhaseSolver<dim>::get_solid_volume_fraction() const
{
  return locally_relevant_solution.block(1);
}



template <int dim>
void
SolidPhaseSolver<dim>::make_output_dir() const
{
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::filesystem::create_directories(parameters.output_folder);
    }
  MPI_Barrier(mpi_communicator);
}



template <int dim>
void
SolidPhaseSolver<dim>::output_results(const double /*time*/)
{
  TimerOutput::Scope t(computing_timer, "output");

  if (timestep_number % parameters.output_every != 0)
    return;

  if (parameters.output_verbose)
    {
      pcout << "Writing output..." << std::endl;
    }
  std::vector<std::string> names(dim, "u");
  names.emplace_back("alpha");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim + 1, DataComponentInterpretation::component_is_scalar);
  for (unsigned int d = 0; d < dim; ++d)
    {
      interpretation[d] =
        DataComponentInterpretation::component_is_part_of_vector;
    }
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);


  data_out.add_data_vector(locally_relevant_solution,
                           names,
                           DataOut<dim>::type_dof_data,
                           interpretation);


  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    {
      subdomain(i) =
        static_cast<float>(triangulation.locally_owned_subdomain());
    }
  data_out.add_data_vector(subdomain,
                           "subdomain",
                           DataOut<dim>::type_cell_data);

  data_out.build_patches();

  data_out.write_vtu_with_pvtu_record(parameters.output_folder,
                                      parameters.output_prefix,
                                      timestep_number,
                                      mpi_communicator,
                                      parameters.digits);
}



template <int dim>
void
SolidPhaseSolver<dim>::setup()
{
  setup_dofs();

  setup_time_assembler();
  setup_core_assembler();

  make_output_dir();

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar alpha(dim);

  const ComponentMask vel_mask   = fe.component_mask(velocities);
  const ComponentMask alpha_mask = fe.component_mask(alpha);

  const double time = simulation_control->get_current_time();

  auto &velocity_ic_functions = *parameters.initial_conditions.velocity;
  velocity_ic_functions.u.set_time(time);
  velocity_ic_functions.v.set_time(time);
  if constexpr (dim == 3)
    velocity_ic_functions.w.set_time(time);

  auto &alpha_ic_function = *parameters.initial_conditions.alpha;
  alpha_ic_function.set_time(time);

  SolidVelocityFunctionDefined<dim> velocity_ic(&velocity_ic_functions.u,
                                                &velocity_ic_functions.v,
                                                (dim == 3 ?
                                                   &velocity_ic_functions.w :
                                                   nullptr));

  SolidAlphaFunctionDefined<dim> alpha_ic(&alpha_ic_function);

  old_solution = 0.0;

  VectorTools::interpolate(dof_handler, velocity_ic, old_solution, vel_mask);
  VectorTools::interpolate(dof_handler, alpha_ic, old_solution, alpha_mask);

  constraints.distribute(old_solution);
  old_solution.compress(VectorOperation::insert);

  older_solution = old_solution;
  older_solution.compress(VectorOperation::insert);

  locally_relevant_old_solution = old_solution;
  locally_relevant_old_solution.update_ghost_values();

  locally_relevant_older_solution = older_solution;
  locally_relevant_older_solution.update_ghost_values();

  solution = old_solution;
  solution.compress(VectorOperation::insert);

  locally_relevant_solution = solution;
  locally_relevant_solution.update_ghost_values();

  picard_solution = solution;
  picard_solution.compress(VectorOperation::insert);

  locally_relevant_picard_solution = picard_solution;
  locally_relevant_picard_solution.update_ghost_values();

  picard_solution_ptr = &locally_relevant_picard_solution;

  timestep_number = 0;
}

template <int dim>
bool
SolidPhaseSolver<dim>::advance_one_step()
{
  // if (finished())
  //   return false;

  ++timestep_number;

  const double dt   = simulation_control->get_time_step();
  const double time = simulation_control->get_current_time();

  time_assembler->update_coefficients(dt);

  const unsigned int max_picard = parameters.picard_max_iterations;
  const double       picard_tol = parameters.picard_tolerance;

  TrilinosWrappers::MPI::BlockVector prev_picard(owned_partitioning,
                                                 mpi_communicator);
  TrilinosWrappers::MPI::BlockVector delta(owned_partitioning,
                                           mpi_communicator);


  const auto print_alpha_bounds =
    [&](const TrilinosWrappers::MPI::Vector &alpha, const std::string &label) {
      double local_min = std::numeric_limits<double>::max();
      double local_max = -std::numeric_limits<double>::max();

      for (const auto i : alpha.locally_owned_elements())
        {
          local_min = std::min(local_min, alpha[i]);
          local_max = std::max(local_max, alpha[i]);
        }

      const double global_min =
        Utilities::MPI::min(local_min, mpi_communicator);

      const double global_max =
        Utilities::MPI::max(local_max, mpi_communicator);

      pcout << label << ": min = " << global_min << " max = " << global_max
            << std::endl;
    };

  const double omega = 0.3;

  for (unsigned int k = 0; k < max_picard; ++k)
    {
      picard_solution = solution;
      picard_solution.compress(VectorOperation::insert);

      locally_relevant_picard_solution = picard_solution;
      locally_relevant_picard_solution.update_ghost_values();

      picard_solution_ptr = &locally_relevant_picard_solution;

      prev_picard = solution;
      prev_picard.compress(VectorOperation::insert);

      if (k < 3 || k % 10 == 0)
        {
          print_alpha_bounds(locally_relevant_picard_solution.block(1),
                             "alpha before solve");
        }



      assemble_system();


      solve();

      if (k < 3 || k % 10 == 0)
        {
          print_alpha_bounds(locally_relevant_solution.block(1),
                             "alpha after solve");
        }



      // const auto &alpha = locally_relevant_solution.block(1);

      // double local_min = std::numeric_limits<double>::max();
      // double local_max = -std::numeric_limits<double>::max();

      // for (const auto i : alpha.locally_owned_elements())
      //   {
      //     local_min = std::min(local_min, alpha[i]);
      //     local_max = std::max(local_max, alpha[i]);
      //   }

      // const double global_min =
      //   Utilities::MPI::min(local_min, mpi_communicator);

      // const double global_max =
      //   Utilities::MPI::max(local_max, mpi_communicator);

      // pcout << "alpha bounds: min = " << global_min << " max = " <<
      // global_max
      //       << std::endl;



      TrilinosWrappers::MPI::BlockVector solved_solution(owned_partitioning,
                                                         mpi_communicator);
      solved_solution = solution;
      solved_solution.compress(VectorOperation::insert);

      solution = prev_picard;
      solution *= (1.0 - omega);
      solution.add(omega, solved_solution);

      constraints.distribute(solution);
      solution.compress(VectorOperation::insert);

      locally_relevant_solution = solution;
      locally_relevant_solution.update_ghost_values();

      if (k < 3 || k % 10 == 0)
        {
          print_alpha_bounds(locally_relevant_solution.block(1),
                             "alpha after relaxation");
        }

      delta = solution;
      delta -= prev_picard;

      const double picard_error =
        delta.l2_norm() / (solution.l2_norm() + 1e-30);

      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0 &&
          parameters.solver_verbose &&
          (k % 10 == 0 || picard_error < picard_tol || k == max_picard - 1))
        {
          // pcout << "alpha bounds: min = " << global_min
          //       << " max = " << global_max << std::endl;
          // print_alpha_bounds(locally_relevant_solution.block(1),
          //                    "alpha after relaxation");

          pcout << " Picard iteration " << k
                << " relative error = " << picard_error << std::endl;
        }

      if (picard_error < picard_tol)
        break;

      if (k == max_picard - 1 &&
          Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          pcout << " Warning: solid Picard did not converge. Final error = "
                << picard_error << std::endl;
        }
    }

  check_solid_mass_and_divergence();

  output_results(time);

  older_solution = old_solution;
  older_solution.compress(VectorOperation::insert);

  old_solution = solution;
  old_solution.compress(VectorOperation::insert);

  locally_relevant_older_solution = older_solution;
  locally_relevant_older_solution.update_ghost_values();

  locally_relevant_old_solution = old_solution;
  locally_relevant_old_solution.update_ghost_values();


  const double max_u_component =
    locally_relevant_solution.block(0).linfty_norm();

  const double max_solid_velocity =
    std::sqrt(static_cast<double>(dim)) * max_u_component;

  double local_h_min = std::numeric_limits<double>::max();

  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      local_h_min = std::min(local_h_min, cell->diameter());
    }

  const double h_min = Utilities::MPI::min(local_h_min, mpi_communicator);

  const double solid_cfl = max_solid_velocity * dt / h_min;

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0 &&
      parameters.solver_verbose)
    {
      pcout << "solid CFL = " << solid_cfl
            << "  max|u_s| ~= " << max_solid_velocity << "  h_min = " << h_min
            << "  dt = " << dt << std::endl;
    }

  if (parameters.print_timer_each_step)
    computing_timer.print_summary();

  return true;
}


template <int dim>
void
SolidPhaseSolver<dim>::finalize()
{
  if (parameters.print_timer_at_end)
    {
      computing_timer.print_summary();
    }
}

template <int dim>
bool
SolidPhaseSolver<dim>::finished() const
{
  return simulation_control->is_at_end();
}

template <int dim>
unsigned int
SolidPhaseSolver<dim>::get_step_number() const
{
  return timestep_number;
}

template <int dim>
double
SolidPhaseSolver<dim>::get_current_time() const
{
  return simulation_control->get_current_time();
}


template <int dim>
void
SolidPhaseSolver<dim>::run()
{
  setup();

  while (advance_one_step())
    {
    }

  finalize();
}



template class SolidPhaseSolver<2>;
template class SolidPhaseSolver<3>;
