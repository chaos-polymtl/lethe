#include <core/manifolds.h>
#include <core/pvd_handler.h>
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

#include <fstream>
#include <iomanip>
#include <limits>


// using namespace dealii;



template <int dim>
void
SolidPhaseSolver<dim>::setup_preconditioner()
{
  this->preconditioner = std::make_shared<SolidPreconditioner<dim>>();


  const bool use_amg = (parameters.amg_sweeps > 0);

  if (use_amg)
    {
      setup_AMG();
    }
  else
    {
      setup_ILU();
    }
}

template <int dim>
void
SolidPhaseSolver<dim>::setup_ILU()
{
  this->preconditioner->mode = SolidPreconditioner<dim>::Mode::ilu;

  // velocity block ILU
  preconditioner->ilu_u = std::make_shared<TrilinosWrappers::PreconditionILU>();
  {
    TrilinosWrappers::PreconditionILU::AdditionalData ilu_data;
    ilu_data.overlap = parameters.ilu_overlap;
    preconditioner->ilu_u->initialize(system_matrix.block(0, 0), ilu_data);
  }

  // alpha block ILU
  preconditioner->ilu_a = std::make_shared<TrilinosWrappers::PreconditionILU>();
  {
    TrilinosWrappers::PreconditionILU::AdditionalData ilu_data;
    ilu_data.overlap = parameters.ilu_overlap;
    preconditioner->ilu_a->initialize(system_matrix.block(1, 1), ilu_data);
  }
}

template <int dim>
void
SolidPhaseSolver<dim>::setup_AMG()
{
  this->preconditioner->mode = SolidPreconditioner<dim>::Mode::amg;

  // AMG on velocity block (0,0)
  preconditioner->amg_u = std::make_shared<TrilinosWrappers::PreconditionAMG>();
  {
    TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
    amg_data.elliptic              = parameters.amg_elliptic;
    amg_data.higher_order_elements = (degree > 1);
    amg_data.smoother_type         = "ILU";
    amg_data.smoother_sweeps       = parameters.amg_sweeps;
    amg_data.aggregation_threshold = parameters.amg_agg_threshold;
    amg_data.output_details        = false;

    preconditioner->amg_u->initialize(system_matrix.block(0, 0), amg_data);
  }

  // ILU on alpha block (1,1)
  preconditioner->ilu_a = std::make_shared<TrilinosWrappers::PreconditionILU>();
  {
    TrilinosWrappers::PreconditionILU::AdditionalData ilu_data;
    ilu_data.overlap = parameters.ilu_overlap;
    preconditioner->ilu_a->initialize(system_matrix.block(1, 1), ilu_data);
  }
}


template <int dim>
class SolidInitialValues : public Function<dim>
{
public:
  SolidInitialValues(const double &alpha0, const double &u0)
    : Function<dim>(dim + 1)
    , alpha_0(alpha0)
    , u_0(u0)
  {}


  virtual double
  value(const Point<dim> & /*p*/, const unsigned int component) const override
  {
    if (component < dim)
      {
        // initial solid velocity: (u_0, 0, 0, ...)
        return (0.0);
      }
    else
      {
        // initial solid volume fraction α_s^0
        return alpha_0;
      }
  }


private:
  const double alpha_0;
  const double u_0;
};



template <int dim>
class SolidBoundaryValues : public Function<dim>
{
public:
  SolidBoundaryValues(const double alpha_inlet, const Tensor<1, dim> &u_inlet)
    : Function<dim>(dim + 1)
    , alpha_in(alpha_inlet)
    , u_in(u_inlet)
  {}

  double
  value(const Point<dim> &, const unsigned int component) const override
  {
    if (component < dim)
      {
        return u_in[component];
      }
    else
      {
        return alpha_in;
      }
  }

private:
  const double         alpha_in;
  const Tensor<1, dim> u_in;
};

template <int dim>
class SolidVelocityBoundaryFromFluid : public Function<dim>
{
public:
  SolidVelocityBoundaryFromFluid(
    const Mapping<dim>                  &fluid_mapping,
    const DoFHandler<dim>               &fluid_dof_handler,
    const TrilinosWrappers::MPI::Vector &fluid_solution)
    : Function<dim>(dim + 1)
    , fluid_field(fluid_dof_handler, fluid_solution, fluid_mapping)
    , n_fluid_components(fluid_dof_handler.get_fe().n_components())
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    values.reinit(dim + 1);
    values = 0.0;

    Vector<double> fluid_values(n_fluid_components);
    fluid_field.vector_value(p, fluid_values);

    for (unsigned int d = 0; d < dim; ++d)
      values[d] = fluid_values[d];

    values[dim] = 0.0;
  }

private:
  Functions::FEFieldFunction<dim, TrilinosWrappers::MPI::Vector> fluid_field;
  const unsigned int n_fluid_components;
};



template <int dim>
SolidPhaseSolver<dim>::SolidPhaseSolver(
  const SolidPhaseParameters                &p,
  parallel::distributed::Triangulation<dim> &tria,
  MPI_Comm                                   comm)
  : parameters(p)
  , mpi_communicator(comm)
  , degree(parameters.degree)
  , fe(FE_Q<dim>(degree), dim, FE_Q<dim>(degree), 1)
  , triangulation(tria)
  , dof_handler(triangulation)
  , rho_s(parameters.rho_s)
  , beta(parameters.beta)
  , time_step(parameters.time_step)
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
  sub.assign(dim, 1);
  sub[0] = parameters.nx;
  if (dim > 1)
    {
      sub[1] = parameters.ny;
    }
  if (dim > 2)
    {
      sub[2] = parameters.nz;
    }

  inlet_velocity    = Tensor<1, dim>();
  inlet_velocity[0] = parameters.u_inlet_x;
  if constexpr (dim >= 2)
    {
      inlet_velocity[1] = parameters.u_inlet_y;
    }
  if constexpr (dim >= 3)
    {
      inlet_velocity[2] = parameters.u_inlet_z;
    }
}


// template <int dim>
// void
// SolidPhaseSolver<dim>::make_grid()
// {
//   Point<dim> p1, p2;


//   for (unsigned int d = 0; d < dim; ++d)
//     {
//       p1[d] = 0.0;
//       p2[d] = 1.0;
//     }



//   // sub.assign(dim, 1);
//   // sub[0] = 10;           // x-direction
//   // if (dim > 1) sub[1] = 10;  // y-direction
//   // if (dim > 2) sub[2] = 10;   // z-direction



//   GridGenerator::subdivided_hyper_subdivided_hyper_rectangle(triangulation,
//   sub, p1, p2);


//   triangulation.refine_global(parameters.global_refinement);


//   const double tol = 1e-12;


//   const auto axis_from = [&](const std::string &d) -> unsigned int {
//     if (d == "x")
//       return 0;
//     if (d == "y")
//       return 1;
//     if (d == "z")
//       return 2;
//     AssertThrow(false, ExcMessage("direction must be x, y, or z"));
//     return 0;
//   };

//   const unsigned int axis1 = axis_from(parameters.direction1);
//   const unsigned int axis2 = axis_from(parameters.direction2);

//   pcout << "dim = " << dim << std::endl;

//   for (const auto &cell : triangulation.active_cell_iterators())
//     {
//       for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
//         {
//           if (cell->face(f)->at_boundary())
//             {
//               const Point<dim> fc = cell->face(f)->center();

//               cell->face(f)->set_boundary_id(0);

//               if (std::fabs(fc[axis1]) < tol)
//                 {
//                   cell->face(f)->set_boundary_id(1);
//                 }

//               else if (std::fabs(fc[axis1] - 1.0) < tol)
//                 {
//                   cell->face(f)->set_boundary_id(2);
//                 }

//               else if (std::fabs(fc[axis2]) < tol)
//                 {
//                   cell->face(f)->set_boundary_id(1);
//                 }

//               else if (std::fabs(fc[axis2] - 1.0) < tol)
//                 {
//                   cell->face(f)->set_boundary_id(2);
//                 }
//             }
//         }
//     }

//   //  if (std::fabs(fc[0]) < tol || (dim >= 2 && std::fabs(fc[1]) < tol))
//   //    cell->face(f)->set_boundary_id(1);      // inlet
//   //  else if (std::fabs(fc[0] - 1.0) < tol || (dim >= 2 && std::fabs(fc[1]
//   //  - 1.0) < tol))
//   //    cell->face(f)->set_boundary_id(2);      // outlet
//   //  else
//   //    cell->face(f)->set_boundary_id(0);      // walls



//   pcout << "active cells " << triangulation.n_global_active_cells()
//         << std::endl;
// }


template <int dim>
void
SolidPhaseSolver<dim>::setup_dofs()
{
  TimerOutput::Scope t(computing_timer, "setup");

  dof_handler.distribute_dofs(fe);

  std::vector<unsigned int> block_component(fe.n_components(), 0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise(dof_handler, block_component);

  const auto dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, block_component);

  const types::global_dof_index n_u = dofs_per_block[0];
  const types::global_dof_index n_a = dofs_per_block[1];

  pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << " ("
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


  SolidBoundaryValues<dim> bc(parameters.alpha_inlet, inlet_velocity);

  ComponentMask vel_mask(fe.n_components(), false);
  for (unsigned int c = 0; c < dim; ++c)
    {
      vel_mask.set(c, true);
    }

  ComponentMask alpha_mask(fe.n_components(), false);
  alpha_mask.set(dim, true);


  VectorTools::interpolate_boundary_values(
    dof_handler, 1, bc, constraints, vel_mask);
  VectorTools::interpolate_boundary_values(
    dof_handler, 1, bc, constraints, alpha_mask);


  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert(0);

  VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                  0,
                                                  no_normal_flux_boundaries,
                                                  constraints);

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

  old_solution.reinit(owned_partitioning, mpi_communicator);
  system_rhs.reinit(owned_partitioning, mpi_communicator);

  locally_relevant_solution.reinit(owned_partitioning,
                                   relevant_partitioning,
                                   mpi_communicator);
  locally_relevant_old_solution.reinit(owned_partitioning,
                                       relevant_partitioning,
                                       mpi_communicator);
}


template <int dim>
void
SolidPhaseSolver<dim>::update_constraints()
{
  AssertThrow(
    has_fluid_velocity_field,
    ExcMessage(
      "Fluid velocity field must be set before updating solid constraints."));
  AssertThrow(fluid_dof_handler_ptr != nullptr,
              ExcMessage("fluid_dof_handler_ptr is null."));
  AssertThrow(fluid_mapping_ptr != nullptr,
              ExcMessage("fluid_mapping_ptr is null."));
  AssertThrow(fluid_solution_ptr != nullptr,
              ExcMessage("fluid_solution_ptr is null."));

  constraints.clear();
  constraints.reinit(locally_owned, locally_relevant);

  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  ComponentMask vel_mask(fe.n_components(), false);
  for (unsigned int c = 0; c < dim; ++c)
    vel_mask.set(c, true);

  ComponentMask alpha_mask(fe.n_components(), false);
  alpha_mask.set(dim, true);

  SolidVelocityBoundaryFromFluid<dim> velocity_bc(*fluid_mapping_ptr,
                                                  *fluid_dof_handler_ptr,
                                                  *fluid_solution_ptr);

  SolidBoundaryValues<dim> alpha_bc(parameters.alpha_inlet, inlet_velocity);

  VectorTools::interpolate_boundary_values(
    dof_handler, 1, velocity_bc, constraints, vel_mask);

  VectorTools::interpolate_boundary_values(
    dof_handler, 1, alpha_bc, constraints, alpha_mask);

  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert(0);

  VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                  0,
                                                  no_normal_flux_boundaries,
                                                  constraints);

  constraints.close();
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
void
SolidPhaseSolver<dim>::assemble_system()
{
  if (parameters.verbose_assembly)
    {
      pcout << "Assembling...\n" << std::endl;
    }

  TimerOutput::Scope t(computing_timer, "assembly");

  system_matrix = 0.0;
  system_rhs    = 0.0;

  const QGauss<dim> quadrature(degree + 1);

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar alpha(dim);

  FEValues<dim> fe_values(fe,
                          quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q           = quadrature.size();

  AssertThrow(
    has_fluid_velocity_field,
    ExcMessage(
      "Fluid velocity field must be set before assembling solid system."));
  AssertThrow(fluid_dof_handler_ptr != nullptr,
              ExcMessage("fluid_dof_handler_ptr is null."));
  AssertThrow(fluid_mapping_ptr != nullptr,
              ExcMessage("fluid_mapping_ptr is null."));
  AssertThrow(fluid_solution_ptr != nullptr,
              ExcMessage("fluid_solution_ptr is null."));

  FEValues<dim> fluid_fe_values(*fluid_mapping_ptr,
                                fluid_dof_handler_ptr->get_fe(),
                                quadrature,
                                update_values);

  const FEValuesExtractors::Vector fluid_velocities(0);


  std::vector<Tensor<1, dim>> fluid_velocity_values(n_q);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


  const double beta  = parameters.beta;
  const double rho_s = parameters.rho_s;


  const double dt = time_step;



  // const Tensor<1, dim> u_f = inlet_velocity;



  // Shape function caches
  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
  std::vector<double>         phi_a(dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi_a(dofs_per_cell);

  // Old solution values
  std::vector<Tensor<1, dim>> u_old(n_q);
  std::vector<Tensor<2, dim>> grad_u_old(n_q);
  std::vector<double>         a_old(n_q);
  std::vector<Tensor<1, dim>> grad_a_old(n_q);


  locally_relevant_old_solution = old_solution;
  locally_relevant_old_solution.update_ghost_values();

  bool printed_u_f_q = false;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          cell_matrix = 0.0;
          cell_rhs    = 0.0;

          typename DoFHandler<dim>::active_cell_iterator fluid_cell(
            &(fluid_dof_handler_ptr->get_triangulation()),
            cell->level(),
            cell->index(),
            fluid_dof_handler_ptr);

          fluid_fe_values.reinit(fluid_cell);
          fluid_fe_values[fluid_velocities].get_function_values(
            *fluid_solution_ptr, fluid_velocity_values);

          fe_values[velocities].get_function_values(
            locally_relevant_old_solution, u_old);
          fe_values[velocities].get_function_gradients(
            locally_relevant_old_solution, grad_u_old);
          fe_values[alpha].get_function_values(locally_relevant_old_solution,
                                               a_old);
          fe_values[alpha].get_function_gradients(locally_relevant_old_solution,
                                                  grad_a_old);

          for (unsigned int q = 0; q < n_q; ++q)
            {
              const double         JxW     = fe_values.JxW(q);
              const Tensor<1, dim> u_k     = u_old[q];
              const double         a_k     = a_old[q];
              const double         div_u_k = trace(grad_u_old[q]);

              const Tensor<1, dim> u_f_q = fluid_velocity_values[q];

              if (!printed_u_f_q &&
                  Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
                {
                  pcout << "u_f_q = ";
                  for (unsigned int d = 0; d < dim; ++d)
                    pcout << u_f_q[d] << ' ';
                  pcout << std::endl;

                  printed_u_f_q = true;
                }

              // shape values at q
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_u[k]      = fe_values[velocities].value(k, q);
                  grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                  phi_a[k]      = fe_values[alpha].value(k, q);
                  grad_phi_a[k] = fe_values[alpha].gradient(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const unsigned int comp_i =
                    fe.system_to_component_index(i).first;

                  // ---------- alpha equation row ----------
                  if (comp_i == dim)
                    {
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                          const unsigned int comp_j =
                            fe.system_to_component_index(j).first;
                          if (comp_j != dim)
                            continue;

                          // (w, alpha^{n+1})/dt
                          cell_matrix(i, j) += (phi_a[i] * phi_a[j] / dt) * JxW;

                          // (w, u_k · ∇alpha^{n+1})
                          cell_matrix(i, j) +=
                            (phi_a[i] * (grad_phi_a[j] * u_k)) * JxW;

                          // (w, alpha^{n+1} div(u_k))

                          cell_matrix(i, j) +=
                            (phi_a[i] * phi_a[j] * div_u_k) * JxW;
                        }

                      // RHS: (w, alpha^n)/dt
                      cell_rhs(i) += (phi_a[i] * a_k / dt) * JxW;
                    }

                  // ---------- momentum equation rows
                  else if (comp_i < dim)
                    {
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                          const unsigned int comp_j =
                            fe.system_to_component_index(j).first;
                          if (comp_j != comp_i)
                            continue;

                          const double v_i =
                            fe_values[velocities].value(i, q)[comp_i];
                          const double u_j =
                            fe_values[velocities].value(j, q)[comp_i];

                          const Tensor<1, dim> grad_u_j =
                            fe_values[velocities].gradient(j, q)[comp_i];

                          const double adv = u_k * grad_u_j;

                          // mass
                          cell_matrix(i, j) +=
                            (rho_s * a_k / dt) * (v_i * u_j) * JxW;

                          // convection
                          cell_matrix(i, j) +=
                            (rho_s * a_k) * (v_i * adv) * JxW;

                          // drag
                          cell_matrix(i, j) += (beta * a_k) * (v_i * u_j) * JxW;
                        }

                      const double v_i =
                        fe_values[velocities].value(i, q)[comp_i];

                      // RHS
                      cell_rhs(i) +=
                        (rho_s * a_k / dt) * (v_i * u_k[comp_i]) * JxW;
                      cell_rhs(i) += (beta * a_k) * (v_i * u_f_q[comp_i]) * JxW;
                    }
                }
            }


          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 cell_rhs,
                                                 local_dof_indices,
                                                 system_matrix,
                                                 system_rhs);
        }
    }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  if (parameters.verbose_assembly)
    {
      pcout << "Done assembling.\n" << std::endl;
    }
}



template <int dim>
void
SolidPhaseSolver<dim>::solve()
{
  TimerOutput::Scope t(computing_timer, "solve");

  if (parameters.solver_verbose)
    {
      pcout << "Solving solid system... " << std::flush;
    }

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


  setup_preconditioner();

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
                   *preconditioner);
    }
  catch (SolverControl::NoConvergence &)
    {
      if (parameters.solver_verbose)
        {
          pcout << "\nNoConvergence: retrying with relaxed settings...\n";
        }
      SolverControl solver_control_refined(system_matrix.m(), tol);

      SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver_refined(
        solver_control_refined,
        mem,
        SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData(
          std::max(parameters.solver_restart, 200u)));

      solver_refined.solve(system_matrix,
                           distributed_solution,
                           system_rhs,
                           *preconditioner);
    }


  constraints.distribute(distributed_solution);
  distributed_solution.compress(VectorOperation::insert);

  solution = distributed_solution;
  solution.compress(VectorOperation::insert);


  locally_relevant_solution = solution;
  locally_relevant_solution.update_ghost_values();

  if (parameters.solver_verbose)
    {
      pcout << solver_control.last_step()
            << " iterations. residual=" << solver_control.last_value()
            << std::endl;
    }
}


template <int dim>
const TrilinosWrappers::MPI::Vector &
SolidPhaseSolver<dim>::get_solid_volume_fraction() const
{
  return locally_relevant_solution.block(1);
}

// template <int dim>
// const TrilinosWrappers::MPI::Vector &
// SolidPhaseSolver<dim>::get_solid_velocity() const
// {
//   return solid_velocity_solution;
// }


template <int dim>
void
SolidPhaseSolver<dim>::make_output_dir() const
{
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      mkdir(parameters.output_folder.c_str(), 0777);
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



// template <int dim>
// void
// SolidPhaseSolver<dim>::run()
// {
//   // make_grid();
//   setup_dofs();

//   make_output_dir();

//   // --- IC
//   SolidInitialValues<dim> ic(parameters.alpha0, 0.0);

//   old_solution = 0.0;
//   VectorTools::interpolate(dof_handler, ic, old_solution);
//   constraints.distribute(old_solution);
//   old_solution.compress(VectorOperation::insert);


//   locally_relevant_old_solution = old_solution;
//   locally_relevant_old_solution.update_ghost_values();

//   solution = old_solution;
//   solution.compress(VectorOperation::insert);

//   locally_relevant_solution = solution;
//   locally_relevant_solution.update_ghost_values();

//   timestep_number = 0;
//   output_results(0.0);

//   const double L = 1.0;
//   const double h = L / static_cast<double>(sub[0]);


//   double max_velocity = 0.0;
//   for (unsigned int d = 0; d < dim; ++d)
//     {
//       max_velocity = std::max(max_velocity, std::abs(inlet_velocity[d]));
//     }
//   for (timestep_number = 1; timestep_number <= parameters.n_steps;
//        ++timestep_number)
//     {
//       const double time = timestep_number * time_step;

//       assemble_system();
//       solve();



//       output_results(time);

//       old_solution = solution;
//       old_solution.compress(VectorOperation::insert);

//       locally_relevant_old_solution = old_solution;
//       locally_relevant_old_solution.update_ghost_values();

//       const double CFL = max_velocity * time_step / h;

//       pcout << "TimeStep " << timestep_number << " time = " << time
//             << " CFL = " << CFL << " ||rhs|| = " << system_rhs.l2_norm()
//             << "\n ";


//       if (parameters.print_timer_each_step)
//         computing_timer.print_summary();
//     }

//   if (parameters.print_timer_at_end)
//     {
//       computing_timer.print_summary();
//     }
// }


template <int dim>
void
SolidPhaseSolver<dim>::setup()
{
  // make_grid();
  setup_dofs();

  make_output_dir();

  // Initial condition
  SolidInitialValues<dim> ic(parameters.alpha0, 0.0);

  old_solution = 0.0;
  VectorTools::interpolate(dof_handler, ic, old_solution);
  constraints.distribute(old_solution);
  old_solution.compress(VectorOperation::insert);

  locally_relevant_old_solution = old_solution;
  locally_relevant_old_solution.update_ghost_values();

  solution = old_solution;
  solution.compress(VectorOperation::insert);

  locally_relevant_solution = solution;
  locally_relevant_solution.update_ghost_values();

  timestep_number = 0;
  output_results(0.0);

  const double L   = 1.0;
  cfl_length_scale = L / static_cast<double>(sub[0]);

  max_inlet_velocity = 0.0;
  for (unsigned int d = 0; d < dim; ++d)
    {
      max_inlet_velocity =
        std::max(max_inlet_velocity, std::abs(inlet_velocity[d]));
    }
}

template <int dim>
bool
SolidPhaseSolver<dim>::advance_one_step()
{
  if (timestep_number >= parameters.n_steps)
    return false;

  ++timestep_number;
  const double time = timestep_number * time_step;

  update_constraints();
  assemble_system();
  solve();

  output_results(time);

  old_solution = solution;
  old_solution.compress(VectorOperation::insert);

  locally_relevant_old_solution = old_solution;
  locally_relevant_old_solution.update_ghost_values();

  const double CFL = max_inlet_velocity * time_step / cfl_length_scale;

  pcout << "TimeStep " << timestep_number << " time = " << time
        << " CFL = " << CFL << " ||rhs|| = " << system_rhs.l2_norm() << "\n ";

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
  return timestep_number >= parameters.n_steps;
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
  return timestep_number * time_step;
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
