// SPDX-FileCopyrightText: Copyright (c) 2019 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Prototype of advection solver using Discontinuous Galerkin method with
// moment-based numerical fluxes. This code is inspired by step-12 and step-59
// of the deal.II library, and implements a 2D advection problem with a given
// velocity field.

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// We want to impose a velocity field for the advection term.
// It is defined as an angular motion around the origin (0., 0.)
// In 2D: u(0) = -x(1)/|x|, u(1) = x(0)/|x|
template <int dim>
class VelocityField : public Function<dim>
{
public:
  VelocityField()
    : Function<dim>(dim)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &value) const override;
};

template <int dim>
void
VelocityField<dim>::vector_value(const Point<dim> &p,
                                 Vector<double>   &values) const
{
  AssertDimension(values.size(), dim);
  // 3D velocity field not yet implemented
  Assert(dim == 2, ExcNotImplemented());
  // 2D velocity field of Couette flow
  if constexpr (dim == 2)
    {
      double invert_norm = 1. / std::max(p.norm(), 1e-10);
      values(0)          = -p(1) * invert_norm; // u(0) = -x(1)/|x|
      values(1)          = p(0) * invert_norm;  // u(1) =  x(0)/|x|
    }
}

// Right-hand side function (zero in the initial case. To be modified later)
template <int dim>
class RHS : public Function<dim>
{
public:
  RHS()
    : Function<dim>(1)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    values(0) = 0;
    values(1) = 0;
  }
};

// Structure to hold FEValues and FEInterfaceValues objects for each cell
template <int dim>
struct ScratchData
{
  ScratchData(const Mapping<dim>        &mapping,
              const FiniteElement<dim>  &fe,
              const Quadrature<dim>     &quadrature,
              const Quadrature<dim - 1> &quadrature_face,
              const UpdateFlags          update_flags = update_values |
                                               update_gradients |
                                               update_quadrature_points |
                                               update_JxW_values,
              const UpdateFlags interface_update_flags =
                update_values | update_gradients | update_quadrature_points |
                update_JxW_values | update_normal_vectors)
    : fe_values(mapping, fe, quadrature, update_flags)
    , fe_interface_values(mapping, fe, quadrature_face, interface_update_flags)
  {}


  ScratchData(const ScratchData<dim> &scratch_data)
    : fe_values(scratch_data.fe_values.get_mapping(),
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                scratch_data.fe_values.get_update_flags())
    , fe_interface_values(scratch_data.fe_interface_values.get_mapping(),
                          scratch_data.fe_interface_values.get_fe(),
                          scratch_data.fe_interface_values.get_quadrature(),
                          scratch_data.fe_interface_values.get_update_flags())
  {}

  FEValues<dim>          fe_values;
  FEInterfaceValues<dim> fe_interface_values;
};

struct CopyDataFace
{
  FullMatrix<double>                   cell_matrix;
  std::vector<types::global_dof_index> joint_dof_indices;
};

struct CopyData
{
  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<CopyDataFace>            face_data;

  template <class Iterator>
  void
  reinit(const Iterator &cell, unsigned int dofs_per_cell)
  {
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);

    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }
};

// Dirichlet boundary condition function
// It imposes a value for our variable at our lower boundary (y=0) of 0.5
// between the origin and x=0.5 and 0 elsewhere.
template <int dim>
class DirichletBC : public Function<dim>
{
public:
  DirichletBC()
    : Function<dim>(1)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    for (unsigned int d = 0; d < dim; ++d)
      {
        if (p(0) < 0.5)
          values(d) = 1.;
        else
          values(d) = 0.;
      }
  }
};

// Main class for the DG advection solver
template <int dim>
class DGAdvectionMoments
{
public:
  DGAdvectionMoments(const unsigned int &fe_degree     = 1,
                     const unsigned int &n_refinements = 4)
    : fe_degree(fe_degree)
    , n_refinements(n_refinements)
    , triangulation()
    , dof_handler(triangulation)
    , fe(fe_degree)
    , quadrature(fe_degree + 1)
    , quadrature_face(fe_degree + 1)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    , computing_timer(MPI_COMM_WORLD,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times){};

  void
  run();

private:
  // Methods:
  void
  make_grid();
  void
  refine_grid();
  void
  setup_dofs();
  void
  assemble_system();
  void
  assemble_cell_terms(const FEValues<dim>      &fe_values,
                      const FullMatrix<double> &cell_matrix,
                      Vector<double>           &local_vector);
  void
  assemble_neumann_boundary_terms(const FEFaceValues<dim>  &fe_face_values,
                                  const FullMatrix<double> &local_matrix,
                                  Vector<double>           &local_vector);
  void
  assemble_dirichlet_boundary_terms(const FEFaceValues<dim>  &fe_face_values,
                                    const FullMatrix<double> &local_matrix,
                                    Vector<double>           &local_vector,
                                    const double             &h);
  void
  assemble_flux_terms(const FEFaceValuesBase<dim> &fe_face_values,
                      const FEFaceValuesBase<dim> &fe_neighbor_face_values,
                      FullMatrix<double>          &vi_ui_matrix,
                      FullMatrix<double>          &vi_ue_matrix,
                      FullMatrix<double>          &ve_ui_matrix,
                      FullMatrix<double>          &ve_ue_matrix,
                      const double                &h);

  void
  distribute_local_flux_to_global(
    const FullMatrix<double>                   &vi_ui_matrix,
    const FullMatrix<double>                   &vi_ue_matrix,
    const FullMatrix<double>                   &ve_ui_matrix,
    const FullMatrix<double>                   &ve_ue_matrix,
    const std::vector<types::global_dof_index> &local_dof_indices,
    const std::vector<types::global_dof_index> &local_neighbor_dof_indices);

  void
  solve();

  void
  output_results() const;

  // Attributes:
  const unsigned int fe_degree;
  const unsigned int n_refinements;
  Triangulation<dim> triangulation;

  DoFHandler<dim> dof_handler;
  FE_DGQ<dim>     fe;
  SparsityPattern sparsity_pattern;

  const QGauss<dim>     quadrature;
  const QGauss<dim - 1> quadrature_face;

  ConditionalOStream pcout;
  TimerOutput        computing_timer;

  VelocityField<dim>     velocity_field;
  SparseMatrix<double>   system_matrix;
  const RHS<dim>         system_rhs;
  const DirichletBC<dim> dirichlet_bc_function;
  Vector<double>         solution;
};

template <int dim>
void
DGAdvectionMoments<dim>::make_grid()
{
  TimerOutput::Scope t(computing_timer, "make grid");

  // Make sure 2D. 3D not yet implemented
  Assert(dim == 2, ExcNotImplemented());

  // Create a simple square mesh
  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(n_refinements);
}

template <int dim>
void
DGAdvectionMoments<dim>::refine_grid()
{
  TimerOutput::Scope t(computing_timer, "refine grid");
}

template <int dim>
void
DGAdvectionMoments<dim>::setup_dofs()
{
  TimerOutput::Scope t(computing_timer, "setup dofs");

  // Distribute degrees of freedom according to our FE
  dof_handler.distribute_dofs(fe);

  // Generate dynamic sparsity parttern from number of dofs
  DynamicSparsityPattern dsp(dof_handler.n_dofs());

  // Add flux components to matrix sparsity pattern
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);

  // Store sparsity pattern in attribute
  sparsity_pattern.copy_from(dsp);

  // Initialize system matrix, solution and rhs vectors with their proper sizes
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  pcout << "Number of active cells: " << triangulation.n_active_cells()
        << std::endl
        << "Total number of cells: " << triangulation.n_cells() << std::endl
        << "Number of degrees of freedom: " << dof_handler.n_dofs()
        << std::endl;
}

template <int dim>
void
DGAdvectionMoments<dim>::assemble_system()
{
  TimerOutput::Scope t(computing_timer, "assemble system");

  // Since we want to use the MeshWorker framework to avoid code repetition, we
  // need to instantiate an iterator for the active cells
  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  // Instantiate Dirichlet boundary condition function to be evaluater for each
  // quadrature point
  const DirichletBC<dim> dirichlet_bc_function;

  // I
  const auto cell_worker = [&](const Iterator   &cell,
                               ScratchData<dim> &scratch_data,
                               CopyData         &copy_data) {
    TimerOutput::Scope t(computing_timer, "assemble cell terms");

    // Get number of DoFs per cell
    const unsigned int dofs_per_cell =
      scratch_data.fe_values.get_fe().dofs_per_cell;

    // Reinitialize copy_data and fe_values for the current cell
    copy_data.reinit(cell, dofs_per_cell);
    scratch_data.fe_values.reinit(cell);

    // Get quadrature points from fe_values
    const std::vector<Point<dim>> &q_points =
      scratch_data.get_quadrature_points();

    // Get references to fe_values and JxW values
    const FEValues<dim>       &fe_values = scratch_data.fe_values;
    const std::vector<double> &JxW       = fe_values.get_JxW_values();

    // Assemble cell internal terms
    for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
      {
        // Get velocity at quadrature point q using velocity field function
        Vector<double> velocity(dim);
        velocity_field.vector_value(q_points[q], velocity);

        // Loop over dofs in cell
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            // Gather shape function value and gradient at quadrature point q
            // The index i corresponds to the test function
            const double         phi_i      = fe_values.shape_value(i, q);
            const Tensor<1, dim> grad_phi_i = fe_values.shape_grad(i, q);

            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                // The index j corresponds to the trial function
                const double         phi_j      = fe_values.shape_value(j, q);
                const Tensor<1, dim> grad_phi_j = fe_values.shape_grad(j, q);

                // Add entry to the local cell matrix
                // Weak form: (\nabla \phi_i, u \phi_j)
                copy_data.cell_matrix(i, j) +=
                  (-grad_phi_i * velocity * phi_j) * JxW[q];
              }
          }
      }
  };

  // Assemble boundary face flux terms
  const auto boundary_face_worker = [&](const Iterator   &cell,
                               ScratchData<dim> &scratch_data,
                               CopyData         &copy_data) {
    TimerOutput::Scope t(computing_timer, "assemble boundary face terms");

    // Reinitialize FEInterfaceValues for the current cell
    FEInterfaceValues<dim> &fe_interface_values =
      scratch_data.fe_interface_values;
    fe_interface_values.reinit(cell);

    // Get quadrature points from fe_interface_values
    const std::vector<Point<dim>> &q_points =
      fe_interface_values.get_quadrature_points();

    // Get velocity in direction normal to the face
    Vector<double> velocity(dim);
    velocity_field.vector_value(q_points, velocity);


    // Evaluate Dirichlet BC at quadrature point q
    Vector<double> g(dim);
    dirichlet_bc_function.vector_value(q_points, g);

    // Get number of DoFs per face
    const unsigned int dofs_per_face =
      fe_interface_values.get_fe().n_dofs_per_cell();

    // Get normal unit vectors
    const std::vector<Tensor<1, dim>> &normals =
      fe_interface_values.get_normal_vectors();

    // Get JxW values
    const std::vector<double> &JxW = fe_interface_values.get_JxW_values();

    // Loop over quadrature points
    for (unsigned int q = 0; q < q_points.size(); ++q)
      {
        const double normal_velocity = velocity(q) * normals(q);

        // If normal_velocity > 0, flow is leaving the cell. As such, use
        // interior value
        if (normal_velocity > 0)
          {
            for (unsigned int i = 0; i < dofs_per_face; ++i)
              {
                // Gather shape function value at quadrature point q}
                // The index i corresponds to the test function
                const double phi_i = fe_interface_values.shape_value(i, q);

                for (unsigned int j = 0; j < dofs_per_face; ++j)
                  {
                    // The index j corresponds to the trial function
                    const double phi_j = fe_interface_values.shape_value(j, q);

                    // Upwind flux: (phi_i, u * phi_j)_face
                    copy_data.cell_matrix(i, j) +=
                      (phi_i * phi_j * normal_velocity) * JxW[q];
                  }
              }
          }

        else
          {
            for (unsigned int i = 0; i < dofs_per_face; ++i)
              {
                // Gather shape function value at quadrature point q
                const double phi_i = fe_interface_values.shape_value(i, q);

                // Add entry to local cell rhs vector
                // (phi_i, g * u * n) where g is the Dirichlet BC
                // value
                copy_data.cell_rhs(i) +=
                  (-phi_i * g(q) * normal_velocity) * JxW[q];
              }
          }
      }
  };
}



int
main()
{
  try
    {
      DGAdvectionMoments<2> dg_advection_moments(1, 4);
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
