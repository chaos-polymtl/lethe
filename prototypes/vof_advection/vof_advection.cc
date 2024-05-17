#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

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

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

template <int dim>
class AdvectionField : public TensorFunction<1, dim>
{
public:
  virtual Tensor<1, dim> value(const Point<dim> &p) const override;
};


template <int dim>
Tensor<1, dim> AdvectionField<dim>::value(const Point<dim> &p) const
{
  Tensor<1, dim> value;
  value[0] = 0.001;
  for (unsigned int i = 1; i < dim; ++i)
    value[i] = 0;

  return value;
}

template <int dim>
class InitialConditions : public Function<dim>
{
public:
  virtual double value(const Point<dim>  &p,
                       const unsigned int component = 0) const override;

private:
  static const Point<dim> center_point;
};

template <int dim>
double InitialConditions<dim>::value(const Point<dim>  &p,
                                  const unsigned int component) const
{
  (void)component;
  Assert(component == 0, ExcIndexRange(component, 0, 1));
  
  Point<dim> center = Point<dim>();
  Tensor<1,dim> dist = center - p;
  if (p[0]> 0.1 || p[0]<-0.1)
    return 0.0;
  return 1.0;
}

template <int dim>
class AdvectionProblem
{
public:
  AdvectionProblem();
  void run();

private:
  
  using VectorType = TrilinosWrappers::MPI::Vector;
  using MatrixType = TrilinosWrappers::SparseMatrix;
  
  void setup_system();
  
  struct AssemblyScratchData
  {
    // Constructors
    AssemblyScratchData(const FiniteElement<dim> &fe);
    AssemblyScratchData(const AssemblyScratchData &scratch_data);
    
    // Set up FEValues to reuse them because their initialization is expensive.
    FEValues<dim>     fe_values;
    
    // Previous phase values
    std::vector<double>    previous_phase_values;
    
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
  
  void assemble_system();
  void local_assemble_system(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    AssemblyScratchData                                  &scratch,
    AssemblyCopyData                                     &copy_data);    
  void copy_local_to_global(const AssemblyCopyData &copy_data);
  
  void set_initial_conditions();
  void solve();
  void refine_grid();
  void output_results() const;
  
  parallel::distributed::Triangulation<dim> triangulation;
  const MappingQ<dim>                       mapping;
  
  FE_Q<dim>                 fe;
  DoFHandler<dim>           dof_handler;
  AffineConstraints<double> constraints;
  MatrixType                system_matrix;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;
  
  VectorType solution;
  VectorType previous_solution;
  
  VectorType system_rhs;
  
  double dt;
  
  MPI_Comm           mpi_communicator;
  ConditionalOStream pcout;
};

template <int dim>
AdvectionProblem<dim>::AdvectionProblem()
  : triangulation(MPI_COMM_WORLD,
                  typename Triangulation<dim>::MeshSmoothing(
                    Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening))
  , mapping(2)
  , fe(1)
  , dof_handler(triangulation)
  , dt(0.01)
  , mpi_communicator(MPI_COMM_WORLD)
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
{}
  
template <int dim>
void AdvectionProblem<dim>::setup_system()
{
  std::cout << "boop 0" << std::endl;
  
  dof_handler.distribute_dofs(fe);
  
  std::cout << "boop 1" << std::endl;
  
  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
  
  std::cout << "boop 2" << std::endl;
  
  solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  previous_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
                     
  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler,
                                          constraints);
  constraints.close();

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs =*/false);
  
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

}  
  
template <int dim>
void AdvectionProblem<dim>::assemble_system()
{
  WorkStream::run(dof_handler.begin_active(),
                  dof_handler.end(),
                  *this,
                  &AdvectionProblem::local_assemble_system,
                  &AdvectionProblem::copy_local_to_global,
                  AssemblyScratchData(fe),
                  AssemblyCopyData());
}      

template <int dim>
AdvectionProblem<dim>::AssemblyScratchData::AssemblyScratchData(
  const FiniteElement<dim> &fe)
  : fe_values(fe,
              QGauss<dim>(fe.degree + 1),
              update_values | update_gradients | update_quadrature_points |
                update_JxW_values)
  , advection_directions(fe_values.get_quadrature().size())
{
  const unsigned int n_q_points =
    fe_values.get_quadrature().size();
    
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
  const unsigned int n_q_points =
    fe_values.get_quadrature().size();
    
  this->previous_phase_values = std::vector<double>(n_q_points);
}

template <int dim>
void AdvectionProblem<dim>::local_assemble_system(
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
  scratch_data.fe_values.get_function_values(previous_solution, scratch_data.previous_phase_values);
  
  scratch_data.advection_field.value_list(
    scratch_data.fe_values.get_quadrature_points(),
    scratch_data.advection_directions);
    
  const double cell_size = cell->diameter();

  
  double dt_inv = 1.0/this->dt;
  
  const auto &sd = scratch_data;
  for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
    {
      const Tensor<1, dim> velocity = sd.advection_directions[q_point];
      
      const double u_mag = std::max(velocity.norm(), 1e-12);
      
      const double tau =   1. /
                std::sqrt(Utilities::fixed_power<2>(dt_inv) +
                          Utilities::fixed_power<2>(2. * u_mag / cell_size));
                          
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        
        
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            // LHS
            
            // Time integration
            copy_data.cell_matrix(i, j) += dt_inv*sd.fe_values.shape_value(i, q_point)*sd.fe_values.shape_value(j, q_point)* sd.fe_values.JxW(q_point);
            // Advective term
            copy_data.cell_matrix(i, j) += sd.fe_values.shape_value(i, q_point)*sd.advection_directions[q_point]*sd.fe_values.shape_grad(j, q_point)*sd.fe_values.JxW(q_point);
            // Stabilization term
            
            copy_data.cell_matrix(i, j) += tau*sd.advection_directions[q_point]*sd.fe_values.shape_grad(i, q_point)*(sd.advection_directions[q_point]*sd.fe_values.shape_grad(j, q_point) + sd.fe_values.shape_value(j, q_point)*dt_inv);
            
          }
          //RHS
          
          // Stabilization term
          copy_data.cell_rhs(i) += tau*sd.advection_directions[q_point]*sd.fe_values.shape_grad(i, q_point)*dt_inv*scratch_data.previous_phase_values[q_point]*sd.fe_values.JxW(q_point);
          
          copy_data.cell_rhs(i) += dt_inv*sd.fe_values.shape_value(i, q_point)*scratch_data.previous_phase_values[q_point]*sd.fe_values.JxW(q_point);
      }
      
    }
    cell->get_dof_indices(copy_data.local_dof_indices);
    
}

template <int dim>
void
AdvectionProblem<dim>::copy_local_to_global(const AssemblyCopyData &copy_data)
{
  constraints.distribute_local_to_global(
    copy_data.cell_matrix,
    copy_data.cell_rhs,
    copy_data.local_dof_indices,
    system_matrix,
    system_rhs);
}

template <int dim>
void
AdvectionProblem<dim>::set_initial_conditions()
{
  VectorTools::interpolate(this->dof_handler,
                           InitialConditions<dim>(),
                           this->previous_solution);
  this->solution = this->previous_solution;                 
}

template <int dim>
void AdvectionProblem<dim>::solve()
{
  SolverControl               solver_control(200,
                               1e-10 * system_rhs.l2_norm());
  TrilinosWrappers::SolverGMRES solver(solver_control);
  
  TrilinosWrappers::PreconditionILU                 preconditioner;
  TrilinosWrappers::PreconditionILU::AdditionalData data_ilu;
  
  preconditioner.initialize(system_matrix, data_ilu);

  solver.solve(system_matrix,
              solution,
              system_rhs,
              preconditioner);

  VectorType residual(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  
  system_matrix.vmult(residual, solution);
  // residual -= system_rhs;
  std::cout << "   Iterations required for convergence: "
            << solver_control.last_step() << '\n'
            << "   Max norm of residual:                "
            << residual.linfty_norm() << '\n';

  constraints.distribute(solution);
}

template <int dim>
void AdvectionProblem<dim>::output_results() const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches(8);

  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::CompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);

  const std::string filename = "solution.vtu";
  std::ofstream     output(filename);
  data_out.write_vtu(output);
  std::cout << "Solution written to " << filename << std::endl;
}

template <int dim>
void AdvectionProblem<dim>::run()
{
  std::cout << "Bonjour from run" << std::endl;
  
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(3);
            
  std::cout << "Bonjour from after triangulation" << std::endl;
  // initial time step
  setup_system();
  
  std::cout << "Bonjour from after setup_system()" << std::endl;
  
  
  set_initial_conditions();
  
  std::cout << "Bonjour from after set_initial_conditions()" << std::endl;
  
  assemble_system();
  
  std::cout << "Bonjour from after assemble_system()" << std::endl;

  solve();
  
  std::cout << "Bonjour from after solve()" << std::endl;
  
  
  output_results();
  
  // 
  // previous_solution = solution;
  

}

int main(int argc, char *argv[])
{  
Utilities::MPI::MPI_InitFinalize       mpi_init(argc, argv, 1);
try
  {
    
    AdvectionProblem<2> advection_problem_2d;
    advection_problem_2d.run();
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


