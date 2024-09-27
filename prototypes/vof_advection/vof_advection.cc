#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>

#include <deal.II/hp/q_collection.h>

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

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

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
#include <deal.II/base/data_out_base.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/non_matching/fe_immersed_values.h>
#include <deal.II/non_matching/fe_values.h>
#include <deal.II/non_matching/mesh_classifier.h>

#include <fstream>
#include <iostream>

using namespace dealii;

template <int dim>
inline double compute_point_2_interface_min_distance(const Point<dim> &point_0,const Point<dim> &point_1,const Point<dim> &y)
{
  const Tensor<1,dim> d = point_1 - point_0;
  
  const double t_bar = d*(y-point_0)/(d.norm()*d.norm());
  
  double D;
  
  if (t_bar <= 0.0)
  {
    const Tensor<1,dim> y_minus_p0 = y-point_0;
    D = y_minus_p0.norm();
  }
  else if (t_bar >= 1.0)
  {
    const Tensor<1,dim> y_minus_p1 = y-point_1;
    D = y_minus_p1.norm();
  }
  else
  {
    const Tensor<1,dim> projection = y-(point_0+t_bar*d);
    D = projection.norm();
  }
  
  return D;
}

template <int dim>
inline void get_dof_opposite_faces(unsigned int local_face_id, std::vector<unsigned int> &local_opposite_faces)
{
  unsigned int local_face_id_2d = local_face_id%4;
  
  local_opposite_faces[0] = (local_face_id_2d + 1)%2;
  local_opposite_faces[1] = 3-local_face_id_2d/2;
  
  if constexpr (dim==3)
    local_opposite_faces[2] = 5 - local_face_id/4;
}

template <int dim>
inline void get_face_transformation_jacobians(const DerivativeForm<1, dim, dim> &cell_transformation_jac, const unsigned int local_face_id, DerivativeForm<1, dim-1, dim> &face_transformation_jac)
{
  for (unsigned int i = 0; i < dim; ++i)
  {
    unsigned int k = 0;
    for (unsigned int j = 0; j < dim; ++j)
    {
      if (local_face_id/2 == j)
        continue;
      face_transformation_jac[i][k] = cell_transformation_jac[i][j];
      k += 1;
    }
  }
}

template <int dim>
inline void compute_residual(const Tensor<1,dim> &x_n_to_x_J_real, const Tensor<1, dim> &distance_gradient, const DerivativeForm<1, dim-1, dim> transformation_jac, Tensor<1, dim-1> &residual_ref)
{
  Tensor<1, dim> residual_real = distance_gradient - (1.0/x_n_to_x_J_real.norm())*x_n_to_x_J_real;
  
  DerivativeForm<1, dim, dim-1> transformation_jac_transpose = transformation_jac.transpose();
  
  for (unsigned int i = 0; i < dim-1; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      residual_ref[i] += transformation_jac_transpose[i][j]*residual_real[j];
}

template <int dim>
inline std::vector<Point<dim>> compute_numerical_jacobian_stencil(const Point<dim> x_ref, const unsigned int local_face_id, const double perturbation)
{
  
  std::vector<Point<dim>> stencil( 2*dim - 1);
  
  for (unsigned int i = 0; i < 2*dim - 1; ++i)
  {
    stencil[i] = x_ref;
    // unsigned int k = 0;
    // for (unsigned int j = 0; j < dim; ++j)
    // {
    //   if (local_face_id/2 == j)
    //     {
    //       stencil[i][j] = double(local_face_id%2);
    //       continue;
    //     }
    //   stencil[i][j] = x_ref[j];
    //   k += 1;
    // }
  }
  
  for (unsigned int i = 1; i < dim; ++i)
  {
    for (unsigned int j = 0; j < dim; ++j)
    {
      if (local_face_id/2 == j)
        continue;
      stencil[2*i-1][j] -= perturbation;
      stencil[2*i][j] += perturbation;
    }
  }
  
  return stencil;
}

template <int dim, typename VectorType>
inline Tensor<1,dim> transform_ref_face_correction_to_ref_cell(const VectorType &correction_ref_face, const unsigned int local_face_id)
{
  Tensor<1,dim> correction_ref_cell;
  
  unsigned int j = 0;
  for (unsigned int i = 0; i < dim; ++i)
  {
    if (local_face_id/2 == i)
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
inline Point<dim> transform_ref_face_point_to_ref_cell(const Point<dim-1> &x_ref_face, const unsigned int local_face_id)
{
  Point<dim> x_ref_cell;
  
  unsigned int j = 0;
  for (unsigned int i = 0; i < dim; ++i)
  {
    if (local_face_id/2 == i)
      {
        x_ref_cell[i] = double(local_face_id%2);
        continue;
      }
    x_ref_cell[i] = x_ref_face[j];
    j += 1;
  }
  
  return x_ref_cell;
}

template <int dim>
inline void compute_numerical_jacobians(const std::vector<Point<dim>> &stencil_real, const Point<dim> &x_J_real, const std::vector<Tensor<1, dim>> &distance_gradients, const std::vector<DerivativeForm<1, dim-1, dim>> &transformation_jacobians, const double perturbation, LAPACKFullMatrix<double> &jacobian_matrix)
{  
  for (unsigned int i = 0; i < dim-1; ++i)
  {
    const Tensor<1,dim> x_n_to_x_J_real_m1 = x_J_real - stencil_real[2*i+1];
  
    Tensor<1, dim-1> residual_ref_m1;
    compute_residual(x_n_to_x_J_real_m1, distance_gradients[2*i+1], transformation_jacobians[2*i+1], residual_ref_m1);
    
  
    const Tensor<1,dim> x_n_to_x_J_real_p1 = x_J_real - stencil_real[2*i+2]; 
  
    Tensor<1, dim-1> residual_ref_p1;
    compute_residual(x_n_to_x_J_real_p1, distance_gradients[2*i+2], transformation_jacobians[2*i+2], residual_ref_p1);
    
  
    for (unsigned int j = 0; j < dim-1; ++j)
    {  
      jacobian_matrix.set(j,i,(residual_ref_p1[j]-residual_ref_m1[j])/(2.0*perturbation));
    }
  }
}

template <int dim>
inline double compute_distance(const Tensor<1,dim> &x_n_to_x_J_real, const double distance)
{
  return distance + x_n_to_x_J_real.norm();
}

template <typename T>
int
sgn(T val)
{
  return (static_cast<T>(0) < val) - (val < static_cast<T>(0));
}

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
  value[0] = 1;
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
  
  double sign_dist = 1.0;
  if (dist.norm() > 0.25)
  {
    sign_dist = -1.0;
  }
  
  return 0.5+0.5*std::tanh((0.25-dist.norm())/0.02);
  
  // if (p[0]> 0.25 || p[0]<-0.25)
  //   return 0.0;
  // return 0.5 + 0.5*std::cos(4*p[0]*M_PI);
}

template <int dim>
class Visualization : public dealii::DataOutInterface<0, dim>
{
public:
  Visualization();
  
  void  
  build_patches(DoFHandler<dim> &dof_handler, NonMatching::FEValues<dim> &non_matching_fe_values);
  
private:
  /**
   * Implementation of the corresponding function of the base class.
   */
  virtual const std::vector<DataOutBase::Patch<0, dim>> &
  get_patches() const override;
  
  virtual std::vector<std::string>
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
Visualization<dim>::build_patches(DoFHandler<dim> &dof_handler, NonMatching::FEValues<dim> &non_matching_fe_values)
{

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
  
    
    if (cell->is_locally_owned())
    {
  
    
    
    non_matching_fe_values.reinit(cell);
  
    const std::optional<NonMatching::FEImmersedSurfaceValues<dim>>
      &surface_fe_values = non_matching_fe_values.get_surface_fe_values();
      
    
    
    if (surface_fe_values)
    {
      
      for (const unsigned int q :
           surface_fe_values->quadrature_point_indices())
        {
          const Point<dim> &point = surface_fe_values->quadrature_point(q);
          // pcout << "x = " << point[0] << "\t y = " << point[1] << std::endl;
  
          DataOutBase::Patch<0, dim> temp;
          temp.vertices[0] = point;
          patches.push_back(temp);
        }
      }
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
  
  void compute_level_set_from_phase_fraction();
  void compute_phase_fraction_from_level_set();
  void compute_sign_distance();
  
  
  void refine_grid();
  void output_results(const int time_iteration, const double time) const;
  
  parallel::distributed::Triangulation<dim> triangulation;
  const MappingQ<dim>                       mapping;
  
  GridTools::Cache<dim> grid_tools_cache;
  
  const FE_Q<dim>           fe;
  DoFHandler<dim>           dof_handler;
  hp::FECollection<dim>     fe_collection;
  
  const FE_Q<dim> fe_level_set;
  DoFHandler<dim> level_set_dof_handler;
  VectorType      level_set;
      
  AffineConstraints<double> constraints;
  MatrixType                system_matrix;
  
  NonMatching::MeshClassifier<dim> mesh_classifier;

  IndexSet locally_owned_dofs;
  types::global_dof_index n_locally_owned_dofs;
  
  IndexSet locally_relevant_dofs;
  
  VectorType solution;
  VectorType previous_solution;
  VectorType location;
  VectorType signed_distance;
  VectorType distance;
  
  std::map<types::global_cell_index,std::vector<Point<dim>>> intersection_point;
  
  std::map<types::global_cell_index,Point<dim>> intersection_cell;
  
  std::set<types::global_cell_index> intersection_halo_cell;
  std::set<types::global_dof_index> dofs_location_status;
  
    
  VectorType system_rhs;
  
  double dt;
  
  MPI_Comm           mpi_communicator;
  ConditionalOStream pcout;
};

enum ActiveFEIndex
{
  lagrange = 0,
  nothing  = 1
};
    
template <int dim>
AdvectionProblem<dim>::AdvectionProblem()
  : triangulation(MPI_COMM_WORLD,
                  typename Triangulation<dim>::MeshSmoothing(
                    Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening))
  , mapping(1)
  , grid_tools_cache(triangulation, mapping)
  , fe(1)
  , dof_handler(triangulation)
  , fe_level_set(1)
  , level_set_dof_handler(triangulation)
  , mesh_classifier(dof_handler, level_set)
  , dt(0.001)
  , mpi_communicator(MPI_COMM_WORLD)
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
{}
  
template <int dim>
void AdvectionProblem<dim>::setup_system()
{
  fe_collection.push_back(FE_Q<dim>(1));
      // fe_collection.push_back(FE_Nothing<dim>());
  
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // const NonMatching::LocationToLevelSet cell_location =
      //   mesh_classifier.location_to_level_set(cell);
      cell->set_active_fe_index(ActiveFEIndex::lagrange);
    }
        
  pcout << "boop 0" << std::endl;
  dof_handler.distribute_dofs(fe_collection);
  
  pcout << "boop 1" << std::endl;
  
  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  n_locally_owned_dofs    = dof_handler.n_locally_owned_dofs();
  
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
  
  pcout << "boop 2" << std::endl;
  
  solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  previous_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  level_set.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  location.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  signed_distance.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  distance.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  
  
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
  system_matrix = 0;
  system_rhs = 0;
  
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
    
  double cell_size;
  
  if (dim == 2)
    {
      cell_size = std::sqrt(4. * cell->measure() / M_PI);
    }
  else if (dim == 3)
    {
      cell_size = std::pow(6 * cell->measure() / M_PI, 1. / 3.);
    }
  
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
            
            copy_data.cell_matrix(i, j) += tau*sd.advection_directions[q_point]*sd.fe_values.shape_grad(i, q_point)*(sd.advection_directions[q_point]*sd.fe_values.shape_grad(j, q_point) + sd.fe_values.shape_value(j, q_point)*dt_inv)*sd.fe_values.JxW(q_point);
            
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
  SolverControl               solver_control(1000,
                               1e-10 * system_rhs.l2_norm() );
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
  pcout << "   Iterations required for convergence: "
            << solver_control.last_step() << '\n'
            << "   Max norm of residual:                "
            << residual.linfty_norm() << '\n';

  constraints.distribute(solution);
}

template <int dim>
void 
AdvectionProblem<dim>::compute_level_set_from_phase_fraction()
{
  VectorType level_set_owned(this->locally_owned_dofs,
                                           mpi_communicator);
  for (auto p : this->locally_owned_dofs)
    {
      level_set_owned[p] = 2.0*solution[p] - 1.0;
    }
  level_set = level_set_owned;  
}

template <int dim>
void 
AdvectionProblem<dim>::compute_sign_distance()
{   
  
  pcout << "In redistancation" << std::endl;
  
  unsigned int n;
  
  if constexpr (dim == 3)
    n = 4;
  if constexpr (dim == 2)
    n = 2;
  
  // Store (or precompute) vertex to cells connectivity
  const auto vetex_2_cells = grid_tools_cache.get_vertex_to_cell_map();
  
  // Local distance vetors initialization
  VectorType distance_owned(this->locally_owned_dofs,
                                           mpi_communicator);
  VectorType signed_distance_owned(this->locally_owned_dofs,
                                           mpi_communicator);           
                                       
  // Should be initialized to the max distance that we want to redistantiate
  distance_owned = 2;
                              
  // 
  std::unordered_set<types::global_dof_index> dofs_in_interface_halo;
                                           
  // Dummy quadrature to have the intersection points 
  const QGaussLobatto<1> quadrature_1D(n);
  
  // Update flag for the non-matching fe values
  NonMatching::RegionUpdateFlags region_update_flags;
  region_update_flags.surface = update_quadrature_points |
                                update_normal_vectors;
  
  // non-matching fe values
  NonMatching::FEValues<dim> non_matching_fe_values(fe_collection,
                                                    quadrature_1D,
                                                    region_update_flags,
                                                    mesh_classifier,
                                                    dof_handler,
                                                    level_set);
  // Interface rescontruction visualization
  Visualization<dim> intersection_data_out;  
  
  intersection_data_out.build_patches(dof_handler, non_matching_fe_values);
  
  intersection_data_out.write_vtu_with_pvtu_record("output/",
                                      "interface",
                                      1,
                                      MPI_COMM_WORLD,
                                      3);
  
  // DoF coordinates
  std::map< types::global_dof_index, Point<dim >> dof_support_points = DoFTools::map_dofs_to_support_points(mapping,
                                       dof_handler);
  
  
  pcout << "In intersection" << std::endl;
  
  // Loop to identify intersected cells and compute the intersection points (interface reconstruction)
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

      cell->get_dof_indices(dof_indices);
      
      non_matching_fe_values.reinit(cell);
    
      // Surface (interface reconstruction in the intersected cell) fe values (non-empty if the cell is indeed intersected)
      const std::optional<NonMatching::FEImmersedSurfaceValues<dim>>
        &surface_fe_values = non_matching_fe_values.get_surface_fe_values();
      
      // If the cell is intersected, reconstruct the interface in it 
      if (surface_fe_values)
      {
        const unsigned int cell_index = cell->global_active_cell_index();
        
        std::vector<Point<dim>> cells_intersection_point;
                  
        cells_intersection_point.reserve(n);
        
        // Rescontruct interface using Gauss Lobatto quadrature first points (because Gauss Lobatto includes extremities as the first points -> first 2 in 2D to define the intersection line, first 3 in 3D to define the intersection plane)
        for (const unsigned int q :
             surface_fe_values->quadrature_point_indices())
          {
            const Point<dim> point = surface_fe_values->quadrature_point(q);
            
            cells_intersection_point.push_back(point);
            // intersection_point_set.insert(point);
          }
          
        // Store the interface reconstruction 
        intersection_point[cell_index] = cells_intersection_point;
      }
    }
  }
  
  //  Loop to compute distance for the Dofs in of the intersected cells and the first neighbor cells (cells that have a vertice shared with an intersected cell)
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      const unsigned int cell_index = cell->global_active_cell_index();
      
      const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
      
      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

      // If the cell is not stored in the intersection_point map, it means it is not intersected. So no distance computation for now.
      if (intersection_point.find(cell_index) == intersection_point.end())
      {
        continue;
      }
      
      // Pre-store intersection points
      std::vector<Point<dim>> cells_intersection_point =  intersection_point.at(cell_index);
      
      // This will not work in 3D
      const Point<dim> point_0 = cells_intersection_point[0];
      const Point<dim> point_1 = cells_intersection_point[1];
      
      const unsigned int vertices_per_cell =
            GeometryInfo<dim>::vertices_per_cell;
      
      // Loop on the vertices of the cell
      for (unsigned int i = 0; i < vertices_per_cell; i++)
      {
        unsigned int vextex_index = cell->vertex_index(i);
        
        // Loop vertice's neighbor cells
        for (const auto &neighbor_cell_acc : vetex_2_cells[vextex_index])
        {
          const auto neighbor_cell = neighbor_cell_acc->as_dof_handler_iterator(dof_handler);
          
          const unsigned int neighbor_cell_index = neighbor_cell->global_active_cell_index();
      
          const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
      
          std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
      
          neighbor_cell->get_dof_indices(dof_indices);
          
          intersection_halo_cell.insert(neighbor_cell_index);
      
          // Compute distance of the cell's dof 
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            dofs_in_interface_halo.insert(dof_indices[i]);
            
            const Point<dim> y = dof_support_points.at(dof_indices[i]);
          
            double D = compute_point_2_interface_min_distance(point_0, point_1, y);
            const double previous_D = distance_owned[dof_indices[i]];
          
            distance_owned[dof_indices[i]] = std::min(std::abs(previous_D), std::abs(D));
            dofs_location_status.insert(dof_indices[i]);
          }
      
        }
      }
      
    }
  }
  
  // Compute signed distance
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      const unsigned int cell_index = cell->global_active_cell_index();
      
      const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
      
      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
      
      cell->get_dof_indices(dof_indices);
      
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        const double level_set_value = level_set[dof_indices[i]];
      
        signed_distance_owned[dof_indices[i]] = distance_owned[dof_indices[i]]*sgn(level_set_value);
      }
    }
  }
  
  // Exchange
  signed_distance = signed_distance_owned;
  distance = distance_owned;
  
  // for (const auto &cell : dof_handler.active_cell_iterators())
  // {
  //   if (cell->is_locally_owned())
  //   {
  //     const unsigned int cell_index = cell->global_active_cell_index();
  // 
  //     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  // 
  //     std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  // 
  //     cell->get_dof_indices(dof_indices);
  // 
  // 
  // 
  //     for (unsigned int i = 0; i < dofs_per_cell; ++i)
  //     {
  //       const Point<dim> x_J_real = dof_support_points.at(dof_indices[i]);
  // 
  //       distance_owned[dof_indices[i]] = 0.5;
  //     }
  //   }
  // }
  // 
  // distance = distance_owned;
  
  
  pcout << "In redistancation from the rest of the mesh" << std::endl;
  const unsigned int n_opposite_faces_per_dofs = dim;

  
  // Compute the rest of the mesh
  bool change = true;
  while (change)
  {
    FEPointEvaluation<1, dim> fe_point_evaluation(mapping, fe, update_values | update_gradients | update_jacobians);
    
    change = false;
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
      {
        
        const unsigned int cell_index = cell->global_active_cell_index();
        
        const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
        
        std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
        
        cell->get_dof_indices(dof_indices);
        
        std::vector<double> cell_dof_values(dofs_per_cell);
        cell->get_dof_values(distance, cell_dof_values.begin(), cell_dof_values.end());
        
        
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          if (dofs_in_interface_halo.find(dof_indices[i]) != dofs_in_interface_halo.end())
          {
            continue;
          }

          // Get opposite faces 
          std::vector<unsigned int> dof_opposite_faces(n_opposite_faces_per_dofs);
          
          get_dof_opposite_faces<dim>(i, dof_opposite_faces);
          
          const Point<dim> x_J_real = dof_support_points.at(dof_indices[i]);
          
          
          // Loop on opposite faces
          for (unsigned int j = 0; j < n_opposite_faces_per_dofs; ++j)
          {
            const auto opposite_face = cell->face(dof_opposite_faces[j]);
            
            
            
            Point<dim> x_n_ref= transform_ref_face_point_to_ref_cell<dim>(Point<dim-1>(0.5),dof_opposite_faces[j]);
            Point<dim> x_n_real;
            
            double correction_norm = 1.0;
            int outside_check = 0;
            
            while (correction_norm > 1e-11 && outside_check<3)
            {
              const double perturbation = 0.1;
                      
              std::vector<Point<dim>> stencil_ref = compute_numerical_jacobian_stencil<dim>(x_n_ref, dof_opposite_faces[j], perturbation);
              
              fe_point_evaluation.reinit(cell, stencil_ref);
              
              fe_point_evaluation.evaluate(cell_dof_values, EvaluationFlags::gradients);
              
              std::vector<Point<dim>> stencil_real(2*dim - 1);
              std::vector<Tensor<1, dim>> distance_gradients(2*dim - 1);
              std::vector<DerivativeForm<1, dim, dim>> cell_transformation_jacobians(2*dim - 1);
              std::vector<DerivativeForm<1, dim-1, dim>> face_transformation_jacobians(2*dim-1);
              
              for (unsigned int k = 0; k < 2*dim - 1; k++)
              {
                stencil_real[k] = fe_point_evaluation.quadrature_point(k);
                distance_gradients[k] = fe_point_evaluation.get_gradient(k);
                cell_transformation_jacobians[k] = fe_point_evaluation.jacobian(k);
                get_face_transformation_jacobians(cell_transformation_jacobians[k], dof_opposite_faces[j], face_transformation_jacobians[k]);
              }
              // FEFaceValues<dim> fe_face_values(mapping, fe, Quadrature<dim-1>(Point<dim-1>(0.5)), update_gradients | update_jacobians | update_jacobian_grads | update_normal_vectors | update_quadrature_points);
              // 
              // fe_face_values.reinit(cell, opposite_face);
              // 
              // std::vector<Point<dim>> stencil_real = fe_face_values. get_quadrature_points();
              // 
              // std::vector<Tensor<1, dim>> distance_gradients(2*dim - 1);
              // std::vector<Tensor<1, dim>> distance_gradient(1);
              
              // fe_face_values.get_function_gradients(distance,distance_gradient);

              
              // 
              // const std::vector<DerivativeForm<1, dim, dim>> cell_transformation_jacobians = fe_face_values.get_jacobians();
              // 
              // const std::vector<Tensor<1, dim>> normal_vectors = fe_face_values.get_normal_vectors();
              // 
              // std::vector<DerivativeForm<1, dim-1, dim>> face_transformation_jacobians(2*dim-1);
              // 
              // for (unsigned int k = 0; k < 2*dim-1; ++k)
              // {
              //   const DerivativeForm<1, dim, dim> cell_transformation_jacobian = cell_transformation_jacobians[k]; 
              // 
              //   get_face_transformation_jacobians(cell_transformation_jacobians[k], dof_opposite_faces[j], face_transformation_jacobians[k]);
              // }
              // 
              LAPACKFullMatrix<double> jacobian_matrix(dim-1,dim-1);
              compute_numerical_jacobians(stencil_real, x_J_real, distance_gradients, face_transformation_jacobians, perturbation, jacobian_matrix);
              
              
              const Tensor<1,dim> x_n_to_x_J_real = x_J_real - stencil_real[0];
              
              Tensor<1, dim-1> residual_n;
              
              compute_residual(x_n_to_x_J_real, distance_gradients[0], face_transformation_jacobians[0], residual_n);
              
              Vector<double> residual_n_vec(dim-1);
              
              residual_n.unroll(residual_n_vec);
              
              jacobian_matrix.set_property(LAPACKSupport::general);
              jacobian_matrix.compute_lu_factorization();
              
              residual_n_vec *= -1.0;
              jacobian_matrix.solve(residual_n_vec);
              
              
              
              // for (unsigned int k = 0; k < dim-1; ++k)
              // {
              //   correction[k] = residual_n_vec[k];
              // }
              correction_norm = residual_n_vec.l2_norm();
              
              Tensor<1,dim> correction = transform_ref_face_correction_to_ref_cell<dim, Vector<double>>(residual_n_vec,dof_opposite_faces[j]);
              
              Point<dim> x_n_p1_ref = stencil_ref[0] + correction;
              
              
              // Tensor<1,dim> correction_real;
              
              
              // fe_point_evaluation.evaluate(cell_dof_values, EvaluationFlags::gradients);
              
              // 
              // for (unsigned int k = 0; k < dim; ++k)
              // {
              //   for (unsigned int l = 0; l < dim-1; ++l)
              //   {
              //     correction_real[k] = correction[l]*face_transformation_jacobians[0][k][l];
              //   }
              // }
              
              // 
              bool check = true;
              
              double relaxation = 1.0;
              while (check)
              {
                x_n_p1_ref = stencil_ref[0] + relaxation*correction;
                
                check = false;
                for (unsigned int k = 0; k < dim; ++k)
                {
                  if (x_n_p1_ref[k] > 1.0 + 1e-8|| x_n_p1_ref[k] < 0.0 - 1e-8)
                  {
                    check = true;
                  }
                }
                if (check)
                  outside_check +=1;
                
                relaxation *= 0.999;
              }
              
              
              std::vector<Point<dim>> x_n_p1_ref_vec = {x_n_p1_ref}; 
              
              fe_point_evaluation.reinit(cell, x_n_p1_ref_vec);
              
              Point<dim> x_n_p1_real = fe_point_evaluation.quadrature_point(0);
              
              
              // x_n_p1_ref = stencil_ref[0] + relaxation*correction;
              // 
              // for (unsigned int k = 0; k < dim; ++k)
              // {
              //   for (unsigned int l = 0; l < dim-1; ++l)
              //   {
              //     correction_real[k] = relaxation*correction[l]*face_transformation_jacobians[0][k][l];
              //   }
              // }
              // 
              // x_n_p1_real = stencil_real[0] + correction_real;
              // 
              x_n_ref = x_n_p1_ref;
              x_n_real = x_n_p1_real;
              
              
              // if (check)
              // {
                
                
                // Tensor<1,dim> x_n_p1_to_x_J_real = x_J_real - x_n_p1_real;
                // Tensor<1, dim-1> residual_n_p1;
                // 
                // // this works only for Q1 in 2D and linear mapping
                // compute_residual(x_n_p1_to_x_J_real, distance_gradients[0],face_transformation_jacobians[0],residual_n_p1);
                // 
                // // outside the face, adaptively relax the scheme
                // while (residual_n_p1.norm() > residual_n.norm())
                // {
                //   relaxation *= 0.5;
                //   std::cout << "relaxation " << relaxation << std::endl;
                // 
                //   x_n_p1_real = stencil_real[0] + relaxation*correction_real;
                // 
                //   std::cout << "x_n_p1_real " << x_n_p1_real << std::endl;
                // 
                //   x_n_p1_to_x_J_real = x_J_real - x_n_p1_real;
                // 
                //   // this works only for Q1 in 2D and linear mapping
                //   compute_residual(x_n_p1_to_x_J_real, distance_gradients[0],face_transformation_jacobians[0],residual_n_p1);
                // 
                // }
              // }
              // x_n_p1_ref = stencil_ref[0] + relaxation*correction;
              // x_n_p1_real = stencil_real[0] + relaxation*correction_real;
          }
          
          const Tensor<1,dim> x_n_to_x_J_real = x_J_real - x_n_real;
          
          fe_point_evaluation.evaluate(cell_dof_values, EvaluationFlags::values);
          
          double distance_value_at_x_n = fe_point_evaluation.get_value(0);
          
          // FEFaceValues<dim> fe_face_values(mapping, fe, Quadrature<dim-1>(x_n_ref), update_values);
          
          // fe_face_values.reinit(cell, opposite_face);
          
          // std::vector<double> distance_value_at_x_n(1);
          // fe_face_values.get_function_values(distance, distance_value_at_x_n);
          

          double approx_distance = compute_distance(x_n_to_x_J_real, distance_value_at_x_n);
          
          if (cell_dof_values[i] > (approx_distance + 1e-8))
          {
            change = true;
            distance_owned[dof_indices[i]] = approx_distance;
          }
        }
        }
      }
    }
    
    // Compute signed distance
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
      {
        const unsigned int cell_index = cell->global_active_cell_index();
        
        const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
        
        std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
        
        cell->get_dof_indices(dof_indices);
        
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const double level_set_value = level_set[dof_indices[i]];
        
          signed_distance_owned[dof_indices[i]] = distance_owned[dof_indices[i]]*sgn(level_set_value);
        }
      }
    }
    
    // Exchange
    signed_distance = signed_distance_owned;
    distance = distance_owned;
  }
  
  // const Tensor<1,dim> d = point_1 - point_0;
  // 
  // // compute sign distance of intersected cell's dof 
  // for (unsigned int i = 0; i < dofs_per_cell; ++i)
  // {
  //   const Point<dim> y = dof_support_points.at(dof_indices[i]);
  // 
  //   double D = compute_point_2_interface_min_distance(point_0, point_1, y);
  //   const double previous_D = distance_owned[dof_indices[i]];
  //   const double level_set_value = level_set[dof_indices[i]];
  // 
  //   distance_owned[dof_indices[i]] = std::min(std::abs(previous_D), std::abs(D))*sgn(level_set_value);
  //   dofs_location_status.insert(dof_indices[i]);
  // }

}

template <int dim>
void AdvectionProblem<dim>::output_results(const int time_iteration, const double time) const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.add_data_vector(level_set, "level_set");
  data_out.add_data_vector(signed_distance, "signed_distance");
  
  // data_out.add_data_vector(location, "location");
  
  // data_out.set_cell_selection(
  //       [this](const typename Triangulation<dim>::cell_iterator &cell) {
  //         return cell->is_active() &&
  //                mesh_classifier.location_to_level_set(cell) ==
  //                  NonMatching::LocationToLevelSet::intersected;
  //       });
  data_out.build_patches();

  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::CompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);

  const std::string filename = "solution.vtu";
  // std::ofstream     output(filename);
  // data_out.write_vtu(output);
  data_out.write_vtu_with_pvtu_record("output/",
                                      "solution",
                                      time_iteration,
                                      MPI_COMM_WORLD,
                                      3);
  // pcout << "Solution written to " << filename << std::endl;
  std::ofstream out("grid-2.svg");
  GridOut       grid_out;
    grid_out.write_svg(triangulation, out);
}

template <int dim>
void AdvectionProblem<dim>::run()
{
  pcout << "Bonjour from run" << std::endl;
  
  Point<dim> p_0 = Point<dim>();
  p_0[0] = -1;
  for (unsigned int i = 1; i < dim; ++i)
    p_0[i] = -1;
    
  Point<dim> p_1 = Point<dim>();
  p_1[0] = 1;
  for (unsigned int i = 1; i < dim; ++i)
    p_1[i] = 1;
    
  std::vector< unsigned int > repetitions(dim);
  repetitions[0] = 1;
  for (unsigned int i = 1; i < dim; ++i)
    repetitions[i] = 1;
    
  GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, p_0, p_1);
  triangulation.refine_global(7);
            
  pcout << "Bonjour from after triangulation" << std::endl;
  // initial time step
  setup_system();
  
  pcout << "Bonjour from after setup_system()" << std::endl;
  
  
  set_initial_conditions();
  
  
  unsigned int it = 0;
  double time = 0.0;
  double final_time = 0.01;
  
  compute_level_set_from_phase_fraction();
  
  mesh_classifier.reclassify();


  compute_sign_distance();


    
  output_results(it,time);
  
  
  pcout << "Bonjour from after set_initial_conditions()" << std::endl;
  
  while (time < final_time) 
  {
    it += 1;
    time += this->dt;
    
    
    
    assemble_system();
    
    pcout << "time = " << time << std::endl;

    solve();
    
    compute_level_set_from_phase_fraction();
    
    output_results(it, time);
    
    // 
    // previous_solution = solution;
    this->previous_solution = this->solution; 
  } 
  
  

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


