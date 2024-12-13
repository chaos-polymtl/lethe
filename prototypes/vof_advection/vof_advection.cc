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
#include <deal.II/distributed/solution_transfer.h>

#include <deal.II/base/work_stream.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table_handler.h>


#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/non_matching/fe_immersed_values.h>
#include <deal.II/non_matching/fe_values.h>
#include <deal.II/non_matching/mesh_classifier.h>

#include <fstream>
#include <iostream>

using namespace dealii;

namespace ChangeVectorTypes
{
  template <typename number>
  void copy(TrilinosWrappers::MPI::Vector                    &out,
            const LinearAlgebra::distributed::Vector<number> &in)
  {
    LinearAlgebra::ReadWriteVector<double> rwv(out.locally_owned_elements());
    rwv.import_elements(in, VectorOperation::insert);
    out.import_elements(rwv, VectorOperation::insert);
  }

  template <typename number>
  void copy(LinearAlgebra::distributed::Vector<number> &out,
            const TrilinosWrappers::MPI::Vector        &in)
  {
    LinearAlgebra::ReadWriteVector<double> rwv;
    rwv.reinit(in);
    out.import_elements(rwv, VectorOperation::insert);
  }
}

template <int dim>
inline double compute_point_2_interface_min_distance(const std::vector<Point<dim>> triangle,const Point<dim> &point)
{
  double D;
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
      
      const double a   = e_0.norm_square();
      const double b   = scalar_product(e_0, e_1);
      const double c   = e_1.norm_square();
      const double d = scalar_product(e_0, vector_to_plane);
      const double e = scalar_product(e_1, vector_to_plane);
      
      const double det = abs(a * c - b * b);
      
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
      const Tensor<1,dim> d = point_1 - point_0;
      
      const double t_bar = d*(point-point_0)/(d.norm()*d.norm());
      
      if (t_bar <= 0.0)
      {
        const Tensor<1,dim> point_minus_p0 = point-point_0;
        D = point_minus_p0.norm();
      }
      else if (t_bar >= 1.0)
      {
        const Tensor<1,dim> point_minus_p1 = point-point_1;
        D = point_minus_p1.norm();
      }
      else
      {
        const Tensor<1,dim> projection = point-(point_0+t_bar*d);
        D = projection.norm();
      }
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
  const double period = 2.0;
  Tensor<1, dim> value;
  value[0] = -(Utilities::pow(sin(numbers::PI*p[0]),2)*sin(2*numbers::PI*p[1])*cos(numbers::PI*this->get_time()/period));
  value[1] = Utilities::pow(sin(numbers::PI*p[1]),2)*sin(2*numbers::PI*p[0])*cos(numbers::PI*this->get_time()/period);
  
  // const double period = 1.0;
  // value[0] = cos(numbers::PI*this->get_time()/period);
  
  return value;
}

template <int dim, typename VectorType = Vector<double>>
class LocalCellWiseFunction : public Function<dim>
{
public:
  LocalCellWiseFunction();
  
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
  FE_Q<dim> element;
  
  unsigned int n_local_dof;
  
  Vector<typename VectorType::value_type> local_dof_values;
};

template <int dim, typename VectorType>
LocalCellWiseFunction<dim, VectorType>::LocalCellWiseFunction()
  : element(1)
{
  this->n_local_dof = element.dofs_per_cell;
} 

template <int dim, typename VectorType>
void
LocalCellWiseFunction<dim, VectorType>::set_active_cell(const VectorType &in_local_dof_values)
{
  n_local_dof = element.dofs_per_cell;
  local_dof_values = in_local_dof_values;
} 

template <int dim, typename VectorType>
double
LocalCellWiseFunction<dim, VectorType>::value(
  const Point<dim>  &point, const unsigned int component) const
{
  double value = 0;
  for (unsigned int i = 0; i < n_local_dof; ++i)
    value += local_dof_values[i] *
             element.shape_value_component(i, point, component);

  return value;
} 

template <int dim, typename VectorType>
Tensor<1, dim>
LocalCellWiseFunction<dim, VectorType>::gradient(
  const Point<dim>  &point, const unsigned int component) const
{
  Tensor<1, dim> gradient;
  for (unsigned int i = 0; i < n_local_dof; ++i)
    gradient += local_dof_values[i] *
             element.shape_grad_component(i, point, component);

  return gradient;
} 

template <int dim, typename VectorType>
SymmetricTensor<2, dim>
LocalCellWiseFunction<dim, VectorType>::hessian(
  const Point<dim>  &point, const unsigned int component) const
{
  Tensor<2, dim> hessian;
  for (unsigned int i = 0; i < n_local_dof; ++i)
    hessian += local_dof_values[i] *
             element.shape_grad_grad_component(i, point, component);

  return symmetrize(hessian);
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
  
  Point<dim> center = Point<dim>(0.5,0.75);
  Tensor<1,dim> dist = center - p;
  
  return 0.5+0.5*std::tanh((0.15-dist.norm())/0.016);
  
  
  // Point<dim> bl = Point<dim>(0.4,0.4);
  // Point<dim> tr = Point<dim>(0.6,0.6);
  // 
  // Tensor<1,dim> dist_bl = bl - p;
  // Tensor<1,dim> dist_tr = p - tr;
  // 
  // Tensor<1,dim> dist_border;
  // dist_border[0] = std::max(dist_bl[0],dist_tr[0]);
  // dist_border[1] = std::max(dist_bl[1],dist_tr[1]);
  // 
  // Tensor<1,dim> dist_corner;
  // dist_corner[0] = std::max(0.0,dist_border[0]);
  // dist_corner[1] = std::max(0.0,dist_border[1]);
  // 
  // double dist = dist_corner.norm() + std::min(0.0, std::max(dist_border[0],dist_border[1]));
  // return 0.5+0.5*std::tanh(-1.0*dist/0.016);
  
  // return 0.5+0.5*std::tanh((0.1-abs(dist[0]))/0.016);
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
  void compute_sign_distance(unsigned int time_iteration, double global_volume);
  void compute_volume();
  double compute_cell_wise_volume(const typename DoFHandler<dim>::active_cell_iterator &cell, Vector<double> cell_dof_values, const double corr, const BoundingBox<dim> &unit_box, LocalCellWiseFunction<dim> &level_set_function, NonMatching::QuadratureGenerator<dim> &quadrature_generator);
  
  void refine_grid(const unsigned int max_grid_level, const unsigned int min_grid_level);
  void output_results(const int time_iteration) const;
  
  
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
  IndexSet locally_active_dofs;
  
  
  VectorType locally_relevant_solution;
  VectorType previous_solution;
  VectorType location;
  // VectorType cell_wise_volume;
  
  LinearAlgebra::distributed::Vector<double> signed_distance;
  LinearAlgebra::distributed::Vector<double> distance;
  LinearAlgebra::distributed::Vector<double> distance_with_ghost;
  LinearAlgebra::distributed::Vector<double> volume_correction;
  
  
  std::map<types::global_cell_index,std::vector<Point<dim>>> intersection_point;
  
  std::map<types::global_cell_index,std::vector<Point<dim>>> interface_reconstruction_vertices;
  
  std::map<types::global_cell_index,std::vector<CellData<dim-1>>> interface_reconstruction_cells;
  
  std::map<types::global_cell_index,Point<dim>> intersection_cell;
  
  std::set<types::global_cell_index> intersection_halo_cell;
  std::set<types::global_dof_index> dofs_location_status;
  
    
  VectorType system_rhs;
  
  double dt;
  double time = 0.0;
  
  MPI_Comm           mpi_communicator;
  ConditionalOStream pcout;
  TableHandler table_volume_monitoring;
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
  , dt(0.5/300.0)
  , mpi_communicator(MPI_COMM_WORLD)
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
{
  fe_collection.push_back(FE_Q<dim>(1));
  
}
  
template <int dim>
void AdvectionProblem<dim>::setup_system()
{
  // fe_collection.push_back(FE_Q<dim>(1));
      // fe_collection.push_back(FE_Nothing<dim>());
  
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // const NonMatching::LocationToLevelSet cell_location =
      //   mesh_classifier.location_to_level_set(cell);
      if (!cell->is_locally_owned())
        continue;
      cell->set_active_fe_index(ActiveFEIndex::lagrange);
    }
        
  dof_handler.distribute_dofs(fe_collection);
  
  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  n_locally_owned_dofs    = dof_handler.n_locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
  locally_active_dofs = DoFTools::extract_locally_active_dofs(dof_handler);
  
  locally_relevant_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  previous_solution.reinit(locally_owned_dofs,locally_relevant_dofs, mpi_communicator);
  
  level_set.reinit(locally_owned_dofs,locally_relevant_dofs, mpi_communicator);
  location.reinit(locally_owned_dofs, mpi_communicator);
  // cell_wise_volume.reinit(locally_owned_dofs, mpi_communicator);
  
  
  signed_distance.reinit(locally_owned_dofs,locally_active_dofs, mpi_communicator);
  distance.reinit(locally_owned_dofs,locally_active_dofs, mpi_communicator);
  volume_correction.reinit(locally_owned_dofs,locally_active_dofs, mpi_communicator);
  distance_with_ghost.reinit(locally_owned_dofs,locally_active_dofs, mpi_communicator);
  
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
                     
  constraints.clear();
  constraints.reinit(locally_owned_dofs,locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler,
                                          constraints);
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
    this->system_matrix,
    this->system_rhs);
    
}

template <int dim>
void
AdvectionProblem<dim>::set_initial_conditions()
{
  
  VectorType completely_distributed_solution(this->locally_owned_dofs, mpi_communicator);
  VectorTools::interpolate(this->mapping, this->dof_handler,
                           InitialConditions<dim>(),
                           completely_distributed_solution);
  this->constraints.distribute(completely_distributed_solution);
  
  this->locally_relevant_solution = completely_distributed_solution;
  this->previous_solution = this->locally_relevant_solution; 
}

template <int dim>
void AdvectionProblem<dim>::solve()
{
  SolverControl               solver_control(1000,
                               1e-12 * system_rhs.l2_norm() );
  TrilinosWrappers::SolverGMRES solver(solver_control);
  
  TrilinosWrappers::PreconditionILU                 preconditioner;
  TrilinosWrappers::PreconditionILU::AdditionalData data_ilu;
  
  preconditioner.initialize(system_matrix, data_ilu);
  
  VectorType completely_distributed_solution(this->locally_owned_dofs, mpi_communicator);  
  
  solver.solve(system_matrix,
              completely_distributed_solution,
              system_rhs,
              preconditioner);

  VectorType residual(locally_owned_dofs,  mpi_communicator);
  
  system_matrix.vmult(residual, completely_distributed_solution);
  // residual -= system_rhs;
  pcout << "   Iterations required for convergence: "
            << solver_control.last_step() << '\n'
            << "   Max norm of residual:                "
            << residual.linfty_norm() << '\n';

  constraints.distribute(completely_distributed_solution);
  locally_relevant_solution = completely_distributed_solution;
}

template <int dim>
void 
AdvectionProblem<dim>::refine_grid(const unsigned int max_grid_level, const unsigned int min_grid_level)
{
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_trans(dof_handler);
    
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  
  KellyErrorEstimator<dim>::estimate(mapping, dof_handler,
                                         QGauss<dim-1>(fe.degree+1),
                                         typename std::map<types::boundary_id, const Function<dim, double> *>(),
                                         locally_relevant_solution,
                                         estimated_error_per_cell,
                                         ComponentMask(),
                                         nullptr,
                                         0,
                                         triangulation.locally_owned_subdomain());
                                         
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                        estimated_error_per_cell,
                                                        0.9,
                                                        0.01);         
  if (triangulation.n_levels() > max_grid_level)
    for (auto &cell :
         triangulation.active_cell_iterators_on_level(max_grid_level))
      cell->clear_refine_flag();
  for (auto &cell :
       triangulation.active_cell_iterators_on_level(min_grid_level))
    cell->clear_coarsen_flag();
    

  
  triangulation.prepare_coarsening_and_refinement();
  
  solution_trans.prepare_for_coarsening_and_refinement(locally_relevant_solution);
  
  triangulation.execute_coarsening_and_refinement();
  
  setup_system();
  
  VectorType tmp_solution(this->locally_owned_dofs, mpi_communicator);
  
  solution_trans.interpolate(tmp_solution);
  
  constraints.distribute(tmp_solution);
  
  locally_relevant_solution = tmp_solution;
  
}

template <int dim>
void 
AdvectionProblem<dim>::compute_phase_fraction_from_level_set()
{
  VectorType solution_owned(this->locally_owned_dofs,
                                           mpi_communicator);
  for (auto p : this->locally_owned_dofs)
    {
      const double signed_dist = signed_distance[p];
      solution_owned[p] = 0.5+0.5*std::tanh(signed_dist/0.016);
    }
    constraints.distribute(solution_owned);
    
    locally_relevant_solution = solution_owned;
    
}

template <int dim>
void 
AdvectionProblem<dim>::compute_level_set_from_phase_fraction()
{
  VectorType level_set_owned(this->locally_owned_dofs,
                                           mpi_communicator);
  for (auto p : this->locally_owned_dofs)
    {
      const double phase = locally_relevant_solution[p];
      // level_set_owned[p] = log10(std::max(phase,1e-8)/std::max(1.0-phase,1e-8));
      double phase_sign = sgn(phase-0.5);
      level_set_owned[p] = 0.016*std::atanh(phase_sign*std::min(abs(phase-0.5)/0.5,1.0-1e-16));
      // level_set_owned[p] = phase-0.5;
    }
  constraints.distribute(level_set_owned);
  
  level_set = level_set_owned; 
   
}

template <int dim>
double 
AdvectionProblem<dim>::compute_cell_wise_volume(const typename DoFHandler<dim>::active_cell_iterator &cell, Vector<double> cell_dof_values, const double corr, const BoundingBox<dim> &unit_box, LocalCellWiseFunction<dim> &level_set_function, NonMatching::QuadratureGenerator<dim> &quadrature_generator)
{
  for (unsigned int j = 0; j < fe.n_dofs_per_cell(); j++)
  {
    cell_dof_values[j] += corr;
  }
  
  level_set_function.set_active_cell(cell_dof_values);
  quadrature_generator.generate(level_set_function, unit_box);
  
  const Quadrature<dim> inside_quadrature = quadrature_generator.get_outside_quadrature();
  
  if (inside_quadrature.size() == 0)
    return 0.0;
    
  FEValues<dim> inside_fe_values(mapping, fe, inside_quadrature, update_values | update_quadrature_points | update_JxW_values);
    
  inside_fe_values.reinit(cell);
  std::vector<double> inside_JxW = inside_fe_values.get_JxW_values();

  double inside_cell_volume = 0.0;
  for (unsigned int q = 0; q < inside_quadrature.size(); q++)
  {
    inside_cell_volume += inside_JxW[q];
  }
  
  return inside_cell_volume;
}

template <int dim>
void 
AdvectionProblem<dim>::compute_sign_distance(unsigned int time_iteration, double global_volume)
{   
  
  pcout << "In redistancation" << std::endl;
  
  intersection_point.clear();
  interface_reconstruction_vertices.clear();
  interface_reconstruction_cells.clear();
  intersection_cell.clear();
  intersection_halo_cell.clear();
  dofs_location_status.clear();
  
  unsigned int n;
  
  if constexpr (dim == 2)
    n = 2;
  if constexpr (dim == 3)
    n = 4;

  
  std::unordered_set<types::global_dof_index> dofs_in_interface_halo;
                                           
  // Dummy quadrature to have the intersection points 
  const QGaussLobatto<1> quadrature_1D(n);
  
  // Update flag for the non-matching fe values
  NonMatching::RegionUpdateFlags region_update_flags;
  region_update_flags.surface = update_quadrature_points | 
                                update_normal_vectors;
  
  region_update_flags.inside = update_quadrature_points | 
                                update_JxW_values;
  // non-matching fe values
  NonMatching::FEValues<dim> non_matching_fe_values(fe_collection,
                                                    quadrature_1D,
                                                    region_update_flags,
                                                    mesh_classifier,
                                                    dof_handler,
                                                    level_set);
                                                    
                                                    
  GridTools::MarchingCubeAlgorithm<dim, VectorType> marching_cube(mapping,
                                                         fe, // todo
                                                         1,
                                                         1e-10);
  
  // double global_volume = 0.0;
  // for (const auto &cell : dof_handler.active_cell_iterators())
  // {
  //   if (cell->is_locally_owned())
  //   {
  // 
  //     non_matching_fe_values.reinit(cell);
  // 
  //     const std::optional<FEValues<dim>>
  //       &inside_fe_values = non_matching_fe_values.get_inside_fe_values();
  // 
  // 
  //     if (!inside_fe_values)
  //     {
  //       continue;
  //     }
  // 
  //     std::vector<double> JxW_inside = inside_fe_values->get_JxW_values();
  // 
  //     for (const unsigned int q :
  //          inside_fe_values->quadrature_point_indices())
  //       {
  //         global_volume += inside_fe_values->JxW(q);
  //       }
  // 
  //   }
  // }
  // global_volume = Utilities::MPI::sum(global_volume, mpi_communicator);
  // 
  // pcout << "global_volume = " << global_volume << std::endl;
  
  // Interface rescontruction visualization
  Visualization<dim> intersection_data_out;  
  
  intersection_data_out.build_patches(dof_handler, non_matching_fe_values);
  
  intersection_data_out.write_vtu_with_pvtu_record("output/",
                                      "interface",
                                      time_iteration,
                                      MPI_COMM_WORLD,
                                      3);
  
  // DoF coordinates
  std::map< types::global_dof_index, Point<dim >> dof_support_points = DoFTools::map_dofs_to_support_points(mapping,
                                       dof_handler);
  
  pcout << "In signed distance computation" << std::endl;
  
  // Local distance vetors initialization
  for (auto p : this->locally_active_dofs)
    {
      distance(p) = 0.08;
      distance_with_ghost(p) = 0.08;
    }
    
  const unsigned int n_quad_points = fe.degree + 1;
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  
  // Loop to identify intersected cells and compute the interface reconstruction
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
        
      std::vector<Point<dim>>    surface_vertices;
      std::vector<CellData<dim-1>> surface_cells;
              
      marching_cube.process_cell(cell, level_set, 0.0, surface_vertices, surface_cells);
      
      // If the cell is intersected, reconstruct the interface in it
      if (surface_vertices.size() != 0)
      {
        const unsigned int cell_index = cell->global_active_cell_index();
        
        interface_reconstruction_vertices[cell_index] = surface_vertices;
        interface_reconstruction_cells[cell_index] = surface_cells;
        
        
        std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
        cell->get_dof_indices(dof_indices);
        
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          dofs_location_status.insert(dof_indices[i]);
        }
        
        Triangulation<dim-1,dim> surface_triangulation;
        surface_triangulation.create_triangulation(surface_vertices, surface_cells, {});
        

        for (const auto &surface_cell : surface_triangulation.active_cell_iterators())
          {
            unsigned int surface_cell_n_vertices = surface_cell->n_vertices();
            std::vector<Point<dim>> surface_cell_vertices(surface_cell_n_vertices);
            
            for (unsigned int p = 0; p < surface_cell_n_vertices; p++)
            {
              surface_cell_vertices[p] = surface_cell->vertex(p);
            }
            
            intersection_point[cell_index] = surface_cell_vertices;
          }
      }
    }
  }
  

  //  Loop to compute distance for the Dofs in of the intersected cells and the first neighbor cells (cells that have a vertice shared with an intersected cell)
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      const unsigned int cell_index = cell->global_active_cell_index();
      
      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
      cell->get_dof_indices(dof_indices);
    
      // If the cell is not stored in the intersection_point map, it means it is not intersected. So no distance computation for now.
      if (interface_reconstruction_vertices.find(cell_index) == interface_reconstruction_vertices.end())
      {
        continue;
      }
      
      std::vector<Point<dim>>    surface_vertices = interface_reconstruction_vertices.at(cell_index);
      std::vector<CellData<dim-1>> surface_cells = interface_reconstruction_cells.at(cell_index);
      
      // Create interface recontruction triangulation
      Triangulation<dim-1,dim> surface_triangulation;
      surface_triangulation.create_triangulation(surface_vertices, surface_cells, {});
      
      // Compute minimum distance of the cell's dof to the interface reconstruction
      for (const auto &surface_cell : surface_triangulation.active_cell_iterators())
        {
          unsigned int surface_cell_n_vertices = surface_cell->n_vertices();
          std::vector<Point<dim>> surface_cell_vertices(surface_cell_n_vertices);
          
          for (unsigned int p = 0; p < surface_cell_n_vertices; p++)
          {
            surface_cell_vertices[p] = surface_cell->vertex(p);
          }
          
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              dofs_in_interface_halo.insert(dof_indices[i]);
        
              const Point<dim> y = dof_support_points.at(dof_indices[i]);
              double D = compute_point_2_interface_min_distance(surface_cell_vertices, y);
              
              distance(dof_indices[i]) = std::min(std::abs(distance(dof_indices[i])), std::abs(D));
            }
        }
    }
  }
  
  // Exchange and "select" min value between the processes
  distance.compress(VectorOperation::min);
  
  // Update local ghost (distance becomes read only)
  distance.update_ghost_values();    
  
  // Copy distance to distance_with_ghost
  distance_with_ghost = distance;
  
  // Zero out ghost DOFs to regain write functionalities in distance (it becomes write only, that is why we need distance_with_ghost - to read the values in it.)
  distance.zero_out_ghost_values();    

  // Copy the ghost values back in distance (zero_out_ghost_values() puts zeros in ghost DOFs)
  for (auto p : this->locally_active_dofs)
  {
    distance(p) = distance_with_ghost(p);

    const double level_set_value = level_set(p);
    signed_distance(p) = distance(p)*sgn(level_set_value);
  }
    
  // Correct distance to conserve volume
  const BoundingBox<dim> unit_box = create_unit_bounding_box<dim>();
  LocalCellWiseFunction<dim> level_set_function = LocalCellWiseFunction<dim>();
  
  hp::QCollection<1> q_collection;
  q_collection.push_back(QGauss<1>(n_quad_points));
  
  NonMatching::QuadratureGenerator<dim> quadrature_generator = NonMatching::QuadratureGenerator<dim>(q_collection);
  // 
  global_volume = 0.0;
  // 
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
  
      Vector<double> cell_dof_values(dofs_per_cell);
  
      cell->get_dof_values(level_set, cell_dof_values.begin(), cell_dof_values.end());
  
      global_volume += compute_cell_wise_volume(cell, cell_dof_values, 0.0, unit_box, level_set_function, quadrature_generator);
  
    }
  }
  global_volume = Utilities::MPI::sum(global_volume, mpi_communicator);
  
  pcout << "global_volume = " << global_volume << std::endl;
  
  pcout << "Coucou" << std::endl;
  volume_correction = 0.0;
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      const unsigned int cell_index = cell->global_active_cell_index();
      
      // The cell is not intersected, no need to correct the mass
      if (interface_reconstruction_vertices.find(cell_index) == interface_reconstruction_vertices.end())
      {
        continue;
      }
      
      Vector<double> cell_dof_values(dofs_per_cell);
      Vector<double> cell_level_set_dof_values(dofs_per_cell);
      
      cell->get_dof_values(distance_with_ghost, cell_dof_values.begin(), cell_dof_values.end());
      cell->get_dof_values(level_set, cell_level_set_dof_values.begin(), cell_level_set_dof_values.end());
      
      double cell_volume = compute_cell_wise_volume(cell, cell_level_set_dof_values, 0.0, unit_box, level_set_function, quadrature_generator);;  
      
      double cell_size;
      if (dim == 2)
        {
          cell_size = std::sqrt(4. * cell->measure() / M_PI);
        }
      else if (dim == 3)
        {
          cell_size = std::pow(6 * cell->measure() / M_PI, 1. / 3.);
        }
      
      
      for (unsigned int j = 0; j < dofs_per_cell; j++)
      {
        const double level_set_value = cell_level_set_dof_values[j];
        cell_dof_values[j] *= sgn(level_set_value);
      }
      
      double local_corr_nm1 = 0.0;
      double inside_cell_volume_nm1 = compute_cell_wise_volume(cell, cell_dof_values, local_corr_nm1, unit_box, level_set_function, quadrature_generator);
      
      double delta_volume_nm1 = abs(cell_volume - inside_cell_volume_nm1);
      
      double local_corr_n = cell_size*delta_volume_nm1/cell->measure();
      
      // pcout << "Bonjour, the volume error is " << delta_volume_nm1 << std::endl;
      while (abs(delta_volume_nm1) > 1e-10)
      {
        
        double inside_cell_volume_n = compute_cell_wise_volume(cell, cell_dof_values, local_corr_n, unit_box, level_set_function, quadrature_generator);
        
        double delta_volume_n = abs(cell_volume - inside_cell_volume_n);
        
        double delta_volume_prime = (delta_volume_n - delta_volume_nm1)/(local_corr_n - local_corr_nm1);
        
        double local_corr_np1 = local_corr_n - delta_volume_n/delta_volume_prime;
        
        local_corr_nm1 = local_corr_n;
        local_corr_n = local_corr_np1;
        inside_cell_volume_nm1 = inside_cell_volume_n;
        delta_volume_nm1 = delta_volume_n;
        
        // pcout << "Now, the corr is "<< local_corr_np1 << " and the volume error is " << delta_volume_nm1 << std::endl;
        
      }
      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
      cell->get_dof_indices(dof_indices);
      
      double n_cells_per_dofs = 4.0;
      if constexpr (dim == 3)
      {
        n_cells_per_dofs = 8.0;
      }
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        volume_correction(dof_indices[i]) += local_corr_n/n_cells_per_dofs;
        
      }
      
      // pcout << "Bye" << std::endl;
    }
  }
  
  volume_correction.compress(VectorOperation::add);
  volume_correction.update_ghost_values();
  
  // for (auto p : this->locally_active_dofs)
  // {
  //   const double level_set_value = level_set(p);
  // 
  //   distance(p) += volume_correction(p)*sgn(level_set_value);
  // }
  // 
  // distance.update_ghost_values();    
  // distance_with_ghost = distance;
  // distance.zero_out_ghost_values();    
  // 
  // for (auto p : this->locally_active_dofs)
  // {
  //   distance(p) = distance_with_ghost(p);
  // }
  
  double corr_global_volume_nm1 = 0.0;
  
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      
      Vector<double> cell_dof_values(dofs_per_cell);
      Vector<double> cell_level_set_dof_values(dofs_per_cell);
      Vector<double> cell_volume_correction_dof_values(dofs_per_cell);
      
      
      cell->get_dof_values(distance_with_ghost, cell_dof_values.begin(), cell_dof_values.end());
      cell->get_dof_values(level_set, cell_level_set_dof_values.begin(), cell_level_set_dof_values.end());
      cell->get_dof_values(volume_correction, cell_volume_correction_dof_values.begin(), cell_volume_correction_dof_values.end());
      
      
      for (unsigned int j = 0; j < dofs_per_cell; j++)
      {
        const double level_set_value = cell_level_set_dof_values[j];
        cell_dof_values[j] *= sgn(level_set_value);
        cell_dof_values[j] += cell_volume_correction_dof_values[j];
      }
      
      corr_global_volume_nm1 += compute_cell_wise_volume(cell, cell_dof_values, 0.0, unit_box, level_set_function, quadrature_generator);
  
    }
  }
  corr_global_volume_nm1 = Utilities::MPI::sum(corr_global_volume_nm1, mpi_communicator);
  
  pcout << "corr_global_volume_nm1 = " << corr_global_volume_nm1 << std::endl;

  double global_delta_volume_nm1 = abs(corr_global_volume_nm1 - global_volume);
  double C_nm1 = 1.0;
  
  double corr_global_volume_n = corr_global_volume_nm1;
  double C_n = global_delta_volume_nm1/global_volume;
  
  double global_delta_volume_n = 0.0;
  while (global_delta_volume_nm1 > 1e-10)
  {
    corr_global_volume_n = 0.0;
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
      {
        
        Vector<double> cell_dof_values(dofs_per_cell);
        Vector<double> cell_level_set_dof_values(dofs_per_cell);
        Vector<double> cell_volume_correction_dof_values(dofs_per_cell);
        
        
        cell->get_dof_values(distance_with_ghost, cell_dof_values.begin(), cell_dof_values.end());
        cell->get_dof_values(level_set, cell_level_set_dof_values.begin(), cell_level_set_dof_values.end());
        cell->get_dof_values(volume_correction, cell_volume_correction_dof_values.begin(), cell_volume_correction_dof_values.end());
        
        
        for (unsigned int j = 0; j < dofs_per_cell; j++)
        {
          const double level_set_value = cell_level_set_dof_values[j];
          cell_dof_values[j] *= sgn(level_set_value);
          cell_dof_values[j] += C_n*cell_volume_correction_dof_values[j];
        }
        
        corr_global_volume_n += compute_cell_wise_volume(cell, cell_dof_values, 0.0, unit_box, level_set_function, quadrature_generator);
    
      }
    }
    corr_global_volume_n = Utilities::MPI::sum(corr_global_volume_n, mpi_communicator);
    
    pcout << "corr_global_volume_n = " << corr_global_volume_n << std::endl;
    
    global_delta_volume_n = abs(corr_global_volume_n - global_volume);
    
    double global_delta_volume_prime = (global_delta_volume_n - global_delta_volume_nm1)/(C_n - C_nm1);
    
    double C_np1 = C_n - global_delta_volume_n/global_delta_volume_prime;
    
    C_nm1 = C_n;
    C_n = C_np1;
    
    corr_global_volume_nm1 = corr_global_volume_n;
    global_delta_volume_nm1 = global_delta_volume_n;
    
  }
  
  
  for (auto p : this->locally_active_dofs)
  {
    const double level_set_value = level_set(p);
    signed_distance(p) +=  C_n*volume_correction(p);
    distance(p) = abs(signed_distance(p));
  }
  
  distance.update_ghost_values();    
  distance_with_ghost = distance;
  distance.zero_out_ghost_values();    
  
  for (auto p : this->locally_active_dofs)
  {
    distance(p) = distance_with_ghost(p);
  }
  
  unsigned int n_opposite_faces_per_dofs = dim;
  
  // Compute the rest of the mesh
  bool change = true;
  int count = 0;
  while (change)
  {
    count +=1;
    FEPointEvaluation<1, dim> fe_point_evaluation(mapping, fe, update_values | update_gradients | update_jacobians);
    
    change = false;
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
      {
        const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
        
        std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
        
        cell->get_dof_indices(dof_indices);
        std::vector<double> cell_dof_values(dofs_per_cell);
        cell->get_dof_values(distance_with_ghost, cell_dof_values.begin(), cell_dof_values.end());
        
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
            // const auto opposite_face = cell->face(dof_opposite_faces[j]);
            
            Point<dim> x_n_ref= transform_ref_face_point_to_ref_cell<dim>(Point<dim-1>(0.5),dof_opposite_faces[j]);
            Point<dim> x_n_real;
            
            double correction_norm = 1.0;
            int outside_check = 0;
            
            while (correction_norm > 1e-10 && outside_check<3)
            {
              const double perturbation = 0.01;
                      
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
              
              correction_norm = residual_n_vec.l2_norm();
              
              Tensor<1,dim> correction = transform_ref_face_correction_to_ref_cell<dim, Vector<double>>(residual_n_vec,dof_opposite_faces[j]);
              
              Point<dim> x_n_p1_ref = stencil_ref[0] + correction;
              
              double relaxation = 1.0;
              
              bool check = false;
              for (unsigned int k = 0; k < dim; ++k)
              {
                if (x_n_p1_ref[k] > 1.0 + 1e-12|| x_n_p1_ref[k] < 0.0 - 1e-12)
                {
                  check = true;
                  if (correction[k] > 1e-12)
                  {
                    relaxation = std::min((1.0 - x_n_ref[k])/(correction[k]+1e-12), relaxation);
                  }
                  else if (correction[k] < -1e-12)
                  {
                    relaxation = std::min((0.0 - x_n_ref[k])/(correction[k]+1e-12), relaxation);
                  }
              
                }
              }
              
              if (check)
                outside_check +=1;
              
              x_n_p1_ref = stencil_ref[0] + relaxation*correction;
              
              std::vector<Point<dim>> x_n_p1_ref_vec = {x_n_p1_ref}; 
              
              fe_point_evaluation.reinit(cell, x_n_p1_ref_vec);
              
              Point<dim> x_n_p1_real = fe_point_evaluation.quadrature_point(0);
              
              x_n_ref = x_n_p1_ref;
              x_n_real = x_n_p1_real;
          }
          
          const Tensor<1,dim> x_n_to_x_J_real = x_J_real - x_n_real;
          
          fe_point_evaluation.evaluate(cell_dof_values, EvaluationFlags::values);
          
          double distance_value_at_x_n = fe_point_evaluation.get_value(0);

          double approx_distance = compute_distance(x_n_to_x_J_real, distance_value_at_x_n);
          
          if (distance(dof_indices[i])> (approx_distance + 1e-8))
          {    
            change = true;
            distance(dof_indices[i]) = approx_distance;
          }
        }
        }
      }
    }
    
    // Exchange and "select" min value between the processes
    distance.compress(VectorOperation::min);

    // Update local ghost (distance becomes read only)
    distance.update_ghost_values();

    // Copy distance to distance_with_ghost
    distance_with_ghost = distance;

    // Zero out ghost DOFs to regain write functionalities in distance (it becomes write only, that is why we need distance_with_ghost - to read the values in it.)
    distance.zero_out_ghost_values();  
    
    // Copy the ghost values back in distance (zero_out_ghost_values() puts zeros in ghost DOFs)
    // distance = distance_with_ghost;
    for (auto p : this->locally_active_dofs)
    {
      distance(p) = distance_with_ghost(p);
    }
    
    change = Utilities::MPI::logical_or(change, mpi_communicator);
    
    pcout << "count = "<< count << std::endl;
    
  }
  
  constraints.distribute(distance);
  
  for (auto p : this->locally_active_dofs)
  {
    const double level_set_value = level_set(p);
    signed_distance(p) = distance(p)*sgn(signed_distance(p));
  }
  signed_distance.update_ghost_values();
  
}

template <int dim>
void AdvectionProblem<dim>::compute_volume()
{
  
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_quad_points = fe.degree + 1;
  
  const BoundingBox<dim> unit_box = create_unit_bounding_box<dim>();
  LocalCellWiseFunction<dim> level_set_function = LocalCellWiseFunction<dim>();
  
  hp::QCollection<1> q_collection;
  q_collection.push_back(QGauss<1>(n_quad_points));

  NonMatching::QuadratureGenerator<dim> quadrature_generator = NonMatching::QuadratureGenerator<dim>(q_collection);
  
  FEValues<dim> fe_values(fe,
                          QGauss<dim>(fe.degree + 1),
                          update_values | update_JxW_values);
  
  const unsigned int n_q_points    = fe_values.n_quadrature_points;
  std::vector<double> phase_values(n_q_points);
  
  double volume_sharp = 0.0;
  double volume_phase = 0.0;
  
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      
      Vector<double> cell_dof_values(dofs_per_cell);
      
      cell->get_dof_values(level_set, cell_dof_values.begin(), cell_dof_values.end());
      
      volume_sharp += compute_cell_wise_volume(cell, cell_dof_values, 0.0, unit_box, level_set_function, quadrature_generator);
      
      fe_values.reinit(cell);
      fe_values.get_function_values(locally_relevant_solution, phase_values);

      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          volume_phase += fe_values.JxW(q) * phase_values[q];
        }
    }
  }
  volume_sharp = Utilities::MPI::sum(volume_sharp, mpi_communicator);
  volume_phase = Utilities::MPI::sum(volume_phase, mpi_communicator);
  
  pcout << "volume_sharp = " << volume_sharp << std::endl;
  pcout << "volume_phase = " << volume_phase << std::endl;
  
  // 
  // double volume = 0.0;
  // FEValues<dim> fe_values(fe,
  //                         QGauss<dim>(fe.degree + 1),
  //                         update_values | update_JxW_values);
  // 
  // const unsigned int n_q_points    = fe_values.n_quadrature_points;
  // std::vector<double> phase_values(n_q_points);
  // 
  // for (const auto &cell : dof_handler.active_cell_iterators())
  //   {
  //     if (cell->is_locally_owned())
  //       {
  //         fe_values.reinit(cell);
  //         fe_values.get_function_values(locally_relevant_solution, phase_values);
  // 
  //         for (unsigned int q = 0; q < n_q_points; ++q)
  //           {
  //             volume += fe_values.JxW(q) * phase_values[q];
  //           }
  //       }
  //     }
  // 
  //     volume = Utilities::MPI::sum(volume, mpi_communicator);
      
      table_volume_monitoring.add_value("time", time);
      table_volume_monitoring.set_scientific("time", true);
      
      table_volume_monitoring.add_value("volume_sharp", volume_sharp);
      table_volume_monitoring.set_scientific("volume_sharp", true);
      
      table_volume_monitoring.add_value("volume_phase", volume_phase);
      table_volume_monitoring.set_scientific("volume_phase", true);
      // pcout << "Volume conservation" << std::endl;
      // table_volume_monitoring.write_text(std::cout);
              
      std::ofstream output("output/volume.dat");
      table_volume_monitoring.write_text(output);
      
}

template <int dim>
void AdvectionProblem<dim>::output_results(const int time_iteration) const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  
  data_out.add_data_vector(locally_relevant_solution, "locally_relevant_solution");
  data_out.add_data_vector(previous_solution, "previous_solution");
  
  
  data_out.add_data_vector(level_set, "level_set");
  
  signed_distance.update_ghost_values();
  data_out.add_data_vector(signed_distance, "signed_distance");
  signed_distance.zero_out_ghost_values();  
  
  distance.update_ghost_values();
  data_out.add_data_vector(distance, "distance");
  distance.zero_out_ghost_values();  
  
  data_out.add_data_vector(distance_with_ghost, "distance_with_ghost");

  volume_correction.update_ghost_values();
  data_out.add_data_vector(volume_correction, "volume_correction");
  volume_correction.zero_out_ghost_values();  
  
  data_out.build_patches();

  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::CompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);

  const std::string filename = "solution.vtu";

  data_out.write_vtu_with_pvtu_record("output/",
                                      "solution",
                                      time_iteration,
                                      MPI_COMM_WORLD,
                                      3);
                                      
                                        
}

template <int dim>
void AdvectionProblem<dim>::run()
{
  pcout << "Bonjour from run" << std::endl;
  
  Point<dim> p_0 = Point<dim>();
  p_0[0] = 0;
  for (unsigned int i = 1; i < dim; ++i)
    p_0[i] = 0;
    
  Point<dim> p_1 = Point<dim>();
  p_1[0] = 1;
  for (unsigned int i = 1; i < dim; ++i)
    p_1[i] = 1;
    
  std::vector< unsigned int > repetitions(dim);
  repetitions[0] = 1;
  for (unsigned int i = 1; i < dim; ++i)
    repetitions[i] = 1;
    
  GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, p_0, p_1);
  triangulation.refine_global(8);
            
  pcout << "Setup system" << std::endl;
  setup_system();
  set_initial_conditions();
  
  unsigned int it = 0;
  double final_time = 2.0;

  refine_grid(8,6);
  refine_grid(8,6);
  
  compute_level_set_from_phase_fraction();

  
  
  previous_solution = locally_relevant_solution;
  
  const unsigned int n_quad_points = fe.degree + 1;
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  
  const BoundingBox<dim> unit_box = create_unit_bounding_box<dim>();
  LocalCellWiseFunction<dim> level_set_function = LocalCellWiseFunction<dim>();
  
  hp::QCollection<1> q_collection;
  q_collection.push_back(QGauss<1>(n_quad_points));

  NonMatching::QuadratureGenerator<dim> quadrature_generator = NonMatching::QuadratureGenerator<dim>(q_collection);
  
  double global_volume = 0.0;
  
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      
      Vector<double> cell_dof_values(dofs_per_cell);
      
      cell->get_dof_values(level_set, cell_dof_values.begin(), cell_dof_values.end());
      
      global_volume += compute_cell_wise_volume(cell, cell_dof_values, 0.0, unit_box, level_set_function, quadrature_generator);
  
    }
  }
  global_volume = Utilities::MPI::sum(global_volume, mpi_communicator);
  
  pcout << "global_volume = " << global_volume << std::endl;
  
  mesh_classifier.reclassify();
  compute_sign_distance(it, global_volume);
  
  compute_phase_fraction_from_level_set();
  
  output_results(it);
  // compute_volume();
  
  
  pcout << "Solve system" << std::endl;
  while (time < final_time) 
  {
    it += 1;
    time += this->dt;
    
    assemble_system();
    
    pcout << "it = " << it << " time = " << time << std::endl;

    solve();
    
    compute_level_set_from_phase_fraction();
    
    if (it%1 == 0)
    {
      mesh_classifier.reclassify();
      compute_sign_distance(it, global_volume);
      
      compute_phase_fraction_from_level_set();
      
      compute_level_set_from_phase_fraction();
    }
    
    compute_volume();
    
    output_results(it);

    refine_grid(8,6);
    
    previous_solution = locally_relevant_solution;
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


