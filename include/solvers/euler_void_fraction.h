#ifndef euler_void_fraction_h
#define euler_void_fraction_h

#include <fem-dem/fluid_dynamics_vans.h>

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
using namespace dealii;

template <int dim>
class EulerEulerVoidFraction
{
public:
  EulerEulerVoidFraction(FluidDynamicsVANS<dim> &fluid_solver,
                         const MPI_Comm         &mpi_communicator,
                         const bool              verbose = true);

  void
  set_solid_volume_fraction(const TrilinosWrappers::MPI::Vector &alpha_s_in);

  void
  calculate_alpha_f();

  void
  pass_alpha_f_to_fluid();

private:
  FluidDynamicsVANS<dim> &fluid_solver;

  MPI_Comm           mpi_communicator;
  ConditionalOStream pcout;
  const bool         verbose;

  TrilinosWrappers::MPI::Vector alpha_s;
  TrilinosWrappers::MPI::Vector alpha_f;

  bool   has_alpha_s       = false;
  bool   has_alpha_f       = false;
  double min_alpha_f_value = 1e-12;
};

#endif