#ifndef euler_void_fraction_h
#define euler_void_fraction_h



#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_vector.h>
using namespace dealii;

template <int dim>
class EulerEulerVoidFraction
{
public:
  EulerEulerVoidFraction(const MPI_Comm &mpi_communicator,
                         const bool      verbose = true);

  void
  set_solid_volume_fraction(const TrilinosWrappers::MPI::Vector &alpha_s_in);

  void
  calculate_alpha_f();

  const TrilinosWrappers::MPI::Vector &
  get_alpha_f() const;

private:
  MPI_Comm           mpi_communicator;
  ConditionalOStream pcout;
  const bool         verbose;

  TrilinosWrappers::MPI::Vector alpha_s;
  TrilinosWrappers::MPI::Vector alpha_f;

  bool has_alpha_s = false;
  bool has_alpha_f = false;

  double min_alpha_f_value = 1e-12;
};

#endif