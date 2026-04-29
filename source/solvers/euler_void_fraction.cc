#include <solvers/euler_void_fraction.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

#include <limits>

namespace
{
  double
  global_min_owned(const TrilinosWrappers::MPI::Vector &vec,
                   const MPI_Comm                      &mpi_communicator)
  {
    double local_min = std::numeric_limits<double>::max();

    for (const auto i : vec.locally_owned_elements())
      local_min = std::min(local_min, vec[i]);

    return Utilities::MPI::min(local_min, mpi_communicator);
  }

  double
  global_max_owned(const TrilinosWrappers::MPI::Vector &vec,
                   const MPI_Comm                      &mpi_communicator)
  {
    double local_max = -std::numeric_limits<double>::max();

    for (const auto i : vec.locally_owned_elements())
      local_max = std::max(local_max, vec[i]);

    return Utilities::MPI::max(local_max, mpi_communicator);
  }
} // namespace

template <int dim>
EulerEulerVoidFraction<dim>::EulerEulerVoidFraction(
  const MPI_Comm &mpi_communicator,
  const bool      verbose)
  : mpi_communicator(mpi_communicator)
  , pcout(std::cout,
          Utilities::MPI::this_mpi_process(mpi_communicator) == 0 && verbose)
  , verbose(verbose)
{}

template <int dim>
void
EulerEulerVoidFraction<dim>::set_solid_volume_fraction(
  const TrilinosWrappers::MPI::Vector &alpha_s_in)
{
  alpha_s.reinit(alpha_s_in);
  alpha_s = alpha_s_in;

  alpha_f.reinit(alpha_s_in.locally_owned_elements(), mpi_communicator);

  has_alpha_s = true;
  has_alpha_f = false;
}

template <int dim>
void
EulerEulerVoidFraction<dim>::calculate_alpha_f()
{
  AssertThrow(has_alpha_s,
              ExcMessage("Solid volume fraction alpha_s is not set."));

  for (const auto &i : alpha_f.locally_owned_elements())
    {
      alpha_f[i] = std::max(1.0 - alpha_s[i], min_alpha_f_value);
    }

  alpha_f.compress(VectorOperation::insert);

  has_alpha_f = true;
}

template <int dim>
const TrilinosWrappers::MPI::Vector &
EulerEulerVoidFraction<dim>::get_alpha_f() const
{
  AssertThrow(has_alpha_f,
              ExcMessage("Fluid volume fraction alpha_f is not calculated."));
  return alpha_f;
}


template class EulerEulerVoidFraction<2>;
template class EulerEulerVoidFraction<3>;