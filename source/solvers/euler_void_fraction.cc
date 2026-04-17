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
  FluidDynamicsVANS<dim> &fluid_solver,
  const MPI_Comm         &mpi_communicator,
  const bool              verbose)
  : fluid_solver(fluid_solver)
  , mpi_communicator(mpi_communicator)
  , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
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

  if (verbose)
    {
      const double alpha_s_min = global_min_owned(alpha_s, mpi_communicator);
      const double alpha_s_max = global_max_owned(alpha_s, mpi_communicator);

      pcout << "alpha_s range: [" << alpha_s_min << ", " << alpha_s_max << "]"
            << std::endl;
    }
}


template <int dim>
void
EulerEulerVoidFraction<dim>::calculate_alpha_f()
{
  AssertThrow(has_alpha_s,
              ExcMessage("alpha_s must be set before calculating alpha_f."));

  if (verbose)
    pcout << "Calculating alpha_f = 1 - alpha_s" << std::endl;

  alpha_f = 0.0;

  for (const auto i : alpha_f.locally_owned_elements())
    {
      const double value = 1.0 - alpha_s[i];
      alpha_f[i]         = std::max(min_alpha_f_value, std::min(1.0, value));
    }

  alpha_f.compress(VectorOperation::insert);
  has_alpha_f = true;

  if (verbose)
    {
      const double alpha_f_min = global_min_owned(alpha_f, mpi_communicator);
      const double alpha_f_max = global_max_owned(alpha_f, mpi_communicator);

      pcout << "alpha_f range: [" << alpha_f_min << ", " << alpha_f_max << "]"
            << std::endl;
    }
}


template <int dim>
void
EulerEulerVoidFraction<dim>::pass_alpha_f_to_fluid()
{
  AssertThrow(has_alpha_f,
              ExcMessage(
                "alpha_f must be calculated before passing to fluid."));

  if (verbose)
    pcout << "Passing alpha_f to Fluid" << std::endl;

  fluid_solver.set_alpha_f(alpha_f);
}


template class EulerEulerVoidFraction<2>;
template class EulerEulerVoidFraction<3>;