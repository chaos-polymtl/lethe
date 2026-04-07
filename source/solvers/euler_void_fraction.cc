#include <solvers/euler_void_fraction.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

using namespace dealii;

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
}


template <int dim>
void
EulerEulerVoidFraction<dim>::pass_alpha_f_to_fluid()
{
  AssertThrow(has_alpha_f,
              ExcMessage(
                "alpha_f must be calculated before passing to fluid."));

  if (verbose)
    pcout << "Passing alpha_f to FluidDynamicsVANS" << std::endl;

  fluid_solver.set_alpha_f(alpha_f);
}


template class EulerEulerVoidFraction<2>;
template class EulerEulerVoidFraction<3>;