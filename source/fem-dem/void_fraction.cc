#include <fem-dem/void_fraction.h>

#include <deal.II/base/timer.h>

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/fe/fe_values.h>


using namespace dealii;

template <int dim>
void VoidFractionBase<dim>::setup_dofs()
{
  // Get a constant copy of the communicator since it is used extensively to establish the solution vectors
  const MPI_Comm mpi_communicator = dof_handler.get_communicator();

  dof_handler.distribute_dofs(fe);
  locally_owned_dofs = dof_handler.locally_owned_dofs();

  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  void_fraction_constraints.clear();
  void_fraction_constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler,
                                          void_fraction_constraints);

/// TODO ADD Periodic constraints
  void_fraction_constraints.close();

  void_fraction_locally_relevant.reinit(locally_owned_dofs,
                                      locally_relevant_dofs,
                                      mpi_communicator);

  // Initialize vector of previous solutions for the void fraction
  for (auto &solution : this->previous_void_fraction)
    {
      solution.reinit(this->locally_owned_dofs_voidfraction,
                      this->locally_relevant_dofs_voidfraction,
                      this->mpi_communicator);
    }

  void_fraction_locally_owned.reinit(locally_owned_dofs,
                                     this->mpi_communicator);


  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  void_fraction_constraints,
                                  false);
  SparsityTools::distribute_sparsity_pattern(
    dsp,
    locally_owned_dofs,
    this->mpi_communicator,
    locally_relevant_dofs);

  system_matrix_void_fraction.reinit(locally_owned_dofs,
                                     locally_owned_dofs,
                                     dsp,
                                     this->mpi_communicator);

  complete_system_matrix_void_fraction.reinit(locally_owned_dofs,
                                              locally_owned_dofs,
                                              dsp,
                                              this->mpi_communicator);

  system_rhs_void_fraction.reinit(locally_owned_dofs,
                                  this->mpi_communicator);

  active_set.set_size(dof_handler.n_dofs());

  mass_matrix.reinit(locally_owned_dofs,
                     locally_owned_dofs,
                     dsp,
                     this->mpi_communicator);

 assemble_mass_matrix_diagonal(mass_matrix);
}

template <int dim>
void VoidFractionBase<dim>::calculate_void_fraction(const double time)
{

  if (parameters.mode !=
      Parameters::VoidFractionMode::function)
    {
      // The void fraction is established using a particle handler.
      // A right-hand side and a linear system of equation are assembled and then solved. The resulting solution yields the nodal values of the void fraction.
      if (parameters.mode ==
          Parameters::VoidFractionMode::pcm)
        ;
      // particle_centered_method();
      else if (parameters.mode ==
               Parameters::VoidFractionMode::qcm)
        ;
      //  quadrature_centered_sphere_method(load_balance_step);
      else if (parameters.mode ==
               Parameters::VoidFractionMode::spm)
        ;
      //  satellite_point_method();

      // solve_L2_system_void_fraction();
      // if (this->cfd_dem_simulation_parameters.void_fraction
      //       ->bound_void_fraction == true)
      //   update_solution_and_constraints();

    }
  else
    {
      // The void fraction is established using a function. This is a straightforward usage of VectorTools.

      // The mapping is assumed to be of equal order to that of the interpolation of the void fraction for consistancy.
      const MappingQ<dim> mapping(fe.degree);

      // The current time of the function is set for time-dependant functions
      parameters.void_fraction.set_time(time);

      // The function is directly interpolate at the nodes. This is not an L2 projection, but a direct evaluation. This may lead to some issues on coarses meshes if a high-order interpolation (>FE_Q(1)) is used.
      VectorTools::interpolate(
        mapping,
        dof_handler,
        parameters.void_fraction,
        void_fraction_locally_owned);

      // Propagate ghost values
      void_fraction_locally_relevant = void_fraction_locally_owned;
    }
}

template <int dim>
void VoidFractionBase<dim>::assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &diagonal_mass_matrix)
{
  Assert(fe.degree == 1, ExcMessage("Constraining the void fraction between lower and upper bound is not supported when using a FE_Q with a degree higher than 1"));
  QGauss<dim>        quadrature_formula(this->number_quadrature_points);
  FEValues<dim>      fe_void_fraction_values(fe,
                                        quadrature_formula,
                                        update_values | update_JxW_values);
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_qpoints     = quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_void_fraction_values.reinit(cell);
          cell_matrix = 0;
          for (unsigned int q = 0; q < n_qpoints; ++q)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_matrix(i, i) += (fe_void_fraction_values.shape_value(i, q) *
                                    fe_void_fraction_values.shape_value(i, q) *
                                    fe_void_fraction_values.JxW(q));
          cell->get_dof_indices(local_dof_indices);
          void_fraction_constraints.distribute_local_to_global(
            cell_matrix, local_dof_indices, mass_matrix);
        }
    }
}