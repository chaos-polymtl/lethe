/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, University of Heidelberg, 2000
 */


// @sect3{Include files}

//BASE
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/table_handler.h>


//LAC
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_solver.h>
// Trilinos includes
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>


//GRID
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>


//DOFS
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

//FE
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

//Numerics
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

// Distributed
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>

#include "exactsolutions.h"
#include "forcingfunctions.h"
#include "boundaryconditions.h"
#include "parameters.h"
#include "simulationcontrol.h"
#include "navierstokessolverparameters.h"

#include <fstream>
#include <iostream>

// Finally, this is as in previous programs:
using namespace dealii;

#ifndef LETHE_GLSNS_H
#define LETHE_GLSNS_H

template <int dim>
class GLSNavierStokesSolver
{

public:
  GLSNavierStokesSolver(NavierStokesSolverParameters<dim> nsparam, const unsigned int degreeVelocity, const unsigned int degreePressure);

  ~GLSNavierStokesSolver();

  void refine_mesh();
  void setup_dofs();
  double calculateL2Error();
  void setInitialCondition(Parameters::InitialConditionType initial_condition_type );
  void postprocess();
  void finishTimeStep();
  void setSolutionVector(double value);

  void newton_iteration(const bool is_initial_step);
  void make_cube_grid(int refinementLevel);

  Function<dim> *exact_solution;
  Function<dim> *forcing_function;

protected:
  MPI_Comm                         mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;
  parallel::distributed::Triangulation<dim> triangulation;





private:
  template <Parameters::SimulationControl::TimeSteppingMethod scheme> void assembleGLS(const bool initial_step,
                                                                                       const bool assemble_matrix);

  void assemble_L2_projection();
  void set_nodal_values();
  void refine_mesh_Kelly();
  void refine_mesh_uniform();
  void projectRHSPressureConstant();
  void projectNewtonUpdatePressureConstant();
  void initializePressureRHSCorrection();
  void calculate_forces();
  void calculate_torques();
  void assemble_system(const bool initial_step);
  void assemble_rhs(const bool initial_step);
  void solve(bool initial_step, double relative_residual, double minimum_residual); // Interface function
  void solveGMRES(bool initial_step, double absolute_residual, double relative_residual);
  void solveBiCGStab(bool initial_step, double absolute_residual, double relative_residual);
  void solveAMG(bool initial_step, double absolute_residual, double relative_residual);
  void write_output_forces();
  void write_output_torques();
  void write_output_results(const std::string folder, const std::string solutionName, const unsigned int cycle, const double time);

  DoFHandler<dim>                  dof_handler;
  FESystem<dim>                    fe;

  IndexSet                         locally_owned_dofs;
  IndexSet                         locally_relevant_dofs;


  AffineConstraints<double>        zero_constraints;
  AffineConstraints<double>        nonzero_constraints;

  SparsityPattern                  sparsity_pattern;
  TrilinosWrappers::SparseMatrix   system_matrix;

  TrilinosWrappers::MPI::Vector    present_solution;
  TrilinosWrappers::MPI::Vector    newton_update;
  TrilinosWrappers::MPI::Vector    system_rhs;
  TrilinosWrappers::MPI::Vector    evaluation_point;
  TrilinosWrappers::MPI::Vector    local_evaluation_point;
  TrilinosWrappers::MPI::Vector    pressure_shape_function_integrals;
  TrilinosWrappers::MPI::Vector    pressure_projection;

  TrilinosWrappers::MPI::Vector    p1_solution;
  TrilinosWrappers::MPI::Vector    p2_solution;


  std::vector<types::global_dof_index> dofs_per_block;

  // Finite element order used
  const  unsigned int            degreeVelocity_;
  const  unsigned int            degreePressure_;
  const  unsigned int            degreeQuadrature_;


  double                         globalVolume_;
  const bool                     SUPG=true;
  const double                   GLS_u_scale=1;
  std::vector<std::pair<double,std::string> > times_and_names_;



protected:

  // Physical Properties
  double                                viscosity_;
  std::vector<double>                   L2ErrorU_;

  ConditionalOStream                    pcout;
  TimerOutput                           computing_timer;

  // Force analysis
  std::vector<Tensor<1,dim>>     forces_;
  std::vector<Tensor<1,3>>       torques_;
  std::vector<TableHandler>      forces_tables;
  std::vector<TableHandler>      torques_tables;


  NavierStokesSolverParameters<dim>     nsparam;

  SimulationControl                               simulationControl;
  Parameters::LinearSolver                        linearSolverParameters;
  Parameters::NonLinearSolver                     nonLinearSolverParameters;
  Parameters::MeshAdaptation                      meshAdaptationParameters;
  Parameters::Mesh                                meshParameters;
  Parameters::PhysicalProperties                  physicalProperties;
  Parameters::Timer                               clock;
  Parameters::FEM                                 femParameters;
  Parameters::Forces                              forcesParameters;
  Parameters::AnalyticalSolution                  analyticalSolutionParameters;
  BoundaryConditions::NSBoundaryConditions<dim>   boundaryConditions;
  Parameters::InitialConditions<dim>              *initialConditionParameters;
};

// Constructor
template<int dim>
GLSNavierStokesSolver<dim>::GLSNavierStokesSolver(NavierStokesSolverParameters<dim> p_nsparam, const unsigned int degreeVelocity, const unsigned int degreePressure):
    mpi_communicator (MPI_COMM_WORLD),
    n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
    this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening)),
    dof_handler (triangulation),
    fe(FE_Q<dim>(degreeVelocity), dim, FE_Q<dim>(degreePressure), 1),
    degreeVelocity_(degreeVelocity),degreePressure_(degreePressure),degreeQuadrature_(degreeVelocity+1),
    pcout (std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    computing_timer (mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times),
    nsparam(p_nsparam)
{
  boundaryConditions           = nsparam.boundaryConditions;
  initialConditionParameters   = nsparam.initialCondition;
  meshParameters               = nsparam.meshParameters;
  linearSolverParameters       = nsparam.linearSolver;
  nonLinearSolverParameters    = nsparam.nonLinearSolver;
  meshAdaptationParameters     = nsparam.meshAdaptation;
  physicalProperties           = nsparam.physicalProperties;
  clock                        = nsparam.clock;
  femParameters                = nsparam.femParameters;
  analyticalSolutionParameters = nsparam.analyticalSolution;
  forcesParameters             = nsparam.forcesParameters;
  simulationControl            = nsparam.simulationControl;

  // Change the behavior of the timer for situations when you don't want outputs
  if (clock.type==Parameters::Timer::none)
    computing_timer.disable_output();

  viscosity_=physicalProperties.viscosity;
  forces_.resize(boundaryConditions.size);
  torques_.resize(boundaryConditions.size);
  forces_tables.resize(boundaryConditions.size);
  torques_tables.resize(boundaryConditions.size);

  pcout << "Running on " << Utilities::MPI::n_mpi_processes(mpi_communicator)<< " MPI rank(s)..." << std::endl;
}

template <int dim>
GLSNavierStokesSolver<dim>::~GLSNavierStokesSolver ()
{
    dof_handler.clear ();
}



template <int dim>
void GLSNavierStokesSolver<dim>::make_cube_grid (int refinementLevel)
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (refinementLevel);
}

template <int dim>
void GLSNavierStokesSolver<dim>::finishTimeStep()
{
  if (simulationControl.getMethod()!=Parameters::SimulationControl::steady)
    p1_solution=present_solution;
  if (this->clock.type==Parameters::Timer::iteration)
  {
    this->computing_timer.print_summary ();
    this->computing_timer.reset ();
  }
}

template <int dim>
void GLSNavierStokesSolver<dim>::setSolutionVector(double value)
{
  present_solution=value;
}




template <int dim>
void GLSNavierStokesSolver<dim>::setup_dofs ()
{
  TimerOutput::Scope t(computing_timer, "setup_dofs");

  system_matrix.clear();

  dof_handler.distribute_dofs(fe);
//  DoFRenumbering::boost::Cuthill_McKee(dof_handler);
  DoFRenumbering::Cuthill_McKee(dof_handler);

  locally_owned_dofs = dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dof_handler,
                                           locally_relevant_dofs);


  std::vector<unsigned int> block_component(dim+1, 0);
  block_component[dim] = 1;
  dofs_per_block.resize (2);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);

  const MappingQ<dim>      mapping (degreeVelocity_,femParameters.qmapping_all);
  FEValuesExtractors::Vector velocities(0);

  {
    nonzero_constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);
    for (unsigned int i_bc=0 ; i_bc < boundaryConditions.size ; ++i_bc)
      {
        if(boundaryConditions.type[i_bc]==BoundaryConditions::noslip)
          {
                VectorTools::interpolate_boundary_values(mapping, dof_handler, i_bc, ZeroFunction<dim>(dim+1), nonzero_constraints,
                                                         fe.component_mask(velocities));
          }
        else if(boundaryConditions.type[i_bc]==BoundaryConditions::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert (i_bc);
            VectorTools::compute_no_normal_flux_constraints (dof_handler, 0,
                                                             no_normal_flux_boundaries,
                                                             nonzero_constraints
                                                             );
          }
        else if(boundaryConditions.type[i_bc]==BoundaryConditions::function)
          {
            VectorTools::interpolate_boundary_values(mapping,
                                                     dof_handler,
                                                     i_bc,
                                                     FunctionDefined<dim>(&boundaryConditions.bcFunctions[i_bc].u,
                                                                          &boundaryConditions.bcFunctions[i_bc].v,
                                                                          &boundaryConditions.bcFunctions[i_bc].w
                                                                          ),
                                                     nonzero_constraints,
                                                     fe.component_mask(velocities));
          }
      }
  }
  nonzero_constraints.close();

  {
    zero_constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);

    for (unsigned int i_bc=0 ; i_bc < boundaryConditions.size ; ++i_bc)
      {
        if(boundaryConditions.type[i_bc]==BoundaryConditions::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert (i_bc);
            VectorTools::compute_no_normal_flux_constraints (dof_handler, 0,
                                                             no_normal_flux_boundaries,
                                                             zero_constraints
                                                             );
          }
        else //if(boundaryConditions.boundaries[i_bc].type==Parameters::noslip || Parameters::function)
          {
            VectorTools::interpolate_boundary_values(mapping,
                                                     dof_handler,
                                                     i_bc,
                                                     ZeroFunction<dim>(dim+1),
                                                     zero_constraints,
                                                     fe.component_mask(velocities));
          }
      }
  }
  zero_constraints.close();


  present_solution.reinit (locally_owned_dofs,locally_relevant_dofs, mpi_communicator);
  p1_solution.reinit (locally_owned_dofs,locally_relevant_dofs, mpi_communicator);
  p2_solution.reinit (locally_owned_dofs,locally_relevant_dofs, mpi_communicator);

  newton_update.reinit (locally_owned_dofs, mpi_communicator);
  system_rhs.reinit (locally_owned_dofs, mpi_communicator);
  local_evaluation_point.reinit (locally_owned_dofs, mpi_communicator);

  DynamicSparsityPattern dsp (locally_relevant_dofs);
  DoFTools::make_sparsity_pattern (dof_handler, dsp, nonzero_constraints,false);
  SparsityTools::distribute_sparsity_pattern (dsp,
                                              dof_handler.n_locally_owned_dofs_per_processor(),
                                              mpi_communicator,
                                              locally_relevant_dofs);
  system_matrix.reinit (locally_owned_dofs,
                        locally_owned_dofs,
                        dsp,
                        mpi_communicator);

  globalVolume_=GridTools::volume(triangulation);


  pcout << "   Number of active cells:       " << triangulation.n_global_active_cells()  << std::endl
        << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
  pcout << "   Volume of triangulation:      " << globalVolume_ << std::endl;

  pressure_shape_function_integrals.reinit (locally_owned_dofs, mpi_communicator);
  pressure_projection.reinit (locally_owned_dofs, mpi_communicator);
  if (linearSolverParameters.rhsCorr)  initializePressureRHSCorrection();
}


template <int dim>
void
GLSNavierStokesSolver<dim>::initializePressureRHSCorrection()
{
  TimerOutput::Scope t(computing_timer, "pressure_RHS_projector");
  pressure_shape_function_integrals=0;
  QGauss<dim>   quadrature_formula(degreeQuadrature_);
  FEValues<dim> fe_values (fe,
                           quadrature_formula,
                           update_values |
                           update_quadrature_points |
                           update_JxW_values );

  const unsigned int                    dofs_per_cell = fe.dofs_per_cell;
  const unsigned int                    n_q_points    = quadrature_formula.size();
  const FEValuesExtractors::Scalar      pressure (dim);
  Vector<double>                        local_pressure_shape_function_integrals(dofs_per_cell);
  std::vector<types::global_dof_index>  local_dof_indices (dofs_per_cell);


  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          local_pressure_shape_function_integrals=0;
          for (unsigned int i=0; i<dofs_per_cell;++i)
            {
              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double phi_p = fe_values[pressure].value(i, q);
                  local_pressure_shape_function_integrals(i) += phi_p * fe_values.JxW(q);
                }
            }
          cell->get_dof_indices (local_dof_indices);
          zero_constraints.distribute_local_to_global(local_pressure_shape_function_integrals,
                                                      local_dof_indices,
                                                      pressure_shape_function_integrals);
        }
    }
  pressure_shape_function_integrals.compress (VectorOperation::add);
}

template <int dim>
template<Parameters::SimulationControl::TimeSteppingMethod scheme>
void GLSNavierStokesSolver<dim>::assembleGLS(const bool initial_step,
                                                const bool assemble_matrix)
{
  TimerOutput::Scope t(computing_timer, "assemble");

  if (assemble_matrix) system_matrix    = 0;
  system_rhs       = 0;
  QGauss<dim>   quadrature_formula(degreeQuadrature_);
  const MappingQ<dim>      mapping (degreeVelocity_,femParameters.qmapping_all);
  FEValues<dim> fe_values (mapping,
                           fe,
                           quadrature_formula,
                           update_values |
                           update_quadrature_points |
                           update_JxW_values |
                           update_gradients |
                           update_hessians);
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);
  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs    (dofs_per_cell);
  std::vector<Vector<double> >      rhs_force (n_q_points, Vector<double>(dim+1));
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  std::vector<Tensor<1, dim> >  present_velocity_values      (n_q_points);
  std::vector<Tensor<2, dim> >  present_velocity_gradients   (n_q_points);
  std::vector<double>           present_pressure_values      (n_q_points);
  std::vector<Tensor<1, dim> >  present_pressure_gradients   (n_q_points);
  std::vector<Tensor<1, dim> >  present_velocity_laplacians  (n_q_points);
  std::vector<Tensor<2, dim> >  present_velocity_hess        (n_q_points);

  Tensor<1, dim>  force;

  std::vector<double>           div_phi_u                 (dofs_per_cell);
  std::vector<Tensor<1, dim> >  phi_u                     (dofs_per_cell);
  std::vector<Tensor<3, dim> >  hess_phi_u                (dofs_per_cell);
  std::vector<Tensor<1, dim> >  laplacian_phi_u           (dofs_per_cell);
  std::vector<Tensor<2, dim> >  grad_phi_u                (dofs_per_cell);
  std::vector<double>           phi_p                     (dofs_per_cell);
  std::vector<Tensor<1, dim> >  grad_phi_p                (dofs_per_cell);

  // Time stepping. Will not be used if steady.
  const double dt=simulationControl.getTimeStep();
  const double sdt=1./dt;

  // Values at previous time step for backward Euler scheme
  std::vector<Tensor<1, dim> >  p1_velocity_values           (n_q_points);

  // Element size
  double h ;

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
  {
      if (cell->is_locally_owned())
      {
          if (dim==2) h = std::sqrt(4.* cell->measure() / M_PI) / degreeVelocity_;
          else if (dim==3) h = pow(6*cell->measure()/M_PI,1./3.) / degreeVelocity_;

          fe_values.reinit(cell);
          local_matrix = 0;
          local_rhs    = 0;
          fe_values[velocities].get_function_values(evaluation_point, present_velocity_values);
          fe_values[velocities].get_function_gradients(evaluation_point, present_velocity_gradients);
          fe_values[pressure].get_function_values(evaluation_point,present_pressure_values);
          fe_values[pressure].get_function_gradients(evaluation_point,present_pressure_gradients);
          fe_values[velocities].get_function_laplacians(evaluation_point,present_velocity_laplacians);
          forcing_function->vector_value_list(fe_values.get_quadrature_points(), rhs_force);

          if (scheme != Parameters::SimulationControl::steady)
            fe_values[velocities].get_function_values(p1_solution,p1_velocity_values);


          for (unsigned int q=0; q<n_q_points; ++q)
          {
              const double u_mag= std::max(present_velocity_values[q].norm(),1e-3*GLS_u_scale);
               double tau = 1./ std::sqrt(std::pow(2.*u_mag/h,2)+9*std::pow(4*viscosity_/(h*h),2));
//              if (scheme == Parameters::SimulationControl::backward)
//                tau = 1./ std::sqrt(std::pow(2*sdt,2)+std::pow(2.*u_mag/h,2)+9*std::pow(4*viscosity_/(h*h),2));

              for (unsigned int k=0; k<dofs_per_cell; ++k)
              {
                  div_phi_u[k]  =  fe_values[velocities].divergence (k, q);
                  grad_phi_u[k] =  fe_values[velocities].gradient(k, q);
                  phi_u[k]      =  fe_values[velocities].value(k, q);
                  hess_phi_u[k] =  fe_values[velocities].hessian(k, q);
                  phi_p[k]      =  fe_values[pressure]  .value(k, q);
                  grad_phi_p[k] =  fe_values[pressure]  .gradient(k, q);

                  for( int d=0; d<dim; ++d )
                      laplacian_phi_u[k][d] = trace( hess_phi_u[k][d] );
              }

              // Establish the force vector
              for( int i=0; i<dim; ++i )
              {
                  const unsigned int component_i = fe.system_to_component_index(i).first;
                  force[i] = rhs_force[q](component_i);
              }

              for (unsigned int i=0; i<dofs_per_cell; ++i)
              {

                  if (assemble_matrix)
                  {
                      for (unsigned int j=0; j<dofs_per_cell; ++j)
                      {
                          local_matrix(i, j) += (  viscosity_*scalar_product(grad_phi_u[j], grad_phi_u[i])
                                                   + present_velocity_gradients[q]*phi_u[j]*phi_u[i]
                                                   + grad_phi_u[j]*present_velocity_values[q]*phi_u[i]
                                                   - div_phi_u[i]*phi_p[j]
                                                   + phi_p[i]*div_phi_u[j]
                                                   )
                                  * fe_values.JxW(q);

                          // Mass matrix
                          if (scheme == Parameters::SimulationControl::backward)
                              local_matrix(i, j) += phi_u[j]*phi_u[i] *sdt * fe_values.JxW(q);

                          //PSPG GLS term
                          local_matrix(i, j) += tau*
                              (  present_velocity_gradients[q]*phi_u[j]*grad_phi_p[i]
                                 + grad_phi_u[j]*present_velocity_values[q]*grad_phi_p[i]
                                 + grad_phi_p[j]*grad_phi_p[i]
                                 - viscosity_* laplacian_phi_u[j] * grad_phi_p[i]
                                 )
                                  * fe_values.JxW(q);

                          if (scheme == Parameters::SimulationControl::backward)
                              local_matrix(i, j) += -tau * sdt * phi_u[j] * grad_phi_p[i] *  fe_values.JxW(q);

                          // Jacobian is currently incomplete
                          if (SUPG)
                          {
                              local_matrix(i, j) +=
                                     tau*
                                      (  present_velocity_gradients[q]*phi_u[j]*(present_velocity_values[q]* grad_phi_u[i])
                                         + grad_phi_u[j]*present_velocity_values[q]*(present_velocity_values[q]* grad_phi_u[i])
                                         + grad_phi_p[j]*(present_velocity_values[q]* grad_phi_u[i])
                                         - viscosity_* laplacian_phi_u[j] * (present_velocity_values[q]* grad_phi_u[i])
                                         )
                                      * fe_values.JxW(q)
                                      + tau*
                                      (    present_velocity_gradients[q]*present_velocity_values[q]*(phi_u[j]*grad_phi_u[i])
                                           + present_pressure_gradients[q]*(phi_u[j]*grad_phi_u[i])
                                           - viscosity_* present_velocity_laplacians[q] * (phi_u[j]*grad_phi_u[i])
                                           - force * (phi_u[j]*grad_phi_u[i])
                                           )
                                      * fe_values.JxW(q)
                                      ;

                              if (scheme == Parameters::SimulationControl::backward)
                                  local_matrix(i, j) +=
                                      -tau*sdt*phi_u[j]*(present_velocity_values[q]* grad_phi_u[i])* fe_values.JxW(q)
                                      - tau*sdt*(present_velocity_values[q]-p1_velocity_values[q])*(phi_u[j]*grad_phi_u[i])* fe_values.JxW(q);

                          }
                      }
                  }
//                  const unsigned int component_i = fe.system_to_component_index(i).first;

                  double present_velocity_divergence =  trace(present_velocity_gradients[q]);
                  local_rhs(i) += ( - viscosity_*scalar_product(present_velocity_gradients[q],grad_phi_u[i])
                                    - present_velocity_gradients[q]*present_velocity_values[q]*phi_u[i]
                                    + present_pressure_values[q]*div_phi_u[i]
                                    - present_velocity_divergence*phi_p[i]
                                    + force * phi_u[i]
                                    )
                          * fe_values.JxW(q);

                  if (scheme == Parameters::SimulationControl::backward)
                   local_rhs(i) += -sdt* (present_velocity_values[q]-p1_velocity_values[q])*phi_u[i]* fe_values.JxW(q);


                  // PSPG GLS term
                  local_rhs(i) +=  tau*
                          (
                              - present_velocity_gradients[q]*present_velocity_values[q]* grad_phi_p[i]
                              - present_pressure_gradients[q]*grad_phi_p[i]
                              + viscosity_* present_velocity_laplacians[q] * grad_phi_p[i]
                              + force * grad_phi_p[i]
                              )
                          * fe_values.JxW(q);

                  if (scheme == Parameters::SimulationControl::backward)
                   local_rhs(i) += tau * sdt* (present_velocity_values[q]-p1_velocity_values[q])*grad_phi_p[i]* fe_values.JxW(q);

                  ////             SUPG GLS term
                  if (SUPG)
                  {
                      local_rhs(i) += tau*
                              (
                                  - present_velocity_gradients[q] *present_velocity_values[q]* (present_velocity_values[q]* grad_phi_u[i])
                                  - present_pressure_gradients[q]* (present_velocity_values[q]* grad_phi_u[i])
                                  + viscosity_* present_velocity_laplacians[q] * (present_velocity_values[q]* grad_phi_u[i])
                                  + force *(present_velocity_values[q]* grad_phi_u[i])
                                  )* fe_values.JxW(q);

                      if (scheme == Parameters::SimulationControl::backward)
                       local_rhs(i) += tau * sdt* (present_velocity_values[q]-p1_velocity_values[q])* (present_velocity_values[q]* grad_phi_u[i])* fe_values.JxW(q);
                  }
                  //local_rhs(i) += fe_values.shape_value(i,q) *
                  //        rhs_force[q](component_i) *
                  //        fe_values.JxW(q);

              }
          }


          cell->get_dof_indices (local_dof_indices);
          const AffineConstraints<double> &constraints_used = initial_step ? nonzero_constraints : zero_constraints;
          if (assemble_matrix)
          {
              constraints_used.distribute_local_to_global(local_matrix,
                                                          local_rhs,
                                                          local_dof_indices,
                                                          system_matrix,
                                                          system_rhs);
          }
          else
          {
              constraints_used.distribute_local_to_global(local_rhs,
                                                          local_dof_indices,
                                                          system_rhs);
          }

      }
  }
  system_matrix.compress (VectorOperation::add);
  system_rhs.compress (VectorOperation::add);
}

/**
 * Set the initial condition using a L2 or a viscous solver
**/
template <int dim>
void GLSNavierStokesSolver<dim>::setInitialCondition (Parameters::InitialConditionType initial_condition_type)
{
  if (initial_condition_type == Parameters::InitialConditionType::L2projection)
  {
    assemble_L2_projection();
    solve(true,1e-12,1e-12);
    present_solution=newton_update;
    finishTimeStep();
    postprocess();
  }
  else if (initial_condition_type == Parameters::InitialConditionType::nodal)
  {
    set_nodal_values();
    finishTimeStep();
    postprocess();
  }

  else if (initial_condition_type == Parameters::InitialConditionType::viscous)
  {
    viscosity_=1;
    Parameters::SimulationControl::TimeSteppingMethod previousControl =  simulationControl.getMethod();
    simulationControl.setMethod(Parameters::SimulationControl::steady);
    newton_iteration(true);
    simulationControl.setMethod(previousControl);
    finishTimeStep();
    postprocess();
    simulationControl.setMethod(previousControl);
    viscosity_=this->physicalProperties.viscosity;
  }
  else
  {
    finishTimeStep();
    postprocess();
  }
}

template <int dim>
void GLSNavierStokesSolver<dim>::assemble_L2_projection()
{
  system_matrix    = 0;
  system_rhs       = 0;
  QGauss<dim>              quadrature_formula(degreeQuadrature_);
  const MappingQ<dim>      mapping (degreeVelocity_,femParameters.qmapping_all);
  FEValues<dim> fe_values (mapping,
                           fe,
                           quadrature_formula,
                           update_values |
                           update_quadrature_points |
                           update_JxW_values
                           );
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs    (dofs_per_cell);
  std::vector<Vector<double> >      initial_velocity (n_q_points, Vector<double>(dim+1));
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);


  Tensor<1, dim>  rhs_initial_velocity_pressure;
  double  rhs_initial_pressure;


  std::vector<Tensor<1, dim> >  phi_u                     (dofs_per_cell);
  std::vector<double>           phi_p                     (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
  {
    if (cell->is_locally_owned())
    {
      fe_values.reinit(cell);
      local_matrix = 0;
      local_rhs    = 0;
      initialConditionParameters->uvwp.vector_value_list(fe_values.get_quadrature_points(), initial_velocity);
      for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
        {
          phi_p[k]  =  fe_values[pressure]  .value(k, q);
          phi_u[k]  =  fe_values[velocities].value(k, q);
        }

        // Establish the rhs tensor operator
        for( int i=0; i<dim; ++i )
        {
            const unsigned int component_i = fe.system_to_component_index(i).first;
            rhs_initial_velocity_pressure[i] = initial_velocity[q](component_i);
        }
        rhs_initial_pressure=initial_velocity[q](dim);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          // Matrix assembly
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            local_matrix(i, j) += (phi_u[j]*phi_u[i]) * fe_values.JxW(q) ;
            local_matrix(i, j) += (phi_p[j]*phi_p[i]) * fe_values.JxW(q) ;
          }
          local_rhs(i) +=  (phi_u[i] * rhs_initial_velocity_pressure + phi_p[i]*rhs_initial_pressure) * fe_values.JxW(q);
        }
      }


      cell->get_dof_indices (local_dof_indices);
      const AffineConstraints<double> &constraints_used = nonzero_constraints;
      constraints_used.distribute_local_to_global(local_matrix,
                                                  local_rhs,
                                                  local_dof_indices,
                                                  system_matrix,
                                                  system_rhs);
    }
  }
  system_matrix.compress (VectorOperation::add);
  system_rhs.compress (VectorOperation::add);
}

template <int dim>
void GLSNavierStokesSolver<dim>::set_nodal_values()
{
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);
  const MappingQ<dim>      mapping (degreeVelocity_,femParameters.qmapping_all);
  VectorTools::interpolate(mapping,
                           dof_handler,
                           initialConditionParameters->uvwp,
                           newton_update,
                           fe.component_mask(velocities));
  VectorTools::interpolate(mapping,
                           dof_handler,
                           initialConditionParameters->uvwp,
                           newton_update,
                           fe.component_mask(pressure));
  present_solution=newton_update;
}

template <int dim>
void GLSNavierStokesSolver<dim>::assemble_system(const bool initial_step)
{
  if (simulationControl.getMethod()==Parameters::SimulationControl::backward)
     assembleGLS<Parameters::SimulationControl::backward>(initial_step,true);
  else if (simulationControl.getMethod()==Parameters::SimulationControl::steady)
      assembleGLS<Parameters::SimulationControl::steady>(initial_step,true);

}
template <int dim>
void GLSNavierStokesSolver<dim>::assemble_rhs(const bool initial_step)
{
  if (simulationControl.getMethod()==Parameters::SimulationControl::backward)
    assembleGLS<Parameters::SimulationControl::backward>(initial_step,false);
  else if (simulationControl.getMethod()==Parameters::SimulationControl::steady)
    assembleGLS<Parameters::SimulationControl::steady>(initial_step,false);
}

template <int dim>
void GLSNavierStokesSolver<dim>::projectRHSPressureConstant()
{
  // calculate projection of system_rhs on pressure_shape_function_integrals
  for (unsigned int dof_i=0 ; dof_i < pressure_shape_function_integrals.size() ; ++dof_i)
    pressure_projection[dof_i] = pressure_shape_function_integrals[dof_i]*system_rhs[dof_i];

  // Add it to RHS
  system_rhs.add(-1.,pressure_projection);
}

template <int dim>
void GLSNavierStokesSolver<dim>::projectNewtonUpdatePressureConstant()
{
  // calculate projection of system_rhs on pressure_shape_function_integrals
  for (unsigned int dof_i=0 ; dof_i < pressure_shape_function_integrals.size() ; ++dof_i)
    pressure_projection[dof_i] = pressure_shape_function_integrals[dof_i]*newton_update[dof_i];

  // Add it to RHS
  newton_update.add(-1.,pressure_projection);
}

template <int dim>
void GLSNavierStokesSolver<dim>::solve(const bool initial_step, double relative_residual, double minimum_residual)
{
  if (linearSolverParameters.solver==linearSolverParameters.gmres)         solveGMRES(initial_step,minimum_residual,relative_residual);
  else if (linearSolverParameters.solver==linearSolverParameters.bicgstab) solveBiCGStab(initial_step,minimum_residual,relative_residual);
  else if (linearSolverParameters.solver==linearSolverParameters.amg)      solveAMG(initial_step,minimum_residual,relative_residual);
  else throw("This solver is not allowed");
}


template <int dim>
void GLSNavierStokesSolver<dim>::solveGMRES (const bool initial_step, double absolute_residual, double relative_residual)
{
  TimerOutput::Scope t(computing_timer, "solve");
  if (linearSolverParameters.rhsCorr) projectRHSPressureConstant();
  const AffineConstraints<double> &constraints_used = initial_step ? nonzero_constraints : zero_constraints;
  const double linear_solver_tolerance = std::max(relative_residual*system_rhs.l2_norm(),absolute_residual);

  if (linearSolverParameters.verbosity!=Parameters::LinearSolver::quiet)
  {
    pcout << "  -Tolerance of iterative solver is : " << std::setprecision(linearSolverParameters.residual_precision) << linear_solver_tolerance << std::endl;
  }
  TrilinosWrappers::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);

  SolverControl solver_control (linearSolverParameters.max_iterations, linear_solver_tolerance,true,true);
  TrilinosWrappers::SolverGMRES solver(solver_control);

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill=linearSolverParameters.ilu_fill;
  const double ilu_atol=linearSolverParameters.ilu_atol ;
  const double ilu_rtol=linearSolverParameters.ilu_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(ilu_fill,ilu_atol,ilu_rtol,0);
  TrilinosWrappers::PreconditionILU preconditioner;

  preconditioner.initialize(system_matrix,preconditionerOptions);

  solver.solve (system_matrix,
               completely_distributed_solution,
               system_rhs,
               preconditioner
               );

  if (linearSolverParameters.verbosity!=Parameters::LinearSolver::quiet)
  {
    pcout << "  -Iterative solver took : " << solver_control.last_step() << " steps " << std::endl;
  }

  constraints_used.distribute(completely_distributed_solution);
  newton_update = completely_distributed_solution;
  if (linearSolverParameters.rhsCorr) projectNewtonUpdatePressureConstant();
}

template <int dim>
void GLSNavierStokesSolver<dim>::solveBiCGStab(const bool initial_step, double absolute_residual, double relative_residual)
{
  TimerOutput::Scope t(computing_timer, "solve");
  if (linearSolverParameters.rhsCorr) projectRHSPressureConstant();

  const AffineConstraints<double> &constraints_used = initial_step ? nonzero_constraints : zero_constraints;
  const double linear_solver_tolerance = std::max(relative_residual*system_rhs.l2_norm(),absolute_residual);
  if (linearSolverParameters.verbosity!=Parameters::LinearSolver::quiet)
  {
    pcout << "  -Tolerance of iterative solver is : " << std::setprecision(linearSolverParameters.residual_precision) << linear_solver_tolerance << std::endl;
  }
 TrilinosWrappers::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);

  SolverControl solver_control (linearSolverParameters.max_iterations, linear_solver_tolerance,true,true);
  TrilinosWrappers::SolverBicgstab solver(solver_control);

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill=linearSolverParameters.ilu_fill;
  const double ilu_atol=linearSolverParameters.ilu_atol ;
  const double ilu_rtol=linearSolverParameters.ilu_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(ilu_fill,ilu_atol,ilu_rtol,0);
  TrilinosWrappers::PreconditionILU preconditioner;

  preconditioner.initialize(system_matrix,preconditionerOptions);

  solver.solve (system_matrix,
               completely_distributed_solution,
               system_rhs,
               preconditioner
               );

  if (linearSolverParameters.verbosity!=Parameters::LinearSolver::quiet)
  {
    pcout << "  -Iterative solver took : " << solver_control.last_step() << " steps " << std::endl;
  }
  constraints_used.distribute(completely_distributed_solution);
  newton_update = completely_distributed_solution;
  if (linearSolverParameters.rhsCorr) projectNewtonUpdatePressureConstant();
}

template <int dim>
void GLSNavierStokesSolver<dim>::solveAMG (const bool initial_step, double absolute_residual, double relative_residual)
{
  TimerOutput::Scope t(computing_timer, "solve");

  const AffineConstraints<double> &constraints_used = initial_step ? nonzero_constraints : zero_constraints;

  const double linear_solver_tolerance = std::max(relative_residual*system_rhs.l2_norm(),absolute_residual);
  if (linearSolverParameters.verbosity!=Parameters::LinearSolver::quiet)
  {
    pcout << "  -Tolerance of iterative solver is : " << std::setprecision(linearSolverParameters.residual_precision) << linear_solver_tolerance << std::endl;
  }
  TrilinosWrappers::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);

  SolverControl solver_control (linearSolverParameters.max_iterations, linear_solver_tolerance,true,true);
  TrilinosWrappers::SolverGMRES solver(solver_control);

  TrilinosWrappers::PreconditionAMG preconditioner;

  std::vector<std::vector<bool> > constant_modes;
  // Constant modes include pressure since everything is in the same matrix
  std::vector<bool>  velocity_components (dim+1,true);
  velocity_components[dim] = true;
  DoFTools::extract_constant_modes (dof_handler, velocity_components,
                                    constant_modes);

  TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
  amg_data.constant_modes = constant_modes;

  const bool            elliptic=false;
  bool            higher_order_elements = false;
  if (degreeVelocity_>1) higher_order_elements = true;
  const unsigned int  	n_cycles = linearSolverParameters.amg_n_cycles;
  const bool            w_cycle = linearSolverParameters.amg_w_cycles;
  const double  	aggregation_threshold = linearSolverParameters.amg_aggregation_threshold;
  const unsigned int  	smoother_sweeps = linearSolverParameters.amg_smoother_sweeps;
  const unsigned int  	smoother_overlap = linearSolverParameters.amg_smoother_overlap;
  const bool            output_details = false;
  const char *  	smoother_type = "ILU";
  const char *  	coarse_type  =  "ILU";
  TrilinosWrappers::PreconditionAMG::AdditionalData preconditionerOptions(
        elliptic,
        higher_order_elements,
        n_cycles,
        w_cycle,
        aggregation_threshold,
        constant_modes,
        smoother_sweeps,
        smoother_overlap,
        output_details,
        smoother_type,
        coarse_type
        );

  Teuchos::ParameterList            parameter_ml;
  std::unique_ptr< Epetra_MultiVector > distributed_constant_modes;
  preconditionerOptions.set_parameters(parameter_ml, distributed_constant_modes, system_matrix);
  const double ilu_fill=linearSolverParameters.ilu_fill;
  const double ilu_atol=linearSolverParameters.ilu_atol ;
  const double ilu_rtol=linearSolverParameters.ilu_rtol;
  parameter_ml.set("smoother: ifpack level-of-fill",ilu_fill);
  parameter_ml.set("smoother: ifpack absolute threshold",ilu_atol);
  parameter_ml.set("smoother: ifpack relative threshold",ilu_rtol);

  parameter_ml.set("coarse: ifpack level-of-fill",ilu_fill);
  parameter_ml.set("coarse: ifpack absolute threshold",ilu_atol);
  parameter_ml.set("coarse: ifpack relative threshold",ilu_rtol);
  preconditioner.initialize(system_matrix,parameter_ml);

  solver.solve (system_matrix,
               completely_distributed_solution,
               system_rhs,
               preconditioner
               );

  if (linearSolverParameters.verbosity!=Parameters::LinearSolver::quiet)
  {
    pcout << "  -Iterative solver took : " << solver_control.last_step() << " steps " << std::endl;
  }

  constraints_used.distribute(completely_distributed_solution);

  newton_update = completely_distributed_solution;
}


template <int dim>
void GLSNavierStokesSolver<dim>::refine_mesh ()
{
  if ( simulationControl.getIter()%meshAdaptationParameters.frequency==0)
  {
    if (meshAdaptationParameters.type==meshAdaptationParameters.kelly) refine_mesh_Kelly();
    if (meshAdaptationParameters.type==meshAdaptationParameters.uniform) refine_mesh_uniform();
  }
}

template <int dim>
void GLSNavierStokesSolver<dim>::refine_mesh_Kelly ()
{
  // Time monitoring
  TimerOutput::Scope t(computing_timer, "refine");

  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  const MappingQ<dim>      mapping (degreeVelocity_,femParameters.qmapping_all);
  const FEValuesExtractors::Vector velocity(0);
  const FEValuesExtractors::Scalar pressure(dim);
  if (meshAdaptationParameters.variable==Parameters::MeshAdaptation::pressure)
  {
    KellyErrorEstimator<dim>::estimate (mapping,
                                        dof_handler,
                                        QGauss<dim-1>(degreeQuadrature_+1),
                                        typename std::map<types::boundary_id, const Function<dim, double> *>(),
                                        present_solution,
                                        estimated_error_per_cell,
                                        fe.component_mask(pressure));
  }
  if (meshAdaptationParameters.variable==Parameters::MeshAdaptation::velocity)
  {
    KellyErrorEstimator<dim>::estimate (mapping,
                                        dof_handler,
                                        QGauss<dim-1>(degreeQuadrature_+1),
                                        typename std::map<types::boundary_id, const Function<dim, double> *>(),
                                        present_solution,
                                        estimated_error_per_cell,
                                        fe.component_mask(velocity));
  }
  if (meshAdaptationParameters.fractionType==Parameters::MeshAdaptation::number)
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(triangulation,estimated_error_per_cell,
                                                    meshAdaptationParameters.fractionRefinement, meshAdaptationParameters.fractionCoarsening,
                                                    meshAdaptationParameters.maxNbElements);

  if (meshAdaptationParameters.fractionType==Parameters::MeshAdaptation::fraction)
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,estimated_error_per_cell,
                                                    meshAdaptationParameters.fractionRefinement, meshAdaptationParameters.fractionCoarsening);


  if (triangulation.n_levels() > meshAdaptationParameters.maxRefLevel)
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active(meshAdaptationParameters.maxRefLevel);
         cell != triangulation.end(); ++cell)
      cell->clear_refine_flag ();
  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active(meshAdaptationParameters.minRefLevel);
       cell != triangulation.end_active(meshAdaptationParameters.minRefLevel); ++cell)
    cell->clear_coarsen_flag ();


  triangulation.prepare_coarsening_and_refinement();
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector > solution_transfer(dof_handler);
  solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
  triangulation.execute_coarsening_and_refinement ();

  setup_dofs();

  // Set up the vectors for the transfer
  TrilinosWrappers::MPI::Vector tmp (locally_owned_dofs, mpi_communicator);

  // Interpolate the solution at time and previous time
  solution_transfer.interpolate( tmp);

  // Distribute constraints
  nonzero_constraints.distribute(tmp);

  // Fix on the new mesh
  present_solution = tmp;
}

template <int dim>
void GLSNavierStokesSolver<dim>::refine_mesh_uniform ()
{
    TimerOutput::Scope t(computing_timer, "refine");
    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector > solution_transfer(dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
    triangulation.refine_global(1);
    setup_dofs();
    TrilinosWrappers::MPI::Vector tmp (locally_owned_dofs, mpi_communicator);
    solution_transfer.interpolate(tmp);
    nonzero_constraints.distribute(tmp);
    present_solution = tmp;
}


template <int dim>
void GLSNavierStokesSolver<dim>::newton_iteration( const bool  is_initial_step)
{
  double current_res;
  double last_res;
  bool   first_step = is_initial_step;
    {
      unsigned int outer_iteration = 0;
      last_res = 1.0;
      current_res = 1.0;
      while ((first_step || (current_res > nonLinearSolverParameters.tolerance)) && outer_iteration < nonLinearSolverParameters.maxIterations)
        {

          if (first_step)
            {
              evaluation_point = present_solution; //local_evaluation_point;
              assemble_system(first_step);
              current_res = system_rhs.l2_norm();
              if (nonLinearSolverParameters.verbosity != nonLinearSolverParameters.quiet)
                pcout  << "Initial Newton iteration: " << outer_iteration << "  - Residual:  " << std::setprecision(nonLinearSolverParameters.display_precision) << current_res << std::endl;
              solve(first_step,linearSolverParameters.relative_residual,linearSolverParameters.minimum_residual);
              local_evaluation_point = newton_update;
              nonzero_constraints.distribute(local_evaluation_point);
              evaluation_point = local_evaluation_point;
              present_solution = evaluation_point;
              first_step = false;
              assemble_rhs(first_step);
              current_res = system_rhs.l2_norm();
              last_res = current_res;
            }
          else
            {
              evaluation_point = present_solution;
              assemble_system(first_step);
              if (outer_iteration==0)
              {
                current_res = system_rhs.l2_norm();
                last_res = current_res;
              }
              if (nonLinearSolverParameters.verbosity != nonLinearSolverParameters.quiet)
                pcout  << "Newton iteration: " << outer_iteration << "  - Residual:  " << current_res << std::endl;
              solve(first_step,linearSolverParameters.relative_residual,linearSolverParameters.minimum_residual);
              for (double alpha = 1.0; alpha > 1e-3; alpha *= 0.5)
                {
                  local_evaluation_point = present_solution;
                  local_evaluation_point.add(alpha, newton_update);
                  nonzero_constraints.distribute(local_evaluation_point);
                  evaluation_point=local_evaluation_point;
                  assemble_rhs(first_step);
                  current_res = system_rhs.l2_norm();
                  if (nonLinearSolverParameters.verbosity != nonLinearSolverParameters.quiet)
                    pcout << "\t\talpha = " << std::setw(6) << alpha << std::setw(0)
                            << " res = " << std::setprecision(nonLinearSolverParameters.display_precision) << current_res << std::endl;
                  if (current_res < 0.9*last_res || last_res < nonLinearSolverParameters.tolerance)
                    break;
                }
              {
                present_solution = evaluation_point;
                last_res = current_res;
              }
            }
          ++outer_iteration;

        }
    }
}

template <int dim>
void GLSNavierStokesSolver<dim>::postprocess()
{
    if (simulationControl.isOutputIteration()) write_output_results(simulationControl.getOutputFolder(),simulationControl.getOuputName(), simulationControl.getIter(), simulationControl.getTime());

    if(forcesParameters.calculate_force)
    {
      if (simulationControl.getIter()%forcesParameters.calculation_frequency==0)
        calculate_forces();
      if (simulationControl.getIter()%forcesParameters.output_frequency==0)
        write_output_forces();
    }

    if(forcesParameters.calculate_torque)
    {
      if (simulationControl.getIter()%forcesParameters.calculation_frequency==0)
        calculate_torques();
      if (simulationControl.getIter()%forcesParameters.output_frequency==0)
        write_output_torques();
    }
}

template <int dim>
void GLSNavierStokesSolver<dim>::write_output_results (const std::string folder, const std::string solutionName, const unsigned int iter, const double time)
{
  TimerOutput::Scope t(computing_timer, "output");
  const MappingQ<dim>      mapping (degreeVelocity_,femParameters.qmapping_all);

  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("pressure");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (present_solution, solution_names, DataOut<dim>::type_dof_data, data_component_interpretation);

  Vector<float> subdomain (triangulation.n_active_cells());
  for (unsigned int i=0; i<subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector (subdomain, "subdomain");


  data_out.build_patches (mapping,simulationControl.getSubdivision());

  const std::string filename = (folder+ solutionName + "." +
                                Utilities::int_to_string (iter, 4) +
                                "." +
                                Utilities::int_to_string
                                (triangulation.locally_owned_subdomain(), 4));
  std::ofstream output ((filename + ".vtu").c_str());
  data_out.write_vtu (output);

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  {
    std::vector<std::string> filenames;
    for (unsigned int i=0;
         i<Utilities::MPI::n_mpi_processes(mpi_communicator);
         ++i)
      filenames.push_back (solutionName + "." +
                           Utilities::int_to_string (iter, 4) +
                           "." +
                           Utilities::int_to_string (i, 4) +
                           ".vtu");

    std::string pvtu_filename = (solutionName + "."+
                                 Utilities::int_to_string (iter, 4) +
                                 ".pvtu");
    std::ofstream master_output ((folder+pvtu_filename).c_str());

    data_out.write_pvtu_record (master_output, filenames);

    const std::string pvdPrefix= (folder + solutionName );
    times_and_names_.push_back (std::pair<double,std::string> (time, pvtu_filename));
    std::ofstream pvd_output (pvdPrefix+".pvd");
    DataOutBase::write_pvd_record (pvd_output, times_and_names_);
  }
}

template <int dim>
void GLSNavierStokesSolver<dim>::write_output_forces ()
{
  TimerOutput::Scope t(computing_timer, "output_forces");
  for (unsigned int boundary_id=0 ; boundary_id<boundaryConditions.size ; ++boundary_id)
  {
    std::string filename = forcesParameters.force_output_name+"."+Utilities::int_to_string (boundary_id, 2)+".dat";
    std::ofstream output (filename.c_str());

    forces_tables[boundary_id].write_text(output);
  }
}

template <int dim>
void GLSNavierStokesSolver<dim>::write_output_torques ()
{
  TimerOutput::Scope t(computing_timer, "output_torques");
  for (unsigned int boundary_id=0 ; boundary_id<boundaryConditions.size ; ++boundary_id)
  {
    std::string filename = forcesParameters.torque_output_name+"."+Utilities::int_to_string (boundary_id, 2)+".dat";
    std::ofstream output (filename.c_str());

    torques_tables[boundary_id].write_text(output);
  }
}

//Find the l2 norm of the error between the finite element sol'n and the exact sol'n
template <int dim>
double GLSNavierStokesSolver<dim>::calculateL2Error()
{
  TimerOutput::Scope t(computing_timer, "error");

  QGauss<dim>  quadrature_formula(degreeQuadrature_+1);
  const MappingQ<dim>      mapping (degreeVelocity_,femParameters.qmapping_all);
  FEValues<dim> fe_values (mapping,
                           fe, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);


  const unsigned int   			dofs_per_cell = fe.dofs_per_cell;         // This gives you dofs per cell
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell); //  Local connectivity

  const unsigned int   n_q_points    = quadrature_formula.size();

  std::vector<Vector<double> > q_exactSol (n_q_points, Vector<double>(dim+1));


  std::vector<Tensor<1,dim> > local_velocity_values (n_q_points);
  std::vector<double > local_pressure_values (n_q_points);

  double l2errorU=0.;

  //loop over elements
  typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
  for (; cell!=endc; ++cell)
  {
    if (cell->is_locally_owned())
    {
      fe_values.reinit (cell);
      fe_values[velocities].get_function_values (present_solution,
                                                 local_velocity_values);
      fe_values[pressure].get_function_values(present_solution,
                                              local_pressure_values);

      //Retrieve the effective "connectivity matrix" for this element
      cell->get_dof_indices (local_dof_indices);

      // Get the exact solution at all gauss points
      exact_solution->vector_value_list(fe_values.get_quadrature_points(),
                                        q_exactSol);

      for(unsigned int q=0; q<n_q_points; q++)
      {
        //Find the values of x and u_h (the finite element solution) at the quadrature points
        double ux_sim=local_velocity_values[q][0];
        double ux_exact=q_exactSol[q][0];

        double uy_sim=local_velocity_values[q][1];
        double uy_exact=q_exactSol[q][1];

        l2errorU += (ux_sim-ux_exact)*(ux_sim-ux_exact) * fe_values.JxW(q);
        l2errorU += (uy_sim-uy_exact)*(uy_sim-uy_exact) * fe_values.JxW(q);

        if (dim==3)
        {
          double uz_sim=local_velocity_values[q][2];
          double uz_exact=q_exactSol[q][2];
          l2errorU += (uz_sim-uz_exact)*(uz_sim-uz_exact) * fe_values.JxW(q);
        }
      }
    }
  }
  l2errorU=Utilities::MPI::sum(l2errorU,mpi_communicator);

  return std::sqrt(l2errorU);
}

// This is a primitive first implementation that could be greatly improved by doing a single pass instead of
// N boundary passes
template<int dim>
void GLSNavierStokesSolver<dim>::calculate_forces()
{
  TimerOutput::Scope t(computing_timer, "calculate_forces");

  QGauss<dim-1>  face_quadrature_formula(degreeQuadrature_+1);
  const MappingQ<dim>      mapping (degreeVelocity_,femParameters.qmapping_all);
  const int n_q_points = face_quadrature_formula.size();
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);
  std::vector<double> pressure_values(n_q_points);
  std::vector<Tensor< 2, dim> > velocity_gradients(n_q_points);
  Tensor<1,dim> normal_vector;
  Tensor<2,dim> fluid_stress;
  Tensor<2,dim> fluid_pressure;
  Tensor<1,dim> force;

  FEFaceValues<dim> fe_face_values (mapping, fe, face_quadrature_formula, update_values
                                    | update_quadrature_points | update_gradients | update_JxW_values
                                    | update_normal_vectors );

  for (unsigned int boundary_id=0 ; boundary_id<boundaryConditions.size ; ++boundary_id)
  {
    force=0;
    typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; face++)
        {
          if (cell->face(face)->at_boundary())
          {
            fe_face_values.reinit (cell, face);
            if (cell->face(face)->boundary_id() == boundary_id)
            {
              std::vector<Point<dim> > q_points = fe_face_values.get_quadrature_points();
              fe_face_values[velocities].get_function_gradients(present_solution, velocity_gradients);
              fe_face_values[pressure].get_function_values(present_solution, pressure_values);
              for (int q = 0; q < n_q_points; q++)
              {
                normal_vector = -fe_face_values.normal_vector(q);
                for (int d=0 ; d<dim ; ++d)
                {
                  fluid_pressure[d][d] = pressure_values[q];
                }
                fluid_stress = viscosity_*(velocity_gradients[q]+transpose(velocity_gradients[q])) - fluid_pressure;
                force += fluid_stress*normal_vector*fe_face_values.JxW(q);
              }
            }
          }
        }
      }
    }
    forces_[boundary_id]=Utilities::MPI::sum(force,mpi_communicator);
  }

  if(forcesParameters.verbosity==forcesParameters.verbose &&  this_mpi_process==0)
  {
    pcout << std::endl;
    TableHandler table;

    for (unsigned int boundary_id=0 ; boundary_id<boundaryConditions.size ; ++boundary_id)
    {
      table.add_value("Boundary ID", boundary_id);
      table.add_value("f_x",forces_[boundary_id][0]);
      table.add_value("f_y",forces_[boundary_id][1]);
      table.set_precision("f_x" ,forcesParameters.display_precision);
      table.set_precision("f_y" ,forcesParameters.display_precision);
      if(dim==3)
      {
        table.add_value("f_z",forces_[boundary_id][2]);
        table.set_precision("f_z" ,forcesParameters.display_precision);
      }
    }
  }

  for (unsigned int boundary_id=0 ; boundary_id<boundaryConditions.size ; ++boundary_id)
  {
    forces_tables[boundary_id].add_value("time",simulationControl.getTime());
    forces_tables[boundary_id].add_value("f_x",forces_[boundary_id][0]);
    forces_tables[boundary_id].add_value("f_y",forces_[boundary_id][1]);
    if(dim==3) forces_tables[boundary_id].add_value("f_z",forces_[boundary_id][2]);
    else forces_tables[boundary_id].add_value("f_z",0.);

    // Precision
    forces_tables[boundary_id].set_precision("f_x" ,forcesParameters.output_precision);
    forces_tables[boundary_id].set_precision("f_y" ,forcesParameters.output_precision);
    forces_tables[boundary_id].set_precision("f_z" ,forcesParameters.output_precision);
    forces_tables[boundary_id].set_precision("time",forcesParameters.output_precision);
  }
}

template<int dim>
void GLSNavierStokesSolver<dim>::calculate_torques()
{
  TimerOutput::Scope t(computing_timer, "calculate_torques");

  QGauss<dim-1>  face_quadrature_formula(degreeQuadrature_+1);
  const MappingQ<dim>      mapping (degreeVelocity_,femParameters.qmapping_all);
  const int n_q_points = face_quadrature_formula.size();
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);
  std::vector<double> pressure_values(n_q_points);
  std::vector<Tensor< 2, dim> > velocity_gradients(n_q_points);
  Tensor<1,dim> normal_vector;
  Tensor<2,dim> fluid_stress;
  Tensor<2,dim> fluid_pressure;

  Tensor<1,dim> force;
  Tensor<1,dim> distance;
  // torque tensor had to be considered in 3D at all time...
  Tensor<1,3> torque;


  FEFaceValues<dim> fe_face_values (mapping, fe, face_quadrature_formula, update_values
                                    | update_quadrature_points | update_gradients | update_JxW_values
                                    | update_normal_vectors );

  for (unsigned int boundary_id=0 ; boundary_id<boundaryConditions.size ; ++boundary_id)
  {
    torque=0;
    Point<dim> center_of_rotation=boundaryConditions.bcFunctions[boundary_id].cor;
    typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; face++)
        {

          if (cell->face(face)->at_boundary())
          {
            fe_face_values.reinit (cell, face);
            if (cell->face(face)->boundary_id() == boundary_id)
            {
              std::vector<Point<dim> > q_points = fe_face_values.get_quadrature_points();
              fe_face_values[velocities].get_function_gradients(present_solution, velocity_gradients);
              fe_face_values[pressure].get_function_values(present_solution, pressure_values);
              for (int q = 0; q < n_q_points; q++)
              {
                normal_vector = -fe_face_values.normal_vector(q);
                for (int d=0 ; d<dim ; ++d)
                {
                  fluid_pressure[d][d] = pressure_values[q];
                }
                fluid_stress = viscosity_*(velocity_gradients[q]+transpose(velocity_gradients[q])) - fluid_pressure;
                force = fluid_stress*normal_vector*fe_face_values.JxW(q);

                distance = q_points[q]-center_of_rotation;
                if (dim==2)
                {
                  torque[0]=0.;
                  torque[1]=0.;
                  torque[2]+=distance[0]*force[1]-distance[1]*force[0];
                }
                else if (dim==3)
                {
                  torque[0] += distance[1] * force[2] - distance[2] * force[1];
                  torque[1] += distance[2] * force[0] - distance[0] * force[2];
                  torque[2] += distance[0] * force[1] - distance[1] * force[0];
                }
              }
            }
          }
        }
      }
    }
    torques_[boundary_id]=Utilities::MPI::sum(torque,mpi_communicator);
  }

  if(forcesParameters.verbosity==forcesParameters.verbose &&  this_mpi_process==0)
  {
    pcout << std::endl;
    TableHandler table;

    for (unsigned int boundary_id=0 ; boundary_id<boundaryConditions.size ; ++boundary_id)
    {
      table.add_value("Boundary ID", boundary_id);
      table.add_value("T_x",torques_[boundary_id][0]);
      table.add_value("T_y",torques_[boundary_id][1]);
      table.set_precision("T_x" ,forcesParameters.display_precision);
      table.set_precision("T_y" ,forcesParameters.display_precision);
      table.add_value("T_z",torques_[boundary_id][2]);
      table.set_precision("T_z" ,forcesParameters.display_precision);
    }

    std::cout << "+------------------------------------------+"  << std::endl;
    std::cout << "|  Torque summary                          |" << std::endl;
    std::cout << "+------------------------------------------+"  << std::endl;
    table.write_text(std::cout);
  }

  for (unsigned int boundary_id=0 ; boundary_id<boundaryConditions.size ; ++boundary_id)
  {
    torques_tables[boundary_id].add_value("time",simulationControl.getTime());
    torques_tables[boundary_id].add_value("T_x",torques_[boundary_id][0]);
    torques_tables[boundary_id].add_value("T_y",torques_[boundary_id][1]);
    torques_tables[boundary_id].add_value("T_z",torques_[boundary_id][2]);

    // Precision
    torques_tables[boundary_id].set_precision("T_x" ,forcesParameters.output_precision);
    torques_tables[boundary_id].set_precision("T_y" ,forcesParameters.output_precision);
    torques_tables[boundary_id].set_precision("T_z" ,forcesParameters.output_precision);
    torques_tables[boundary_id].set_precision("time",forcesParameters.output_precision);
  }
}

#endif
