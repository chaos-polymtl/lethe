// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_subequations_interface_h
#define lethe_subequations_interface_h

#include <core/vector.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/physics_subequations_solver.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/trilinos_vector.h>


using namespace dealii;

template <int dim>
class SubequationsInterface
{
public:
  /**
   * @brief Default constructor of
   *
   * @param nsparam
   * @param p_triangulation
   * @param p_simulation_control
   * @param p_pcout
   */
  SubequationsInterface(
    const SimulationParameters<dim> &sim_param,
    MultiphysicsInterface<dim>      *p_multiphysics,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                       p_triangulation,
    std::shared_ptr<SimulationControl> p_simulation_control,
    ConditionalOStream                &p_pcout);

  /**
   * @brief Default destructor.
   */
  ~SubequationsInterface() = default;


  /**
   * @brief Setup the DofHandler and the degree of freedom associated with the
   * physics.
   */
  void
  setup_dofs()
  {
    for (const auto &subequation : this->subequations)
      {
        subequation.second->setup_dofs();
      }
  };

  /**
   * @brief
   */
  void
  solve()
  {
    for (const auto &subequation : subequations)
      {
        subequation.second->solve();
      }
  }

  /**
   * @brief
   *
   * @return
   */
  std::vector<SubequationsID>
  get_active_subequations()
  {
    return this->active_subequations;
  }

  /**
   * @brief
   *
   * @param[in] subequation_id
   *
   * @return
   */
  DoFHandler<dim> *
  get_dof_handler(const SubequationsID subequation_id)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());

    return this->subequations_dof_handler[subequation_id];
  }

  /**
   * @brief
   *
   * @param[in] subequation_id
   *
   * @return
   */
  GlobalVectorType *
  get_solution(const SubequationsID subequation_id)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());
    return subequations_solutions[subequation_id];
  }

  /**
   * @brief
   *
   * @param[in] subequation_id
   *
   * @param[in] dof_handler
   */
  void
  set_dof_handler(const SubequationsID subequation_id,
                  DoFHandler<dim>     *dof_handler)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());
    subequations_dof_handler[subequation_id] = dof_handler;
  }

  /**
   * @brief
   *
   * @param[in] subequation_id
   *
   * @param[in] solution_vector
   */
  void
  set_solution(const SubequationsID subequation_id,
               GlobalVectorType    *solution_vector)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());
    subequations_solutions[subequation_id] = solution_vector;
  }



private:
  MultiphysicsInterface<dim> *multiphysics;

  std::map<SubequationsID, Parameters::Verbosity> verbosity;
  ConditionalOStream                              pcout;

  // Data structure that stores all enabled physics
  std::vector<SubequationsID> active_subequations;

  // Subequations stored within a map of shared pointer to ensure proper
  // memory management.
  std::map<
    SubequationsID,
    std::shared_ptr<PhysicsLinearSubequationsSolver<dim, GlobalVectorType>>>
    subequations;

  std::map<SubequationsID, DoFHandler<dim> *> subequations_dof_handler;

  // Present solutions
  std::map<SubequationsID, GlobalVectorType *> subequations_solutions;
};



#endif
