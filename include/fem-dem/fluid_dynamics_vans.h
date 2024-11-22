// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fluid_dynamics_vans_h
#define lethe_fluid_dynamics_vans_h

#include <core/bdf.h>
#include <core/dem_properties.h>
#include <core/grids.h>
#include <core/manifolds.h>
#include <core/parameters.h>
#include <core/time_integration_utilities.h>

#include <solvers/fluid_dynamics_matrix_based.h>
#include <solvers/postprocessing_cfd.h>

#include <dem/dem.h>
#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/vans_assemblers.h>
#include <fem-dem/void_fraction.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/property_pool.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>


using namespace dealii;

/**
 * A solver class for the VANS equation using GLS stabilization
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Toni EL Geitani, 2020
 */

template <int dim>
class FluidDynamicsVANS : public FluidDynamicsMatrixBased<dim>
{
public:
  FluidDynamicsVANS(CFDDEMSimulationParameters<dim> &nsparam);

  ~FluidDynamicsVANS();

  virtual void
  solve() override;

private:
  void
  assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix);

  void
  update_solution_and_constraints();

  void
  read_dem();

  /**
   * @brief This function calculates and returns the periodic offset distance of the domain which is needed
   * for the periodic boundary conditions using the QCM or SPM for void fraction
   * with the GLS VANS/CFD-DEM solver. The distance is based on one of the
   * periodic boundaries and all particle location shifted by this distance is
   * according to this periodic boundary.
   *
   * @param boundary_id The id of one of the periodic boundaries
   *
   * @return The periodic offset distance
   */
  inline Tensor<1, dim>
  get_periodic_offset_distance(unsigned int boundary_id) const
  {
    Tensor<1, dim> offset;

    // Iterating over the active cells in the triangulation
    for (const auto &cell : (*this->triangulation).active_cell_iterators())
      {
        if (cell->is_locally_owned() || cell->is_ghost())
          {
            if (cell->at_boundary())
              {
                // Iterating over cell faces
                for (unsigned int face_id = 0; face_id < cell->n_faces();
                     ++face_id)
                  {
                    unsigned int face_boundary_id =
                      cell->face(face_id)->boundary_id();

                    // Check if face is on the boundary, if so, get
                    // the periodic offset distance for one pair of periodic
                    // faces only since periodic boundaries are aligned with the
                    // direction and only axis are currently allowed
                    if (face_boundary_id == boundary_id)
                      {
                        Point<dim> face_center = cell->face(face_id)->center();
                        auto periodic_cell = cell->periodic_neighbor(face_id);
                        unsigned int periodic_face_id =
                          cell->periodic_neighbor_face_no(face_id);
                        Point<dim> periodic_face_center =
                          periodic_cell->face(periodic_face_id)->center();

                        offset = periodic_face_center - face_center;

                        return offset;
                      }
                  }
              }
          }
      }

    // A zero tensor is returned in case no cells are found on the periodic
    // boundaries on this processor. This processor won't handle particle in
    // cells at periodic boundaries, so it won't affect any computation.
    return offset;
  }

protected:
  /**
   * @brief associates the degrees of freedom to each vertex of the finite elements
   * and initialize the void fraction
   */
  virtual void
  setup_dofs() override;

  virtual void
  iterate() override;

  void
  calculate_void_fraction(const double time, bool load_balance_step);

  void
  vertices_cell_mapping();

  virtual void
  monitor_mass_conservation();

  /**
   * @brief finish_time_step
   * Finishes the time step
   * Post-processing and time stepping
   */

  virtual void
  finish_time_step_fd();

  /**
   *  @brief Assembles the matrix associated with the solver
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assembles the rhs associated with the solver
   */
  void
  assemble_system_rhs() override;

  /**
   * @brief Assembles the local matrix for a given cell.
   *
   * This function is used by the WorkStream class to assemble
   * the system matrix. It is a thread safe function.
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim>                         &scratch_data,
    StabilizedMethodsTensorCopyData<dim>                 &copy_data) override;

  /**
   * @brief Assembles the local rhs for a given cell
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim>                         &scratch_data,
    StabilizedMethodsTensorCopyData<dim>                 &copy_data) override;

  /**
   * @brief sets up the vector of assembler functions
   */
  void
  setup_assemblers() override;


  /**
   * @brief Copies local cell information to global matrix
   */

  void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Copies local cell rhs information to global rhs
   */

  void
  copy_local_rhs_to_global_rhs(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief a function for adding data vectors to the data_out object for
   * post_processing additional results
   */
  virtual void
  output_field_hook(DataOut<dim> &data_out) override;



  /**
   * Member Variables
   */

  CFDDEMSimulationParameters<dim> cfd_dem_simulation_parameters;

  MappingQGeneric<dim> particle_mapping;

  IndexSet locally_owned_dofs_voidfraction;
  IndexSet locally_relevant_dofs_voidfraction;

  // Assemblers for the particle_fluid interactions
  std::vector<std::shared_ptr<ParticleFluidAssemblerBase<dim>>>
    particle_fluid_assemblers;


  const bool   PSPG        = true;
  const bool   SUPG        = true;
  const double GLS_u_scale = 1;
  double       pressure_drop;

  bool           has_periodic_boundaries;
  Tensor<1, dim> periodic_offset;
  unsigned int   periodic_direction;

  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_cell;
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_periodic_cell;

protected:
  Particles::ParticleHandler<dim, dim> particle_handler;

  VoidFractionBase<dim> void_fraction_manager;
};

#endif
