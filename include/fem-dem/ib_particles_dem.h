//
// Created by lucka on 2021-10-11.
//

#ifndef LETHE_IB_PARTICLES_DEM_H
#define LETHE_IB_PARTICLES_DEM_H

#include <core/ib_particle.h>
#include <core/ib_stencil.h>
#include <core/lethegridtools.h>

#include <solvers/navier_stokes_base.h>

#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

using namespace dealii;
using namespace std;


/**
 * A solver class for the DEM used in conjunction with IB particles and
 * gls_sharp_navier_stokes. This class defines and use some function of the DEM
 * class that has been modified and simplified to use IB_particles.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Lucka Barbeau, Bruno Blais, Shahab Golshan 2021
 */

template <int dim>
class IBParticlesDEM
{
public:
  /**
   * @brief
   * Initialize the IBParticlesDEM object with the parameter, the mpi
   * communicator and the particles.
   *
   * @param p_nsparam The parameter of the simulation.
   *
   * @param mpi_communicator_input The mpi communicator of the simulation.
   *
   * @param particles The particles vector containing all the IB particles.
   */
  void
  initialize(SimulationParameters<dim>    p_nsparam,
             MPI_Comm &                   mpi_communicator_input,
             std::vector<IBParticle<dim>> particles);


  /**
   * @brief
   * update the boundary cells that are contact candidate for each of the
   * particle.
   *
   * @param particles The particles vector containing all the IB particles.
   */
  void
  update_particles(std::vector<IBParticle<dim>> particles, double time);


  /**
   * @brief
   * Integrate the dynamics of the IB_particle taking into account the contact
   * between particles and between particles and walls.
   * @param dt The CFD time step.
   *
   */
  void
  particles_dem(double dt);

  /**
   * @brief Calculate non-linear (Hertzian) particle-particle contact force
   *
   * @param dt_dem The sub time stepping time step.
   *
   * @param contact_force a vector containing the contact force between particles
   *
   * @param contact_force a vector containing the contact torques between particles
   */
  void
  calculate_pp_contact_force(const double &               dt_dem,
                             std::vector<Tensor<1, dim>> &contact_force,
                             std::vector<Tensor<1, 3>> &  contact_torque);


  /**
   * @brief Calculate non-linear (Hertzian) particle-wall contact force
   *
   * @param dt_dem The sub time stepping time step.
   *
   * @param contact_force a vector containing the contact force between particles
   *
   * @param contact_force a vector containing the contact torques between particles
   */
  void
  calculate_pw_contact_force(const double &               dt_dem,
                             std::vector<Tensor<1, dim>> &contact_force,
                             std::vector<Tensor<1, 3>> &  contact_torque);

  /**
   * @brief  Update the boundary cells that are contact candidate for each of the particle.
   *
   * @param particles The particles vector containing all the IB particles.
   *
   * @param dof_handler The dof handler of the mesh used for the fluid simulation.
   */

  void
  update_particles_boundary_contact(
    std::vector<IBParticle<dim>> &particles,
    DoFHandler<dim> &             dof_handler,
    const Quadrature<dim - 1> &   face_quadrature_formula,
    const Mapping<dim> &          mapping);


  std::vector<IBParticle<dim>> dem_particles;

private:
  // A struct to store contact tangential history
  struct ContactTangentialHistory
  {
    Tensor<1, dim> tangential_relative_velocity;
    Tensor<1, dim> tangential_overlap;
  };
  // Particles contact history

  // A struct to store boundary cells' information
  struct BoundaryCellsInfo
  {
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
      for (unsigned int i = 0; i < dim; ++i)
        {
          ar &normal_vector;
          ar &point_on_boundary;
        }
    }

    Tensor<1, dim> normal_vector;
    Point<dim>     point_on_boundary;
  };

  /** This function is used to find the projection of vector_a on
   * vector_b
   * @param vector_a A vector which is going to be projected on vector_b
   * @param vector_b The projection vector of vector_a
   * @return The projection of vector_a on vector_b
   */
  inline Tensor<1, dim>
  find_projection(const Tensor<1, dim> &vector_a,
                  const Tensor<1, dim> &vector_b)
  {
    Tensor<1, dim> vector_c;
    vector_c = ((vector_a * vector_b) / (vector_b.norm_square())) * vector_b;

    return vector_c;
  }

  SimulationParameters<dim> parameters;

  MPI_Comm mpi_communicator;

  // Particles contact history
  std::map<unsigned int, std::map<unsigned int, ContactTangentialHistory>>
    pp_contact_map;
  std::map<unsigned int, std::map<unsigned int, ContactTangentialHistory>>
    pw_contact_map;

  // A vector of vectors of candidate cells for each of the particle.
  std::vector<std::vector<BoundaryCellsInfo>> boundary_cells;

private:
  double cfd_time;
};


#endif // LETHE_IB_PARTICLES_DEM_H
