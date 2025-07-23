######################
Applications Overview
######################

Lethe has several applications that can be used depending on the problem to be solved. Once Lethe is compiled, all the executables to these applications become available. In the following table, we provide a brief description of each application. The applications are sorted from the most commonly used to the least.

.. list-table::
   :header-rows: 1
   :widths: 40 40 60 

   * - Application
     - Physics
     - Description
   * - ``lethe-fluid``
     - * Single-phase flows
       * Fluid-fluid flows 
       * Heat transfer
       * Passive tracer
     - This application solves the Navier-Stokes equations in a monolithic way. **This is the most robust fluid dynamics solver in Lethe and it is compatible with all multiphysics features**.
   * - ``lethe-fluid-matrix-free``
     - * Single-phase flows
       * Heat transfer
       * Passive tracer
     - This application solves the Navier-Stokes equations in a monolithic way and using a matrix-free approach. **This is the fastest fluid dynamics solver in Lethe. It should be used for high-order simulations. It does not fully support all multiphysics features**.
   * - ``lethe-particles``
     - * Granular flows of spherical particles
     - This application uses the Discrete Element Method to simulate spherical particles.
   * - ``lethe-fluid-particles``
     - * Particle-laden flows
     - This application allows to run unresolved Computational Fluid Dynamics-Discrete Element Method simulations.
   * - ``lethe-fluid-sharp``
     - * Single-phase flows
       * Particle-laden flows
       * Passive tracer  
     - This application solves the Navier-Stokes equations for particle-laden flows. A sharp-edge method is used to impose immersed boundary conditions at the particle walls, removing the requirement of a boundary-conforming mesh. This application allows to perform resolved Computational Fluid Dynamics-Discrete Element Method simulations.
   * - ``lethe-fluid-nitsche``
     - * Single-phase flows
       * Fluid-fluid flows 
       * Heat transfer
       * Passive tracer
     - This application solves the Navier-Stokes equations by imposing immersed boundary conditions without discretizing the boundaries using the Nitsche method.
   * - ``lethe-fluid-vans``
     - * Single phase volume-averaged flows
     - This application solves the Volume-Averaged Navier-Stokes equations in a monolithic way. It supports calculation of the void fraction from particles, but the particles remain static.
   * - ``lethe-fluid-block``
     - * Steady-state single-phase flows
     - This application solves the Navier-Stokes equations in a monolithic way using block preconditioner. **This solver is experimental and is only adequate for steady-state simulations**.
   * - ``lethe-rpt-3d``
     - * Radioactive Particle Tracking
     - This application solves the Beam Monte-Carlo model for the interaction between gamma ray emitting particles and detectors.
   * - ``lethe-rpt-cell-reconstruction-3d``
     - * Radioactive Particle Tracking
     - This application reconstructs the position of particles using a cellular decomposition approach.
   * - ``lethe-rpt-fem-reconstruction-3d``
     - * Radioactive Particle Tracking
     - This application allows to reconstruct positions from ray tracing experimental measurements obtained by Radioactive Particle Tracking.
   * - ``lethe-rpt-l2-projection-3d``
     - * Radioactive Particle Tracking
     - This application calculates a volumetric projection of the Beam model on a finite element mesh using an L2 projection.