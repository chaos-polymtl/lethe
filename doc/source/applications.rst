######################
Applications Overview
######################

Lethe has several applications that can be used depending on the problem to be solved. Once Lethe is compiled, all the executables to these applications become available. In the following table, a brief description of each application is provided:

.. list-table::
    :header-rows: 1
    :widths: 40 40 60 

    * - Physics 
      - Application
      - Description
    * - * Single-phase flows
        * Fluid-fluid flows 
        * Heat transfer
        * Tracer
      - ``lethe-fluid``
      - This application solves the Navier-Stokes equations in a monolithic way. **This is the most robust fluid dynamics solver in Lethe and it is compatible with all multiphysics features**.
    * - * Single-phase flows
        * Heat transfer
        * Tracer
      - ``lethe-fluid-matrix-free``
      - This application solves the Navier-Stokes equations in a monolithic way and using a matrix-free approach.  **This is the fastest fluid dynamics solver in Lethe. It should be used for high-order simulations. It does not fully support all multiphysics features**.
    * - * Single-phase flows
        * Fluid-fluid flows 
        * Heat transfer
        * Tracer
        * Immersed boundaries
      - ``lethe-fluid-nitsche``
      - This application solves the Navier-Stokes equations by imposing immersed boundary conditions without discretizing the boundaries using the Nitsche method. 
    * - * Particle-laden flows
      - ``lethe-fluid-particles``
      - This application allows to run unresolved Computational Fluid Dynamics-Discrete Element Method simulations.
    * - * Particle-laden flows
        * Immersed boundaries
      - ``lethe-fluid-sharp``
      - This application solves the Navier-Stokes equations for particle-laden flows. A sharp-edge method is used to impose immersed boundary conditions at the particle walls, removing the requirement of a boundary-conforming mesh. This application allows to perform resolved Computational Fluid Dynamics-Discrete Element Method simulations.
    * - * Single phase volume-averaged flows
      - ``lethe-fluid-vans``
      - This application solves the Volume-Averaged Navier-Stokes equations in a monolithic way.
    * - * Granular flow of spherical particles
      - ``lethe-particles``
      - This application uses the Discrete Element Method to simulate spherical particles. 
    * - * Steady-state single-phase flows
      - ``lethe-fluid-block``
      - This application solves the Navier-Stokes equations in a monolithic way using block preconditioner. **This solver is experimental and is only adequate for steady-state simulations**.
    * - * Radioactive Particle Tracking
      - ``lethe-rpt-3d``
      - This application solves the Beam Monte-Carlo model for the interaction between gamma ray emitting particles and detectors.
    * - * Radioactive Particle Tracking
      - ``lethe-rpt-cell-reconstruction-3d``
      - This application reconstructs the position of particles using a cellular decomposition approach.
    * - * Radioactive Particle Tracking
      - ``lethe-rpt-fem-reconstruction-3d``
      - This application allows to reconstruct positions from ray tracing experimental measurements obtained by Radioactive Particle Tracking.
    * - * Radioactive Particle Tracking
      - ``lethe-rpt-l2-projection-3d``
      - This application calculates a volumetric projection of the Beam model on a finite element mesh using an L2 projection.