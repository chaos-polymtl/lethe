######################
Applications Overview
######################

Lethe has several applications that can be used depending on the problem to be solved. Once Lethe is compiled, all the executables to these applications become available. In the following table, we provide a brief description of each application. The applications are sorted from the most commonly used to the least.

.. warning::
  
  The radioactive particle-tracking (RPT) applications of Lethe has been migrated to a separate repository which is available `here <https://github.com/chaos-polymtl/lethe-rpt>`_.

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
   * - ``lethe-particles-ray-tracing``
     - * Profilometry extraction of particle surfaces
     -  This application performs particle ray tracing to reconstruct the geometry of the surface of an assembly of particles, similar to a virtual profilometry technique. It emits rays through the simulation domain and detects particle intersections to represent the surface.
   * - ``lethe-fluid-particles``
     - * Particle-laden flows
     - This application runs unresolved Computational Fluid Dynamics-Discrete Element Method simulations.
   * - ``lethe-fluid-particles-matrix-free``
     - * Particle-laden flows
     - This application runs unresolved Computational Fluid Dynamics-Discrete Element Method simulations using a matrix-free approach for the fluid dynamics. **This solver is still experimental**.
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

