######################
Applications Overview
######################

Lethe has several applications that can be used depending on the problem that wants to be solved. Once Lethe is compiled, all the executables to these applications become available. In the following table a brief description of each application is provided:

.. list-table::
    :header-rows: 1
    :widths: 40 60 

    * - Application
      - Description
    * - ``lethe-fluid``
      - This application solves the Navier-Stokes equations in a monolithic way. **This is the most robust fluid dynamics solver in Lethe**.
    * - ``lethe-fluid-block``
      - This application solves the Navier-Stokes equations in a monolithic way using block preconditioner.
    * - ``lethe-fluid-matrix-free``
      - This application solves the Navier-Stokes equations in a monolithic way and using a matrix-free approach. 
    * - ``lethe-fluid-nitsche``
      - This application solves the Navier-Stokes equations by imposing immersed boundary conditions without discretizing the boundaries using the Nitsche method. 
    * - ``lethe-fluid-particles``
      - This application allows to run unresolved Computational Fluid Dynamics-Discrete Element Method simulations.
    * - ``lethe-fluid-sharp``
      - This application solves the Navier-Stokes equations by imposing immersed boundary conditions without discretizing the boundaries using the sharp-edge method. 
    * - ``lethe-fluid-vans``
      - This application solves the Volume-Averaged Navier-Stokes equations in a monolithic way.
    * - ``lethe-particles``
      - This application uses the Discrete Element Method to simulate particles. 
    * - ``lethe-rpt-3d``
      - This application solves the Beam Monte-Carlo model for the interaction between gamma ray emitting particles and detectors.
    * - ``lethe-rpt-cell-reconstruction-3d``
      - This application reconstructs the position of particles using a cellular decomposition approach.
    * - ``lethe-rpt-fem-reconstruction-3d``
      - This application allows to reconstruct positions from ray tracing experimental measurements obtained by Radioactive Particle Tracking.
    * - ``lethe-rpt-l2-projection-3d``
      - This application calculates a volumetric projection of the Beam model on a finite element mesh using an L2 projection.