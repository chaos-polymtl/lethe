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
      - This application solves the Navier-Stokes equations in a monolithic way.
    * - ``lethe-fluid-block``
      - This application solves the Navier-Stokes equations in a monolithic way with a block approach.
    * - ``lethe-fluid-matrix-free``
      - This application solves the Navier-Stokes equations in a monolithic way and using a matrix-free approach. 
    * - ``lethe-fluid-nitsche``
      - This application solves the Navier-Stokes equations by imposing boundary conditions without discretizing the boundaries. 
    * - ``lethe-fluid-particles``
      - This application allows to run unresolved Computational Fluid Dynamics-Discrete Element Method simulations.
    * - ``lethe-fluid-sharp``
      - This application solves the Navier-Stokes equations by imposing boundary conditions without discretizing the boundaries. 
    * - ``lethe-fluid-vans``
      - This application solves the Volume-Averaged Navier-Stokes equations in a monolithic way.
    * - ``lethe-particles``
      - This application uses the Discrete Element Method to simulate particles motion.
    * - ``lethe-rpt-3d``
      - TODO: add description
    * - ``lethe-rpt-cell-reconstruction-3d``
      - TODO: add description
    * - ``lethe-rpt-fem-reconstruction-3d``
      - This application allows to reconstruct positions from ray tracing experimental measurements obtained by Radioactive Particle Tracking.
    * - ``lethe-rpt-l2-projection-3d``
      - TODO: add description