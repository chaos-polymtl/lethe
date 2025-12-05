==================================
Gas-Solid Fluidized Bed
==================================

This example simulates the fluidization of particles in air within a cylindrical bed. It is based on the fluidized bed test case presented in the work of El Geitani *et al* [#ElGeitani2023]_ and was used to validate one of the earliest CFD-DEM implementations in Lethe. Most importantly, this example compares the pressure drop across the bed as a function of the superficial gas velocity with correlations available in the literature.

----------------------------------
Features
----------------------------------

- Solvers: ``lethe-particles``, ``lethe-fluid-particles`` and ``lethe-fluid-particles-matrix-free``, with Q1-Q1
- Three-dimensional problem
- Displays the selection of models and physical properties
- Simulates a solid-gas fluidized bed


---------------------------
Files Used in This Example
---------------------------

All files mentioned below are located in the example's folder (``examples/unresolved-cfd-dem/gas-solid-fluidized-bed-cylindrical``).

- Parameter file for particle generation and packing: ``particle-packing.prm``
- Parameter file for CFD-DEM simulation of the gas-solid fluidized bed: ``mb-fluidized-bed-modelA.prm``, ``mb-fluidized-bed-modelB.prm``, ``mb-fluidized-bed-modelA-project.prm``, ``mb-fluidized-bed-modelB-project.prm`` and ``mf-fluidized-bed-modelA.prm``
- Post-processing Python script: ``plot-pressure.py``


-----------------------
Description of the Case
-----------------------

This example simulates a gas–solid fluidized bed inside a cylindrical column (diameter 0.02 m, height 0.4 m). First, ``lethe-particles`` is used with ``particle-packing.prm`` to generate and pack spherical particles (diameter 0.0005 m, density 2500 kg·m⁻³) inside the column. After packing, the solid–fluid mixture is simulated with two solvers: the matrix-free solver ``lethe-fluid-particle-matrix-free`` and the matrix-based CFD–DEM solver ``lethe-fluid-particles``. For the matrix-based solver we test two VANS models (model A and model B), and for each model we project particle–fluid forces using one of two approaches: a cell-based filter or the Quadrature-Centered Method (QCM) filter. The superficial gas velocity at the inlet is varied from 0.02 to 0.4 m·s⁻¹ and the pressure drop across the bed is recorded. Results are compared with correlations from the literature for validation


-------------------
DEM Parameter File
-------------------

A detailed description of all the DEM parameter subsections can be found in the `parameter section <../../../parameters/dem.html>`_. The subsections in the DEM parameter file ``particle-packing.prm`` that are pertinent to this example are described below. 

Mesh
~~~~~

As mentioned in the example description, the particles are packed inside a cylindrical column. For this reason, the mesh type is set to ``cylinder`` with a ``balanced`` grid. This mesh uses the same input arguments as the ``GridGenerator::subdivided_cylinder`` function of Deal.II, yet leads to more uniform cells across the domain. An initial refinement level of 2 provides enough cells for the CFD solver while keeping the smallest cell size larger than the particle diameter. Finally, the particle–wall contact search expansion is enabled to ensure proper detection of particle–wall interactions in the curved convex geometry.

.. code-block:: text

    subsection mesh
        set type               = cylinder
        set grid type          = balanced
        set grid arguments     = 44:0.01:0.22
        set initial refinement = 2
        set expand particle-wall contact search = true
    end
    
Simulation Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Restart
~~~~~~~~~~~~~~~~~~~



Model Parameters
~~~~~~~~~~~~~~~~~



Lagrangian Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    
Insertion Info
~~~~~~~~~~~~~~~~~~~




Floating Walls
~~~~~~~~~~~~~~~~~~~



---------------------------
Running the DEM Simulation
---------------------------



-----------------------
CFD-DEM Parameter File
-----------------------



Simulation Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~



Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~



Initial Conditions
~~~~~~~~~~~~~~~~~~


 

Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Void Fraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~



CFD-DEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Non-linear Solver
~~~~~~~~~~~~~~~~~


    
Linear Solver
~~~~~~~~~~~~~




------------------------------
Running the CFD-DEM Simulation
------------------------------



--------
Results
--------


    
----------
Reference
----------

.. [#ElGeitani2023] T. El Geitani, S. Golshan, and B. Blais, “Toward High-Order CFD-DEM: DevelopmentandValidation,” *Industrial & Engineering Chemistry Research*, vol. 62, no. 2, pp. 1141–1159, January 2023. Available: `<https://doi.org/10.1021/acs.iecr.2c03546>`_.
    
