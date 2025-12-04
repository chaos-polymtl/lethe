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




-------------------
DEM Parameter File
-------------------


Mesh
~~~~~

    
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
    
