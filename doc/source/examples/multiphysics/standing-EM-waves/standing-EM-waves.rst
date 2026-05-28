=============================
Waveguide
=============================

This example showcases an electromagnetic problem that we can solve analytically, allowing us to test the results of our simulation. It is the first example using the DPG method (Discontinuous Petrov-Gallerkin method) of the time-harmonic section of lethe.

-----------------------------
Features
-----------------------------

- Solvers: ``lethe-fluid``
- Steady-state problem



----------------------------
Files Used In This Example
----------------------------

All files mentioned below are located in the example's folder (``examples/multiphysics/waveguide``)

- Base case parameter file (:math:`\mathrm{f=2.45 GHz}`): ``parameters.prm``
- Postprocessing Python script for the :math:`\mathrm{Re}=400` case: ``convergence.py``


---------------------------
Description of the Case
---------------------------

The waveguide at stake has a rectangular cross section.It consists of four perfectly reflective walls and two empty sections: one for the wave’s entry and the other for its exit. The wave propagates along the z-axis. It is described on the next figure:

.. image:: images/schematic.png
    :alt: schematic
    :align: center
    :name: geometry
    :width: 500

    
--------------
Parameter File
--------------

-----------------------
Running the Simulations
-----------------------

-----------------------
Results and Discussion
-----------------------

.. image:: images/model_validation.png
    :alt: final graph
    :align: center
