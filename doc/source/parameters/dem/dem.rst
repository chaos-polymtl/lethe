***********************************************
Discrete Element Method (DEM)
***********************************************
There are two discrete element method (DEM) solvers in Lethe: ``dem_2d`` and ``dem_3d`` for two-dimensional and three-dimensional simulations. In ``simulation_control``, the general simulation parameters (for example, time-step, end time, etc.) are defined. ``timer`` and ``test`` sections are used for timing the classes and functions, and testing the reproducibility of the results. In ``model_parameters``, we define the simulation models, including particle-particle and particle-wall contact models. ``lagrangian_physical_properties`` defines the physical properties of the particles and walls. Insertion information including the dimensions of the insertion box and insertion frequency is defined in the ``insertion_info`` section. In the ``floating_walls`` section, we can add flat walls to the simulation domain. ``mesh`` section defines the simulation triangulation and refinements. In the ``boundary_conditions`` section, we can define boundaries as outlets, or rotate or slide them.

.. toctree::
   :maxdepth: 1

   simulation_control
   timer
   test
   model_parameters
   lagrangian_physical_properties
   insertion_info
   floating_walls
   mesh
   boundary_conditions