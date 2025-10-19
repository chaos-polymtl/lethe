***********************************************
Discrete Element Method (DEM)
***********************************************
The discrete element method (DEM) solver in Lethe is called ``lethe-particles`` and supports two-dimensional and three-dimensional simulations.
In ``simulation_control``, the general simulation parameters (for example, time-step, end time) are defined. The ``timer`` and ``test`` sections are used for timing the classes and functions, and testing the reproducibility of the results.
In ``model parameters``, we define the simulation models, including particle-particle and particle-wall contact models.
In the ``lagrangian physical properties`` section,  particles and walls physical properties are defined.
The ``mesh`` section defines the simulation triangulation and refinements.
The ``insertion info`` section defines every parameters related to each insertion method available in Lethe, i.e. the insertion box and insertion frequency for the volume insertion method.
In the ``floating walls`` section, hyperplanes defined as mathematical planes can be added to the simulation domain to create walls in the simulation domain.
In the ``solid objects`` section, surface meshes can be defined and their motion can be controlled during a simulation.
In the ``boundary conditions`` section, simulation domain boundaries can be defined as walls, outlets or they can get periodicity, rotation or translation.
The ``post-processing``` can activate the post-processing lagrangian and force chains visualization tools.
The ``particle ray tracing`` section defines the parameters used by the ``lethe-particles-ray-tracing`` application to simulate the propagation of photons (or rays) through the simulation domain to reconstruct the surface made by particle using the same principles as profilometry. This applications uses a lot of the same section as for the DEM solver, such as ``mesh``, and ``insertion_info`` which is why it is included here.

.. toctree::
   :maxdepth: 1

   simulation_control
   test
   timer
   model_parameters
   lagrangian_physical_properties
   mesh
   insertion_info
   floating_wall
   solid_objects
   boundary_conditions
   post-processing
   particle_ray_tracing