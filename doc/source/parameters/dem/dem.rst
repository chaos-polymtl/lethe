***********************************************
Discrete Element Method (DEM)
***********************************************
The discrete element method (DEM) solver in Lethe is called ``lethe-particles`` and supports two-dimensional and three-dimensional simulations. In ``simulation_control``, the general simulation parameters (for example, time-step, end time, etc.) are defined. ``timer`` and ``test`` sections are used for timing the classes and functions, and testing the reproducibility of the results. 
In ``model parameters``, we define the simulation models, including particle-particle and particle-wall contact models. ``lagrangian physical properties`` defines the physical properties of the particles and walls. 
Insertion information including the dimensions of the insertion box, insertion frequency, etc. are defined in the ``insertion info`` section. 
The``mesh`` section defines the simulation triangulation and refinements. 
In the ``DEM boundary conditions`` section, boundaries can be defined as outlets or they can get periodicity, rotation or translation.  
In the ``floating walls`` section, hyperplanes can be added to the simulation domain. 
In the ``solid objects`` section, surface meshes can defined and their motion can be controlled.

.. toctree::
   :maxdepth: 1

   boundary_conditions
   floating_wall
   insertion_info
   lagrangian_physical_properties
   mesh
   model_parameters
   post-processing
   simulation_control
   solid_objects
   test
   timer