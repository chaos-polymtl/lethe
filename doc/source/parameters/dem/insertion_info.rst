Insertion Info
-------------------
Particle insertion information is defined in this section. This information includes ``insertion method``, ``inserted number of particles at each insertion step``, ``insertion frequency``, dimensions of the insertion box, approximate initial distance between particles at insertion steps, and information about the randomness of initial positions of particles. If the insertion method is ``uniform``, the particles are uniformly (without randomness in initial location) inserted. On the other hand, in ``non_uniform`` insertion, initial locations of particles in `x` and `y` directions are chosen randomly. Using ``non_uniform`` insertion, insertion ``random number range`` and ``insertion random number seed`` are also required in the parameter handler file.

.. note::
    Insertion in Lethe-DEM starts inserting particles from type 0 and proceeds to the next type when all the particles from the previous type are inserted.

.. note::
    If the insertion box is not large enough to insert ``inserted number of particles at each time step``, Lethe-DEM prints a warning to inform the user about the number of particles that can fit in the insertion box.

.. note::
    ``insertion distance threshold`` determines the initial distance between the particles in the insertion phase. As a result, it must be larger than 1 to avoid any initial collision between the inserted particles. It should also be compatible with the ``random number range``; especially if the ``random number range`` is large, a large value should be defined for ``insertion distance threshold``. Generally, we recommend users to use a value in the range of 1.3-2 (depending on the value of ``random number range``) for ``insertion distance threshold``.

.. code-block:: text

 subsection insertion info
  # Insertion method
  # Choices are uniform|non_uniform
  set insertion method				        = non_uniform

  # Number of inserted particles at each insertion step. This value may change automatically if the insertion box is not adequately large to handle all the inserted particles
  set inserted number of particles at each time step    = 100

  # Insertion frequency
  set insertion frequency            		 	= 20000

  # Insertion box dimensions (xmin, xmax, ymin, ymax, zmin, zmax)
  set insertion box minimum x            	 	= -0.05
  set insertion box minimum y            	        = -0.05
  set insertion box minimum z            	        = -0.03
  set insertion box maximum x            	        = 0.05
  set insertion box maximum y           	 	= 0.05
  set insertion box maximum z            	        = 0.07

  # This value controls the distance between the center of inserted particles (the distance is: [distance threshold] * [diameter of particles]). The distance is modified by a random number if non_uniform insertion is chosen
  set insertion distance threshold			= 2

  # Random number range and seed for non_uniform insertion
  set insertion random number range			= 0.75
  set insertion random number seed			= 19
 end