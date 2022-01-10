Lagrangian Physical Properties
----------------------------------
In this subsection, gravitational acceleration, and the physical properties of the particles and walls are defined. These properties include ``number of particle types``, and for each type, particle ``diameter``, particle ``density``, ``Young's modulus`` of particle and wall, ``Poisson ratio`` of particle and wall, ``restitution coefficient`` of particle and wall, ``friction coefficient`` of particle and wall and ``rolling friction coefficient`` of particle and wall.

.. code-block:: text

 subsection lagrangian physical properties
  # Gravitational acceleration in x direction
  set gx            		 			 = 0.0

  # Gravitational acceleration in y direction
  set gy            		 			 = 0.0

  # Gravitational acceleration in z direction
  set gz            		 			 = -9.81

  # Number of particle types
  set number of particle types             = 1

  # Entering particle type 0
    	subsection particle type 0

  # Size distribution of particle type 0
		set size distribution type              = uniform

  # Particle diameter
        set diameter                            = 0.005

  # Number of particles in type 0
		set number                              = 132300

  # Particle density
        set density particles                   = 2000

  # Young's modulus of particle
        set young modulus particles             = 1000000

  # Poisson ratio of particle
        set poisson ratio particles             = 0.3

  # Coefficient of restitution of particle
        set restitution coefficient particles   = 0.95

  # Coefficient of friction of particle
        set friction coefficient particles      = 0.05

  # Coefficient of rolling friction of particle
        set rolling friction particles          = 0.1
	end

  # Young's modulus of wall
  set young modulus wall                                = 1000000

  # Poisson ratio of wall
  set poisson ratio wall                                = 0.3

  # Coefficient of restitution of wall
  set restitution coefficient wall                      = 0.95

  # Coefficient of friction of wall
  set friction coefficient wall                         = 0.05

  # Coefficient of rolling friction of wall
  set rolling friction wall                             = 0.1
 end

* The ``gx``, ``gy``, and ``gz`` parameters define the gravitational acceleration in `x`, `y`, and `z` directions.

* The ``number of particle types`` parameter specifies the number of particle types in a simulation. Particles with different sizes, size distributions, and physical properties have to be defined as separate particle types.

* For each particle type, we have to define a separate subsection (for instance, ``subsection particle type 0``) to specify its physical properties.

.. note::
    If the particles in a simulation are monodispersed and have the same physical properties, the ``number of particle types`` should be equal to zero. For polydispersed systems, the ``number of particle types`` is selected equal to the number of particles types in the simulation. For each particle type, a separate subsection ``particle type n`` should be defined (n starts from zero to ``number of particle types`` - 1) which contains all the physical properties related to that particle type.

* The ``size distribution type`` parameter specifies the size distribution for each particle type. The acceptable choices are ``uniform`` and ``normal`` distributions.

.. note::
    For each particle type, two ``size distribution type``s can be defined: ``uniform`` and ``normal``. In ``uniform`` size distribution, the diameter of the particles is constant, while in ``normal`` size distribution, the particle diameters are sampled from a normal distribution with an average of ``average diameter`` and standard deviation of ``standard deviation``.

* The ``diameter`` parameter defines the diameter of the particles in a ``uniform`` distribution.

* For a ``normal`` distribution, we need to define ``average diameter`` and ``standard deviation`` parameters.

* The ``number`` parameter defines the number of particles for each type.

* The ``density particles`` defines the density of particles for each type.

* The ``young modulus particles`` defines the Young's modulus for particles in each type.

* The ``poisson ratio particles`` defines the Poisson's ratio for particles in each type.

* The ``restitution coefficient particles`` defines the restitution coefficient for particles in each type.

* The ``friction coefficient particles`` defines the friction coefficient for particles in each type.

* The ``rolling friction particles`` defines the rolling friction coefficient of particles for each type.

* The ``young modulus wall`` defines the Young's modulus of the walls.

* The ``poisson ratio wall`` defines the Poisson's ratio of the walls.

* The ``restitution coefficient wall`` defines the restitution coefficient of the walls.

* The ``friction coefficient wall`` defines the friction coefficient of the walls.

* The ``rolling friction wall`` defines the rolling friction coefficient of the walls.

