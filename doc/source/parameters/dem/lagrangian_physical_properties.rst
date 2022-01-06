Lagrangian Physical Properties
----------------------------------
In this subsection, the physical properties of the simulation are defined. These properties include gravitational acceleration (``gx``, ``gy`` and ``gz``), ``number of particle types``, and for each type, particle ``diameter``, particle ``density``, ``Young's modulus`` of particle and wall, ``Poisson ratio`` of particle and wall, ``restitution coefficient`` of particle and wall, ``friction coefficient`` of particle and wall and ``rolling friction coefficient`` of particle and wall.

.. note::
    If the particles in the simulation are monodispersed, the ``number of particle types`` should be equal to zero. For polydispersed systems, the ``number of particle types`` is selected equal to the number of particles types in the simulation. For each particle type, a separate subsection ``particle type n`` should be defined (n starts from zero to ``number of particle types`` - 1) which contains all the physical properties related to that particle type.

.. note::
    For each particle type, two ``size distribution type``s can be defined: ``uniform`` and ``normal``. In ``uniform`` size distribution, the diameter of the particles is constant, while in ``normal`` size distribution, the particle diameters are sampled from a normal distribution with an average of ``average diameter`` and standard deviation of ``standard deviation``.

.. code-block:: text

 subsection model parameters
  # Gravitational acceleration in x direction
  set gx            		 			 = 0.0

  # Gravitational acceleration in y direction
  set gy            		 			 = 0.0

  # Gravitational acceleration in z direction
  set gz            		 			 = -9.81

  # Number of particle types
  set number of particle types	                         = 1

  # Entering particle type 0
    	subsection particle type 0

  # Size distribution of particle type 0
		set size distribution type		= uniform

  # Particle diameter
                set diameter            	 	= 0.005

  # Number of particles in type 0
		set number				= 132300

  # Particle density
                set density particles           	= 2000

  # Young's modulus of particle
                set young modulus particles         	= 1000000

  # Poisson ratio of particle
                set poisson ratio particles            	= 0.3

  # Coefficient of restitution of particle
                set restitution coefficient particles   = 0.95

  # Coefficient of friction of particle
                set friction coefficient particles      = 0.05

  # Coefficient of rolling friction of particle
                set rolling friction particles          = 0.1
	end

  # Young's modulus of wall
  set young modulus wall            			 = 1000000

  # Poisson ratio of wall
  set poisson ratio wall            			 = 0.3

  # Coefficient of restitution of wall
  set restitution coefficient wall           		 = 0.95

  # Coefficient of friction of wall
  set friction coefficient wall         		 = 0.05

  # Coefficient of rolling friction of wall
  set rolling friction wall         	      	  	 = 0.1
 end