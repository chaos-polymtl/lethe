==============================
Lagrangian Physical Properties
==============================

In this subsection, gravitational acceleration, and the physical properties of the particles and walls are defined. These properties include ``number of particle types``, and for each type, particle ``diameter``, particle ``density``, ``Young's modulus`` of particle and wall, ``Poisson ratio`` of particle and wall, ``restitution coefficient`` of particle and wall, ``friction coefficient`` of particle and wall and ``rolling friction coefficient`` of particle and wall.

.. code-block:: text

  subsection lagrangian physical properties
    # Gravitational acceleration in x direction
    set gx                       = 0.0

    # Gravitational acceleration in y direction
    set gy                       = 0.0

    # Gravitational acceleration in z direction
    set gz                       = -9.81

    # Number of particle types
    set number of particle types = 1

    # Entering particle type 0
    subsection particle type 0

      # Choices are uniform, normal or custom
      set size distribution type            = uniform

      # If distribution type = uniform or normal
      set diameter                          = 0.001

      # If distribution type = custom
      set custom diameters                  = 0.001 , 0.0005
      set custom probabilities              = 0.6   , 0.4

      # If distribution type = normal
      set standard deviation                = 0.0

      # For every distribution types
      set number of particles               = 0
      set density particles                 = 1000
      set young modulus particles           = 1000000
      set poisson ratio particles           = 0.3
      set restitution coefficient particles = 0.1
      set friction coefficient particles    = 0.1
      set rolling friction particles        = 0.1
      set surface energy particles          = 0.0

    end

    # Wall properties
    set young modulus wall           = 1000000
    set poisson ratio wall           = 0.3
    set restitution coefficient wall = 0.1
    set friction coefficient wall    = 0.1
    set rolling friction wall        = 0.1
    set surface energy wall          = 0.0
  end

* The ``gx``, ``gy``, and ``gz`` parameters define the gravitational acceleration in `x`, `y`, and `z` directions.

* The ``number of particle types`` parameter specifies the number of particle types in a simulation. Particles with different sizes, size distributions, and physical properties have to be defined as separate particle types.

* For each particle type, we have to define a separate subsection (for instance, ``subsection particle type 0``) to specify its physical properties.

.. note::
    If the particles in a simulation are monodispersed and have the same physical properties, the ``number of particle types`` should be equal to zero. For polydispersed systems, the ``number of particle types`` is selected equal to the number of particles types in the simulation. For each particle type, a separate subsection ``particle type n`` should be defined (n starts from zero to ``number of particle types`` - 1) which contains all the physical properties related to that particle type.

* The ``size distribution type`` parameter specifies the size distribution for each particle type. For each particle type, three ``size distribution type`` can be defined: ``uniform``, ``normal`` and ``custom``.

  - For the ``uniform`` size distribution, the diameter of the particles is constant.
  - For the ``normal`` size distribution, the particle diameters are sampled from a normal distribution with an average diameter and a standard deviation.
  - For the ``custom`` size distribution, particle diameters are sampled from a list of diameters with a corresponding list of probabilities.

.. note::
    In the ``custom`` size distribution, the probability values are based on the volume fraction taken by all the particles of the associated diameter, not to the total number of particles. For example, if a probability is equal to ``0.5`` , this means that half of the total volume of inserted particles will be occupied by particle with the associated diameter value.

* The ``diameter`` parameter defines the diameter of the particles in a ``uniform`` distribution. In the case of a ``normal`` distribution, this parameter indicates the average diameter.

* For a ``normal`` distribution, the ``standard deviation`` parameter should be defined to indicate the standard deviation on the particle size distribution.

* For a ``custom`` distribution, the ``custom diameters`` parameter defines the different diameter values used when generating particles. The ``custom probabilities`` parameter defines the probabilities corresponding to each diameter value previously declared based on volume fraction, not on the number of particles. Both list must have the same length.

* The ``number of particles`` parameter defines the number of particles for each type.

* The ``density particles`` defines the density of particles for each type.

* The ``young modulus particles`` defines the Young's modulus for particles in each type.

* The ``poisson ratio particles`` defines the Poisson's ratio for particles in each type.

* The ``restitution coefficient particles`` defines the restitution coefficient for particles in each type.

* The ``friction coefficient particles`` defines the friction coefficient for particles in each type.

* The ``rolling friction particles`` defines the rolling friction coefficient of particles for each type.

* The ``surface energy particles`` defines the surface energy of particles for each type. This parameter is used with the JKR force model.

* The ``young modulus wall`` defines the Young's modulus of the walls.

* The ``poisson ratio wall`` defines the Poisson's ratio of the walls.

* The ``restitution coefficient wall`` defines the restitution coefficient of the walls.

* The ``friction coefficient wall`` defines the friction coefficient of the walls.

* The ``rolling friction wall`` defines the rolling friction coefficient of the walls.

* The ``surface energy wall`` defines the surface energy of the walls. This parameter is used with the JKR force model.

