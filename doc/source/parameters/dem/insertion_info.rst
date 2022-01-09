Insertion Info
-------------------
In Lethe, particles are inserted using the Insertion class. We need an insertion box (which is defined by its minimum and maximum values in `x`, `y`, and `z` directions), ``insertion frequency``, the inserted number of particles at each insertion step, and other parameters to specify the relative positions of the particles during insertion.

Particle insertion information is defined in this section. This information includes ``insertion method``, ``inserted number of particles at each insertion step``, ``insertion frequency``, dimensions of the insertion box, approximate initial distance between particles at insertion steps, and information about the randomness of initial positions of particles. If the insertion method is ``uniform``, the particles are uniformly (without randomness in initial location) inserted. On the other hand, in ``non_uniform`` insertion, the initial locations of particles are chosen randomly. Using ``non_uniform`` insertion, insertion ``random number range`` and ``insertion random number seed`` are also required in the parameter handler file.

.. note::
    Insertion in Lethe-DEM starts inserting particles from type 0 and proceeds to the next type when all the particles from the previous type are inserted.


.. code-block:: text

 subsection insertion info
  # Insertion method
  # Choices are uniform|non_uniform
  set insertion method				        = non_uniform

  # Number of inserted particles at each insertion step. This value may change automatically if the insertion box is not adequately large to handle all the inserted particles
  set inserted number of particles at each time step    = 100

  # Insertion frequency
  set insertion frequency                       = 20000

  # Insertion box dimensions (xmin, xmax, ymin, ymax, zmin, zmax)
  set insertion box minimum x                   = -0.05
  set insertion box minimum y                   = -0.05
  set insertion box minimum z                   = -0.03
  set insertion box maximum x                   = 0.05
  set insertion box maximum y                   = 0.05
  set insertion box maximum z                   = 0.07

  # insertion distance threshold controls the distance between the center of inserted particles (the distance is: [distance threshold] * [diameter of particles]). The distance is modified by a random number if non_uniform insertion is chosen
  set insertion distance threshold              = 2

  # Random number range and seed for non_uniform insertion
  set insertion random number range             = 0.75
  set insertion random number seed              = 19
 end

* The ``insertion method`` parameter chooses the type of insertion. Acceptable choices are ``uniform`` and ``non_uniform``.

.. note::
    In ``uniform`` insertion, the insertion locations of the particles inside the insertion box are tidy. In ``non_uniform`` insertion, however, the insertion locations of particles are randomly selected.


* The ``inserted number of particles at each time step`` defines the desired number of particles to be inserted at each insertion step.

.. note::
    If the insertion box is not adequately large to insert ``inserted number of particles at each time step`` particles with the defined arrangement (initial distance between the inserted particles), Lethe prints a warning and inserts the maximum number of particles that fit inside the insertion box at each insertion step.

* The ``insertion frequency`` parameter defines the frequency of insertion. For example, if the ``insertion frequency`` is set equal to 10000, the iterations 0, 10000, 20000, ... will be defined as insertion steps.

.. note::
    The ``insertion frequency`` should be selected adequately large, so that the inserted particles at the previous insertion step have enough time to leave the insertion box, and the insertion box becomes empty for the next insertion. If ``insertion frequency`` does not meet this criterion, particles may have large overlaps during the insertion steps which leads to a large velocity of particles. Some particles may leave the simulation domain in this case.

* The ``insertion box minimum x``, ``insertion box minimum y``, ``insertion box minimum z``, ``insertion box maximum x``, ``insertion box maximum y``, ``insertion box maximum z`` parameters define the insertion box dimensions.

.. note::
    We recommend that the defined insertion box have at least a distance of :math:`{d^{max}_p}` (maximum diameter of particles) from the triangulation boundaries. Otherwise, particles may have an overlap with the triangulation walls in the insertion.

* The ``insertion distance threshold`` parameter determines the initial distance between the particles in the insertion. As a result, it must be larger than 1 to avoid any initial collision between the inserted particles.

* The ``random number range`` and ``insertion random number seed`` parameters determine the random added values to the positions of particles during a ``non_uniform`` insertion. ``random number range`` defines the maximum value for the random displacement in the ``non_uniform`` insertion locations. ``insertion random number seed`` is the seed for the random number generator.

The distance between the inserted particles is equal to:

.. math::
    D_i=\epsilon * d^{max}_p

in an ``uniform`` insertion, and

.. math::
    D_i=(\epsilon + \psi)  d^{max}_p

in a ``non_uniform`` insertion. :math:`{\epsilon}`, :math:`{\psi}`, and :math:`{d^{max}_p}` denote ``insertion distance threshold``, a generated random number (in the range of 0-``random number range``, and from the seed of ``insertion random number seed``), and maximum particle diameter.
 
.. note::
     ``insertion distance threshold`` should also be compatible with the ``random number range``; especially if the ``random number range`` is large, a large value should be defined for ``insertion distance threshold``. Generally, we recommend users to use a value in the range of 1.3-2 (depending on the value of ``random number range``) for the ``insertion distance threshold``.

