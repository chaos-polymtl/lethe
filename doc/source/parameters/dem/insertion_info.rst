==============
Insertion Info
==============

In this subsection, insertion methods which are ``volume``, ``plane`` and ``list`` are defined.

.. note::
    Insertion in Lethe starts inserting particles from type 0 and proceeds to the next type when all the particles from the previous type are inserted.


.. code-block:: text

  subsection insertion info
    # Choices are volume|plane|list
    set insertion method                               = volume

    # Every method
    set insertion frequency                            = 20000

    # If method = volume
    set inserted number of particles at each time step = 100

    set insertion box minimum x                        = -0.05
    set insertion box minimum y                        = -0.05
    set insertion box minimum z                        = -0.03
    set insertion box maximum x                        = 0.05
    set insertion box maximum y                        = 0.05
    set insertion box maximum z                        = 0.07

    set insertion first direction                      = 0
    set insertion second direction                     = 1
    set insertion third direction                      = 2

    set velocity x                                     = 0.0
    set velocity y                                     = 0.0
    set velocity z                                     = 0.0
    set omega x                                        = 0.0
    set omega y                                        = 0.0
    set omega z                                        = 0.0

    set insertion distance threshold                   = 2

    # If method = plane
    set insertion method                               = plane
    set insertion plane point                          = 0, 0, 0
    set insertion plane normal vector                  = 1, 0, 0

    # If method = volume or plane
    set insertion maximum offset                       = 0.75
    set insertion prn seed                             = 19

    # If method = list
    set list x                                         = 0.
    set list y                                         = 0.
    set list z                                         = 0.
    set list velocity x                                = 0.
    set list velocity y                                = 0.
    set list velocity z                                = 0.
    set list omega x                                   = 0.
    set list omega y                                   = 0.
    set list omega z                                   = 0.
    set list diameters                                 = 0.

  end

The ``insertion method`` parameter chooses the type of insertion. Acceptable choices are ``volume``, ``plane`` and ``list``. Different insertion method can share the same parameter.

* The ``insertion frequency`` parameter defines the frequency of insertion. For example, if the ``insertion frequency`` is set equal to 10000, the iterations 1, 10001, 20001, ... will be defined as insertion steps.  The ``insertion frequency`` should be selected adequately depending on the insertion method. For ``volume`` it should be large, so that the inserted particles at the previous insertion step have enough time to leave the insertion box for the next insertion step, otherwise large overlap may occur which leads to a large velocity of particles. For the ``plane`` method, it should be small so that particles are being inserted as soon as a cell is empty.

* The ``insertion maximum offset`` and ``insertion prn seed`` parameters defines the random offset values to the initial positions of particles during a ``volume`` and ``plane`` insertion. The ``insertion maximum offset`` parameter defines the maximum value for an offset. The ``insertion prn seed`` parameter defines the pseudo-random number (PRN) with which offset values are getting generated.

-------
Volume
-------
The volume insertion method uses an insertion box where particles will be inserted. The insertion locations of particles are randomly selected if the ``insertion maximum offset`` is not equal to zero, otherwise, the particles will perfectly aligns with the x, y and z directions.

* The ``inserted number of particles at each time step`` defines the desired number of particles to be inserted at each insertion step. If the insertion box is not adequately large to insert ``inserted number of particles at each time step`` particles with the defined arrangement (initial distance between the inserted particles), Lethe prints a warning and inserts the maximum number of particles that fit inside the insertion box at each insertion step.

* The ``insertion box minimum x``, ``insertion box minimum y``, ``insertion box minimum z``, ``insertion box maximum x``, ``insertion box maximum y``, ``insertion box maximum z`` parameters define the insertion box dimensions.

.. note::
    We recommend that the defined insertion box have at least a distance of :math:`{d^{max}_p}` (maximum diameter of particles) from the triangulation boundaries. Otherwise, particles may have an overlap with the triangulation walls in the insertion.

* The ``insertion first direction``, ``insertion second direction``, and ``insertion third direction`` parameters define the directions of insertion. For example, if ``insertion first direction`` = 0, ``insertion second direction`` = 1, and ``insertion third direction`` = 2, the particles are inserted in priority in the x, in y, and then in z directions. This is the default configuration. This is useful to specify the insertion directions to cover a specific area of the insertion box with the first and second direction parameters.

* The ``velocity x``, ``velocity y``, and ``velocity z`` determine the initial translational velocity (in :math:`\frac{m}{s}`) at which particles are inserted in the x, y, and z directions, respectively.

* The ``omega x``, ``omega y``, and ``omega z`` determine the initial rotational velocity (in :math:`\frac{rad}{s}`) at which particles are inserted in the x, y, and z directions, respectively. 

.. note:: 
    Since the ``insertion info`` subsection is valid for all particle types, by using ``velocity x``, ``velocity y``, ``velocity z``, ``omega x``, ``omega y``, or ``omega z``, the given condition is applied to all particles, indistinctively.

* The ``insertion distance threshold`` parameter determines the initial distance between the particles in the insertion box. As a result, it must be larger than 1 to avoid any initial collision between the inserted particles.

The distance between the inserted particles is equal to:

.. math::
    D_i=(\epsilon + \psi)  d^{max}_p

Where, :math:`{\epsilon}`, :math:`{\psi}`, and :math:`{d^{max}_p}` denote ``insertion distance threshold``, a generated random number (in the range of 0-``insertion maximum offset``, and from the seed of ``insertion prn seed``), and maximum particle diameter.
 
.. note::
    ``insertion distance threshold`` should also be compatible with the ``insertion maximum offset``; especially if the ``insertion maximum offset`` is large, a large value should be defined for ``insertion distance threshold``. Generally, we recommend users to use a value in the range of 1.3-2 (depending on the value of ``insertion maximum offset``) for the ``insertion distance threshold``.

--------------------
Plane
--------------------
The Plane insertion method inserts particles at the centroid of insertion cells. These cells are defined as intersected by a mathematical plane. This plane is define by an ``insertion plane point`` and an ``insertion plane normal vector``. A cell is considered as intersected by the plane if at least one of its vertex is on each side of the plane of if at least one of its vertex is directly on the plane (the normal distance between the vertex and the plane is zero). At each insertion step, a particle will be inserted in a insertion cell if that cell is empty (no particle is present inside it). This guarantee the absence of big overlap with the particles already inserted. This method of inserting is useful when dealing with a domain dense with particles.

* The ``insert plane point`` defines the point coordinates for the plane. Each component of this parameter represent the x, y and z directions, respectively.

* The ``insertion plane normal vector`` define the normal vector component for the plane. of the  Each component of the parameter represent the x, y and z directions, respectively.

--------------------
List
--------------------
The List insertion method insert particles at precis coordinates with specific velocities (translational and angular) and diameters.  This method is preferred for small number of particles.

* The ``list x``, ``list y`` and ``list z`` define the coordinates of every particles in the x, y and z directions, respectively. For example, if you want to insert particles at two locations, ``(0.,0.,0.) and (1.,2.,3.)`` , the list parameters should look like this :

.. code-block:: text

    set list x = 0., 1.
    set list y = 0., 2.
    set list z = 0., 3.

* The ``list velocity x``, ``list velocity y``, ``list velocity z``, ``list omega x``, ``list omega y``, ``list omega z`` and ``list diameters`` define the initial translational velocities, the initial angular velocities and diameters of each particles respectively following the same logic as the insertion coordinates.
