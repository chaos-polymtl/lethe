==============
Insertion Info
==============

In this subsection, insertion methods which are ``volume``, ``plane``, ``list`` and ``file`` are defined.

.. note::
    Insertion in Lethe starts inserting particles from type 0 and proceeds to the next type when all the particles from the previous type are inserted.

.. code-block:: text

  subsection insertion info
    # Choices are volume|plane|list|file
    set insertion method                               = volume

    # Every method
    set insertion frequency                            = 20000

    # If method = volume
    set inserted number of particles at each time step = 100

    set insertion box points coordinates               = -0.05, -0.05, -0.03 : 0.05, 0.05, 0.07

    set insertion order of direction                   = 0, 1, 2

    set initial velocity                               = 0.0, 0.0, 0.0
    set initial angular velocity                       = 0.0, 0.0, 0.0

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

    # If method = file
    set insertion file name                            = particles.input
  end

The ``insertion method`` parameter chooses the type of insertion. Acceptable choices are ``volume``, ``plane``, ``list`` and ``file``. Different insertion method can share the same parameter.

* The ``insertion frequency`` parameter defines the frequency of insertion. For example, if the ``insertion frequency`` is set equal to 10000, the iterations 1, 10001, 20001, ... will be defined as insertion steps.  The ``insertion frequency`` should be selected adequately depending on the insertion method. For ``volume`` it should be large, so that the inserted particles at the previous insertion step have enough time to leave the insertion box for the next insertion step, otherwise large overlap may occur which leads to a large velocity of particles. For the ``plane`` method, it should be small so that particles are being inserted as soon as a cell is empty.

* The ``insertion maximum offset`` and ``insertion prn seed`` parameters defines the random offset values to the initial positions of particles during a ``volume`` and ``plane`` insertion. The ``insertion maximum offset`` parameter defines the maximum value for an offset. The ``insertion prn seed`` parameter defines the pseudo-random number (PRN) with which offset values are getting generated.

-------
Volume
-------
The ``volume`` insertion method uses an insertion box where particles will be inserted. The insertion locations of particles are randomly selected if the ``insertion maximum offset`` is not equal to zero, otherwise, the particles will perfectly aligns with the x, y and z directions.

* The ``inserted number of particles at each time step`` defines the desired number of particles to be inserted at each insertion step. If the insertion box is not adequately large to insert ``inserted number of particles at each time step`` particles with the defined arrangement (initial distance between the inserted particles), Lethe prints a warning and inserts the maximum number of particles that fit inside the insertion box at each insertion step.

* The ``insertion box points coordinates`` parameter defines the insertion box dimensions using two points: ``x1, y1, z1 : x2, y2, z2``. It is the same principle has what is being done for the `CFD <https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/mesh.html>`_ triangulation.

.. note::
    We recommend that the defined insertion box have at least a distance of :math:`{d^{max}_p}` (maximum diameter of particles) from the triangulation boundaries. Otherwise, particles may have an overlap with the triangulation walls in the insertion.

* The ``insertion order of direction`` parameter defines the directions of insertion. For example, if the parameter is equal to ``0, 1, 2``, the particles are inserted in priority in the x, in y, and then in z directions. This is the default configuration. This is useful to specify the insertion directions to cover a specific area of the insertion box with the first and second direction parameters.

* The ``initial velocity`` determine the initial translational velocity (in :math:`\frac{m}{s}`) at which particles are inserted in the x, y, and z directions.

* The ``initial angular velocity`` determine the initial rotational velocity (in :math:`\frac{rad}{s}`) at which particles are inserted in the x, y, and z directions.

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
The ``plane`` insertion method inserts particles at the centroid of insertion cells. These cells are defined as intersected by a mathematical plane. This plane is define by an ``insertion plane point`` and an ``insertion plane normal vector``. A cell is considered as intersected by the plane if at least one of its vertex is on each side of the plane of if at least one of its vertex is directly on the plane (the normal distance between the vertex and the plane is zero). At each insertion step, a particle will be inserted in a insertion cell if that cell is empty (no particle is present inside it). This guarantee the absence of big overlap with the particles already inserted. This method of inserting is useful when dealing with a domain dense with particles.

* The ``insert plane point`` defines the point coordinates for the plane. Each component of this parameter represent the x, y and z directions, respectively.

* The ``insertion plane normal vector`` define the normal vector component for the plane. of the  Each component of the parameter represent the x, y and z directions, respectively.

--------------------
List
--------------------
The ``list`` insertion method insert particles at precis coordinates with specific velocities (translational and angular) and diameters.  This method is preferred for small number of particles.

* The ``list x``, ``list y`` and ``list z`` define the coordinates of every particles in the x, y and z directions, respectively. For example, if you want to insert particles at two locations, ``(0.,0.,0.) and (1.,2.,3.)`` , the list parameters should look like this :

.. code-block:: text

    set list x = 0., 1.
    set list y = 0., 2.
    set list z = 0., 3.

* The ``list velocity x``, ``list velocity y``, ``list velocity z``, ``list omega x``, ``list omega y``, ``list omega z`` and ``list diameters`` define the initial translational velocities, the initial angular velocities and diameters of each particles respectively following the same logic as the insertion coordinates.

---------------------
File
---------------------
The ``file`` insertion method insert particles in a similar way to the ``list`` insertion method. The main difference between these two methods is the option to use an external file provided by the ``insertion file name`` parameter. This parameter is set at ``particles.input`` by default. This file has to follow this structure:

.. code-block:: text

    p_x; p_y; p_z; v_x; v_y; v_z; w_x; w_y; w_z; diameters; fem_force_x; fem_force_y; fem_force_z; fem_torque_x; fem_torque_y; fem_torque_z;
    0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;       0.2;           0;           0;           0;            0;            0;            0;
    1.0; 2.0; 3.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;       0.2;           0;           0;           0;            0;            0;            0;

Each line is associated with a particle and its properties. The ``fem_force`` and ``fem_torque`` properties are only used in the CFD-DEM solver, but must be specified in all cases. The main advantage of using the ``file`` method over the ``list`` method is that the number of inserted particles is not limited to the maximum number of characters on a single line of parameter files. To generate an insertion file, particle positions and properties can be generated manually or with any script. An other option is to use the python code ``extract-particles-properties-from-vtu.py`` in ``lethe/contrib/preprocessing/`` directory. This code extracts particle properties from the last vtu file from a given simulation.

.. note::
    The ``file`` insertion combine with the ``extract-particles-properties-from-vtu.py`` python code can be a useful tool. The loading of particles and the rest of the simulation can be performed in two different triangulations, witch is not the case of the the restart feature. This means that the loading triangulation can have smaller cells and a bigger domain to allow for the use of larger insertion boxes. Then, particles properties can be extracted and the remainder of the simulation can be performed in the appropriate triangulation.

.. warning::
    The critical Rayleigh time step is computed from the parameters in the ``particle type`` subsections, not the ``insertion info`` subsection. It is the user's responsibility to fill the ``particle type`` subsections correctly according to the diameter values stored in the insertion input file, otherwise Rayleigh time percentage displayed at the start of every DEM simulation may not be accurate.