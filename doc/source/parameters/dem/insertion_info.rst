==============
Insertion Info
==============

In this subsection, insertion methods ``volume``, ``plane``, ``list`` and ``file`` are defined.

.. note::
    Insertion in Lethe starts inserting particles from type 0 and proceeds to the next type when all the particles from the previous type are inserted.

.. code-block:: text

  subsection insertion info
    # Choices are volume|plane|list|file
    set insertion method                               = volume

    # Every method
    set insertion frequency                            = 1

    # If method = volume
    set inserted number of particles at each time step = 1
    set insertion box points coordinates               = 0., 0., 0. : 1., 1., 1.
    set insertion insertion direction sequence         = 0, 1, 2
    set insertion distance threshold                   = 1.

    # If method = plane
    set insertion method                               = plane
    set insertion plane point                          = 0, 0, 0
    set insertion plane normal vector                  = 1, 0, 0
    set insertion plane threshold distance             = 0.

    # If method = volume or plane
    set initial velocity                               = 0.0, 0.0, 0.0
    set initial angular velocity                       = 0.0, 0.0, 0.0
    set insertion maximum offset                       = 1.
    set insertion prn seed                             = 1

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
    set list of input files                            = particles.input

    # Box removal
    set remove particles                               = false
    set removal box points coordinates                 = 0., 0., 0.: 1., 1., 1.
  end

The ``insertion method`` parameter defines insertion type, where choices are: ``volume``, ``plane``, ``list`` and ``file``. Different insertion methods can share the same parameter. In order to simplify the explanation on those parameters, they are introduced in each method with their specificities.

-------
Volume
-------
The ``volume`` insertion method uses an insertion box where particles will be inserted. The box can inserts particles in a structured manner or in a random manner according to the  ``insertion maximum offset`` and ``insertion prn seed`` settings.

* ``insertion frequency`` defines the frequency of the insertion of particles in the box. For example, if the ``insertion frequency`` is set to 10000, the iterations 1, 10001, 20001, ... will be defined as insertion steps. The frequency should be large, so that the inserted particles at the previous insertion step have enough time to leave the insertion box for the next insertion step, otherwise large overlap may occur which leads to a large velocity of particles.

* ``insertion box points coordinates`` defines the insertion box dimensions using two points: ``x1, y1, z1 : x2, y2, z2``. It is based on the same principle as for the rectangular mesh `CFD <../../parameters/cfd/mesh.html>`_ triangulation.

.. note::
    The insertion box may fit tightly to the triangulation, but the following condition must be met:

    .. math::
        \psi < \frac{\epsilon - 1}{2}

    Where, :math:`{\psi}`, :math:`{\epsilon}`, and :math:`{d^{max}_p}` denote a generated random number (in the range of 0-``insertion maximum offset``, and from the seed of ``insertion prn seed``), the ``insertion distance threshold``, and the maximum particle diameter, respectively.

    Otherwise, particles may have an overlap with the triangulation walls during the insertion. The minimal insertion box coordinates should be adjusted so that:

    .. math::
        (x_1, y_1, z_1) > (x_{min}, y_{min}, z_{min}) + \left(\frac{1-\epsilon}{2} + \psi\right) d^{max}_p

* ``inserted number of particles at each time step`` defines the desired number of particles to be inserted at each insertion step. If the insertion box is not adequately large to insert this number of particles with the defined arrangement (initial distance between the inserted particles), the maximum number of particles that fits inside the insertion box at each insertion step are inserted and a warning is displayed.

* ``insertion maximum offset`` defines the maximum value of the offset in relation to the structured discrete positions in the box. If the offset is equal to 0.0, the particles will perfectly aligns along the x, y and z directions. Otherwise, the particle insertion locations are randomly selected within the offset of this positions.

* ``insertion prn seed`` seeds the pseudo-random number (PRN) generator. It defines the starting value from which the offset values are generated.

* ``insertion distance threshold`` determines the initial distance between the particles in the insertion box. It must be larger than 1 to avoid any initial collision between the inserted particles.
  The distance between the inserted particles is equal to:

  .. math::
      D_i=(\epsilon + \psi)  d^{max}_p

.. note::
    ``insertion distance threshold`` should also be compatible with the ``insertion maximum offset``. Inserted particles will not overlap if:
    :math:`\epsilon < \psi + 1` See note on the ``insertion box points coordinates`` parameter.

    Generally, we recommend users to use a threshold in the range of 1.3-2.0, depending on the value of offset.

* ``insertion direction sequence`` defines the sequence of directions of insertion in the box. For example, if the parameter is equal to ``0, 1, 2``, the particles are inserted in priority in the x, in y, and then in z directions. This is the default configuration. This is useful to specify the insertion directions to cover a specific area of the insertion box with the first and second direction parameters.

* ``initial velocity`` determine the initial translational velocity (in :math:`\frac{m}{s}`) at which particles are inserted in the x, y, and z directions.

* ``initial angular velocity`` determine the initial rotational velocity (in :math:`\frac{rad}{s}`) at which particles are inserted in the x, y, and z directions.

.. note::
    Since the ``insertion info`` subsection is valid for all particle types, by using ``initial velocity`` or ``initial angular velocity`` the given conditions are applied to all particles, regardless of the type.

--------------------
Plane
--------------------
The ``plane`` insertion method inserts particles at the centroid of insertion cells. These cells are defined as intersected by a mathematical plane. This plane is defined by an ``insertion plane point`` and an ``insertion plane normal vector``. A cell is considered as intersected by the plane if at least one of its vertex is on each side of the plane or if at least one of its vertex is directly on the plane (the normal distance between the vertex and the plane is zero). At each insertion step, a particle will be inserted in a insertion cell if that cell is empty (no particle is present inside it). This guarantee the absence of big overlap with the particles already inserted. This method of inserting is useful when dealing with a domain dense with particles.

* ``insertion frequency`` defines the frequency of the check for particle insertion. The insertion method will check if the cell in empty, and will only insert a particle if so. The frequency should be small so that particles are being inserted as soon as a cell is empty.

* ``insertion maximum offset`` defines the maximum value of the offset in relation to centroid of the cell. The insertion locations of particles are randomly selected if the offset is not equal to zero, otherwise, the particles will be inserted at the centroid.

* ``insertion prn seed`` seeds the pseudo-random number (PRN) generator. It defines the starting value from which the offset values are generated.

* ``insert plane point`` defines the point coordinates for the plane. Each component of this parameter represent the x, y and z directions, respectively.

* ``insertion plane normal vector`` define the normal vector component for the plane. Each component of the parameter represent the x, y and z directions, respectively.

* ``initial velocity`` determine the initial translational velocity (in :math:`\frac{m}{s}`) at which particles are inserted in the x, y, and z directions.

* ``initial angular velocity`` determine the initial rotational velocity (in :math:`\frac{rad}{s}`) at which particles are inserted in the x, y, and z directions.

--------------------
List
--------------------
The ``list`` insertion method insert particles at precis coordinates with specific velocities (translational and angular) and diameters.  This method is preferred for small number of particles.

* ``insertion frequency`` defines the frequency of the insertion of particles based on the list. If the list contains 3 coordinates, 3 new particles will be inserted at the same positions at each insertion step.

* ``list x``, ``list y``, and ``list z``: define the coordinates of every particles in the x, y and z directions, respectively. For example, if you want to insert particles at two locations, ``(0.,0.,0.) and (1.,2.,3.)`` , the list parameters should look like this :

.. code-block:: text

    set list x = 0., 1.
    set list y = 0., 2.
    set list z = 0., 3.

* ``list velocity x``, ``list velocity y``, and ``list velocity z`` define the initial translational velocities of each particles respectively following the same logic as the insertion coordinates.

* ``list omega x``, ``list omega y``, and ``list omega z`` define the initial angular velocities of each particles respectively following the same logic as the insertion coordinates.

* ``list diameters`` defines the diameters of each particles respectively following the same logic as the insertion coordinates.

---------------------
File
---------------------
The ``file`` insertion method inserts particles in a similar way to the ``list`` insertion method. The main difference between these two methods is the use of external files provided by the ``list of input files`` parameter. A single file or a list of files may be specified. At each insertion time step, a different file will be used. If the end of the file list is reached and there are still particles to be inserted, the list returns to the first file. An insertion file must follow this structure:

.. code-block:: text

    p_x; p_y; p_z; v_x; v_y; v_z; w_x; w_y; w_z; diameters;
    0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;       0.2;
    1.0; 2.0; 3.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;       0.2;

Each line is associated with a particle and its properties. The main advantage of using the ``file`` method over the ``list`` method is that the number of inserted particles is not limited to the maximum number of characters on a single line of parameter files. To generate an insertion file, particle positions and properties can be generated manually or with any script. An other option is to use the python code ``extract-particles-properties-from-vtu.py`` in ``lethe/contrib/preprocessing/`` directory. This code extracts particle properties from the last vtu file from a given simulation.

* ``insertion frequency`` defines the frequency of the insertion of particles based on the list in the file(s)

* ``list of input files`` defines the list of files to be used for the insertion. The default value is ``particles.input``.

.. note::
    The ``file`` insertion combine with the ``extract-particles-properties-from-vtu.py`` python code can be a useful tool. The loading of particles and the rest of the simulation can be performed in two different triangulations, witch is not the case of the the restart feature. This means that the loading triangulation can have smaller cells and a bigger domain to allow for the use of larger insertion boxes. Then, particles properties can be extracted and the remainder of the simulation can be performed in the appropriate triangulation.

.. warning::
    The critical Rayleigh time step is computed from the parameters in the ``particle type`` subsections, not the ``insertion info`` subsection. It is the user's responsibility to fill the ``particle type`` subsections correctly according to the diameter values stored in the insertion input file, otherwise Rayleigh time percentage displayed at the start of every DEM simulation may not be accurate.

--------------------
Removal
--------------------
With all insertion methods, it is possible to define a removal box where particles will be removed from the triangulation just before the insertion of new particles.

* ``remove particles`` enables (true) or disables (false) the particle removal.

* ``removal box points coordinates`` defines a removal box where particles will be removed. It uses the same principle as the insertion box.