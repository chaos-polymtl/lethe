Post-processing
-------------------
Lethe has built in post-processing capabilities. The post-processing subsection is according to the following example:

.. code-block:: text

 subsection post-processing
  # Enable the use of calculate granular temperature and calculate particles' average velocity
  set Lagrangian post-processing = false

  # Enable writing grid VTU files
  set write grid = true

  # Enable calculate particles' average velocity
  set calculate particles average velocity = false

  # Enable calculate granular temperature
  set calculate granular temperature = false

  # Choose time step to start calculating particles average velocity and/or granular temperature
  set initial step = 0

  # Choose time step to stop calculating particles average velocity and/or granular temperature
  set end step = 0

  # Choose output frequency of particles average velocity and/or granular temperature calculation
  set output frequency = 100000

  # Name of output file with particle's velocity
  set particles velocity output name = particles_velocity

  # Name of output file with granular temperature
  set granular temperature output name = granular_temperature

 end

.. note::
 The parameters default values are the ones presented above.

* The ``Lagrangian post processing`` enables the use of ``calculate particles average velocity``, ``calculate granular temperature``, and ``write grid``.

* The ``write grid`` parameter enables one to output the grid as VTU files. The output files can be used for post-processing. By default, it is ``true``. The files with the grid will have ``-grid`` in their names.

* The ``calculate particles average velocity`` chooses whether the average velocity of the particles will be output or not. It generates a file with the average velocity of the particles.

* The ``calculate granular temperature`` chooses whether the average velocity of the particles will be output or not. It generates a file with the average granular temperature of the particles.

.. warning::
 ``calculate particles average velocity``, ``calculate granular temperature``, and ``write grid`` will not work without setting ``Lagrangian post processing`` to ``true``.

* The ``initial step`` is used to choose when to start the calculating either the particle's average velocity, the granular temperature, and/or the grid outputting.

* The ``end step`` is used to choose when to stop the calculating either the particle's average velocity, the granular temperature, and/or the grid outputting.

.. warning::
 Since the default value for the ``end step`` is 0, the ``end step`` must be set according to the duration of the simulation or at least is higher than the ``initial step``, otherwise the particles' average velocity, granular temperature, and/or grid will not be output.

* The ``output frequency`` is used to choose the frequency of creation of output files with the particles' average velocity, granular temperature, and/or grid outputting.

* The ``particles velocity output name`` is used to set a name for the output file with the particle's average velocity. By default, the file is named ``particles_velocity``.

* The ``granular temperature output name`` is used to set a name for the output file with the particle's granular temperature. By default, the file is named ``granular_temperature``.

