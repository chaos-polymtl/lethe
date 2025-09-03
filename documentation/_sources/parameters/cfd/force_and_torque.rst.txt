============================
Force and Torque Calculation
============================

This subsection controls the post-processing of the force and the torque on the boundary conditions. Using these options, the user can calculate the force  acting on a wall or the torque acting on a boundary. 

.. code-block:: text

  subsection forces
    # Enable calculation of force
    set calculate force       = true

    # Enable calculation of torque
    set calculate torque      = true

    # Calculation frequency
    set calculation frequency = 1

    # File output force prefix
    set force name            = force

    # File output torque prefix
    set torque name           = torque

    # Output frequency
    set output frequency      = 1

    # Calculation frequency
    set output precision      = 10

    # State whether information from the non-linear solver should be printed (quiet or verbose).
    set verbosity             = quiet
  end

* ``calculate force`` enables the calculation of the force on all boundaries. If multiple walls bear the same ID, the total force for this ID will be calculated.

* ``calculate torque`` enables the calculation of the torque on all boundaries that bear an ID. If multiple walls bear the same ID, the total torque for this ID will be calculated.

* ``calculation frequency`` is an integer that specifies the frequency of the calculation of the force. Setting ``calculation frequency=10`` means that forces and torques will be calculated every 10 iterations. Calculating the forces and the torques on the boundaries is a very cheap operation and, consequently, there is not much to optimize by using a larger frequency.

* ``output frequency`` specifies at which frequency the force and torque files are written to the disk. These files are small ``.dat`` files, but writing them every iteration can significantly slow down the simulation depending on the file system. On clusters, such as the one provided by Compute Canada, it is preferable to write these files every 10-100 iteration to prevent slowing down the file system.

* ``output precision`` specifies the number of significant digits in the output files. There is little to gain in reducing the default value here.

* ``force name`` and ``torque name`` specify the prefix of the files used to output the force and the torque respectively. The files are also numbered by the boundary ID and bear the ``.dat`` suffix. If ``force name=force``, the name of the file for boundary ID 1 will be ``force.01.dat``.

* ``verbosity`` specifies if the force and torque values are printed to the terminal at the end of every iteration. >This does not affect the printing of output files.



