Finite element interpolation (FEM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This subsection controls the restart and checkpointing of the simulation. 
This feature is used to restart a time-dependent simulation using the data written on the last checkpoint before the simulation stopped or crashed.
The default parameters are:

.. code-block:: text

  subsection restart
    set checkpoint = false
    set restart    = false
    set filename   = restart
    set frequency  = 1
  end

* The ``checkpoint`` option enables creating a checkpoint of the simulation when it is set to ``true``. 

.. tip::
  You may want to enable it when launching a time-dependent simulation if you fear you might either run out of memory/simulation time during the simulation.

* The ``restart`` option controls if the simulation will be restarted from a previous checkpoint. Turn this parameter to ``false`` for the first checkpoint creation. The simulation will start from the last checkpoint when the parameter is set to ``true`` and will finish at the provided time end (see `Simulation Control <https://lethe-cfd.github.io/lethe/parameters/cfd/simulation_control.html>`_).

* The ``filename`` parameter fixes the prefix of all the checkpointing files. 

.. warning::

  Restart files are searched for in the output path specified in subsection simulation control.

* The ``frequency`` parameter controls the checkpointing frequency. Checkpointing is very intensive from an I/O point of view, consequently it should not be done at every iteration. The checkpointing frequency is is accordance to the number of ``time step`` (see `Simulation Control <https://lethe-cfd.github.io/lethe/parameters/cfd/simulation_control.html>`_).
