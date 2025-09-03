=======
Restart
=======

This subsection controls the checkpointing and restarting of time-dependent simulations.

This feature is used both to write data at checkpoints (`checkpointing`) and restart a simulation (`restarting`) using data from the last checkpoint.

.. tip::
  It is highly advised to use it if you fear you might either run out of memory/simulation time during the simulation.

The default parameters are:

.. code-block:: text

  subsection restart
    # Checkpointing parameters
    set checkpoint = false
    set frequency  = 1

    # Output/input files
    set filename   = restart

    # Restarting parameter
    set restart    = false
  end

* ``checkpoint``: controls if a checkpoint of the simulation is created. All the files needed to restart the simulation from this checkpoint (such as ``.pvdhandler``, ``.triangulation``) will be written.

* ``frequency``: number of iteration before the first checkpoint, and between each subsequent checkpoint. 

.. tip::
  Checkpointing is very intensive from an I/O point of view, consequently it should not be done at every iteration. The checkpointing frequency is in accordance to the number of ``time step`` (see :doc:`simulation_control`).

.. warning::
  Each checkpoint erases the previous one, so a simulation can only be restarted from the last check-point.

* ``filename``: prefix for the files. This parameter is used to output checkpoint files as well as to load files during a restart.

.. warning::
  Files are saved (when checkpointing) and searched for (when restarting) in the ``output path`` specified in subsection :doc:`simulation_control`.

* ``restart``: controls if the simulation is to be restarted from a previous checkpoint. 
   * ``set restart = false``: for the first checkpoint creation. 
   * ``set restart = true``: the simulation will start from the last checkpoint and will finish at the provided ``time end`` (see :doc:`simulation_control`).



