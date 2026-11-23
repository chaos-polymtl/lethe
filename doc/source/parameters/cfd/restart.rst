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
    set checkpoint          = false
    set checkpoint control  = iteration
    set frequency           = 1
    set time frequency      = 1
    set wall time frequency = 7200

    # Output/input files
    set filename            = restart

    # Restarting parameter
    set restart             = false
  end

* ``checkpoint``: controls if a checkpoint of the simulation is created. All the files needed to restart the simulation from this checkpoint (such as ``.pvdhandler``, ``.triangulation``) will be written.

* ``checkpoint control``: what sets when the checkpoints are written, like ``output control`` for the results (see :doc:`simulation_control`). With ``iteration``, a checkpoint is written every ``frequency`` iterations; with ``time``, every ``time frequency`` of simulated time; with ``wall_time``, every ``wall time frequency`` of elapsed wall time. Whatever the control, a checkpoint is also written at the last time step.

* ``frequency``: number of iteration before the first checkpoint, and between each subsequent checkpoint, for the ``iteration`` control.

* ``time frequency``: simulated time (s) between checkpoints, for the ``time`` control. A checkpoint is written at the first time step that reaches each multiple of this value. This control cannot be used with a ``steady`` simulation.

* ``wall time frequency``: elapsed wall time (s) between checkpoints, for the ``wall_time`` control. It is counted from the start of the simulation, then from the last checkpoint.

.. note::
  In DEM simulations, which usually have a large ``log frequency``, the elapsed wall time of the ``wall_time`` control is only checked at log iterations, to avoid a communication between all the processes at every iteration.

.. tip::
  Checkpointing is very intensive from an I/O point of view, consequently it should not be done at every iteration.

.. note::
  Checkpoints alternate between two slots, so a simulation interrupted while writing a checkpoint can restart from the previous one. Their files end with ``_0`` or ``_1`` (e.g., ``<filename>_0.triangulation``, ``force_00_0.checkpoint``, ...), and ``<filename>.checkpoint_controller`` gives the slot of the last checkpoint, from which a simulation restarts.

* ``filename``: prefix for the files. This parameter is used to output checkpoint files as well as to load files during a restart.

.. warning::
  Files are saved (when checkpointing) and searched for (when restarting) in the ``output path`` specified in subsection :doc:`simulation_control`.

* ``restart``: controls if the simulation is to be restarted from a previous checkpoint. 
   * ``set restart = false``: for the first checkpoint creation. 
   * ``set restart = true``: the simulation will start from the last checkpoint and will finish at the provided ``time end`` (see :doc:`simulation_control`).



