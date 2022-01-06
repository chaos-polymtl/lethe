Simulation Control
-------------------
This subsection contains the general information of the simulation, including ``time step``, ``end time``, ``output path``, ``log frequency``, and ``output frequency``. It is the most commonly modified section for a simulation. ``time step`` in DEM simulations is generally in the range of 1e-7 to 1e-5 seconds. With this small ``time step``, the DEM simulation can capture a single collision in several iterations, which leads to the high accuracy of the simulation. 

.. note::
    Lethe-DEM compares the selected ``time step`` to Rayleigh time-step and prints a warning if the selected ``time step`` is too small or too large. Users should modify the ``time step`` if they see this warning.

.. note::
    It should be mentioned that if the ``output path`` is specified in the parameter handler file, a folder with the specified name should be created before starting simulation in the working directory.

.. code-block:: text

 subsection simulation control
  # DEM time-step 
  set time step                         = 1e-5

  # Simulation end time
  set time end                          = 0.2

  # File output prefix
  set output path                       = ./

  # Log frequency
  set log frequency                     = 1000

  # Output frequency
  set output frequency                  = 10000
 end


