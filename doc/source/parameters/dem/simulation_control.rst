==================
Simulation Control
==================

The Simulation Control subsection of DEM simulations is identical to the `CFD <https://lethe-cfd.github.io/lethe/documentation/parameters/cfd/simulation_control.html>`_ in Lethe.

.. note::
    A constant ``time step`` must be defined in the DEM solver. Lethe-DEM does not support adaptative time step scaling method.

.. note::
    ``time step`` in DEM simulations is generally in the range of 1e-7 to 1e-5 seconds. With this small ``time step``, the DEM simulation can capture a single collision in several iterations, which leads to the high accuracy of the simulation. 

.. note::
    Lethe-DEM compares the selected ``time step`` to Rayleigh time-step and prints a warning if the selected ``time step`` is too small or too large. Users should modify the ``time step`` if they see this warning.

The Rayleigh time-step is defined as:

.. math::
    {\Delta}t_{Ra}=\frac{\pi}{2}{d_p}\sqrt{\frac{\rho_p}{G}}(\frac{1}{0.1631\nu+0.8766})

where :math:`{d_p}`, :math:`{\rho_p}`, :math:`{G}`, :math:`{\nu}` denote particle diameter, particle density, shear modulus, and Poisson's ratio.

.. note::
    If the ``output path`` is specified in the parameter handler file, a folder with the specified name should be created before starting simulation in the working directory.

.. code-block:: text

    subsection simulation control
      # DEM time-step
      set time step        = 1e-5

      # Simulation end time
      set time end         = 0.2

      # File output prefix
      set output path      = ./

      # Log frequency
      set log frequency    = 1000

      # Output frequency
      set output frequency = 10000
    end


