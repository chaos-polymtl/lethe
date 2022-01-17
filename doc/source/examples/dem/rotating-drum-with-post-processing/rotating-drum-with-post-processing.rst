==================================
Rotating Drum with Post-Processing
==================================

This is the fourth example of Lethe-DEM. This is a mini-example that only adds Lagrangian post-processing features to the rotating drum example (example 3). Hence, we only explain the post-processing subsection in this example.

Features
----------------------------------
- Solvers: ``dem_3d``
- Rotational boundary
- Load-balancing
- Lagrangian post-processing


Location of the examples
------------------------
 ``/examples/dem/3d_rotating_drum_with_post_processing/rotating_drum_with_post_processing.prm``


Description of the case
-----------------------

This example is identical the rotating drum example. The only difference is that in this example, we use Lagrangian post-processing to obtain granular temperature and average velocity (averaged in cells) distribution.


Parameter file
--------------

Post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we set the variable ``Lagrangian post processing`` equal to true. This enables Lagrangian post-processing calculations. Then we specify the post-processing features that should be obtained. At the moment, Lethe-DEM supports the calculation of ``particle average velocity``, and ``granular temperature``. We set the ``initial step`` and ``end step`` of the post-processing calculations. In the period between initial and end steps, Lethe-DEM calculates and writes the granular temperature and average velocity of particles in cells at a frequency of ``output frequency``. ``particles velocity output name`` and ``granular temperature output name`` define the names of the written post-processing files for average velocity and granular temperature, respectively.

.. code-block:: text

    subsection post-processing
        set Lagrangian post processing				= true
        set calculate particles average velocity	= true
        set calculate granular temperature			= true
        set initial step            				= 8500000
        set end step       							= 9500000
        set output frequency						= 1000
        set particles velocity output name   		= average_velocity
        set granular temperature output name		= granular_temperature
    end

