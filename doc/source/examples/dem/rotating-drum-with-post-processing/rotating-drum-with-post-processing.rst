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

Mesh
~~~~~

.. code-block:: text

    subsection mesh
        set type                 				= dealii
        set grid type      	     				= cylinder
        set grid arguments       				= 0.12:0.18
        set initial refinement   				= 4
    end


Insertion info
~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection insertion info
        set insertion method								= non_uniform
        set inserted number of particles at each time step  = 37555
        set insertion frequency            		 			= 150000
        set insertion box minimum x            	 			= -0.175
        set insertion box minimum y            	        	= -0.07
        set insertion box minimum z            	        	= 0
        set insertion box maximum x            	        	= 0.175
        set insertion box maximum y           	 			= 0.07
        set insertion box maximum z            	        	= 0.09
        set insertion distance threshold					= 1.575
        set insertion random number range					= 0.05
        set insertion random number seed					= 19
    end


Lagrangian physical properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection lagrangian physical properties
        set gx            		 						= 0.0
        set gy            		 						= 0.0
        set gz            		 						= -9.81
        set number of particle types	                = 1
            subsection particle type 0
            set size distribution type					= uniform
                set diameter            	 			= 0.003
                set number              				= 226080
                set density particles  	 				= 2500
                set young modulus particles         	= 100000000
                set poisson ratio particles          	= 0.24
                set restitution coefficient particles	= 0.97
                set friction coefficient particles      = 0.3
                set rolling friction particles         	= 0.01
        end
        set young modulus wall            				= 100000000
        set poisson ratio wall            				= 0.24
        set restitution coefficient wall           		= 0.85
        set friction coefficient wall         			= 0.35
        set rolling friction wall         	      	  	= 0.01
    end


Model parameters
~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection model parameters
      set contact detection method 		   		 	= dynamic
      set dynamic contact search size coefficient	= 0.8
      set neighborhood threshold				 	= 1.3
      set load balance method				 		= once
  	  set load balance step					 		= 150000
      set particle particle contact force method	= hertz_mindlin_limit_overlap
      set particle wall contact force method        = nonlinear
      set integration method				 		= velocity_verlet
    end


Boundary Condition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection DEM boundary conditions
      set number of boundary conditions         = 1
        subsection boundary condition 0
            set boundary id						= 4
            set type              				= rotational
            set rotational speed				= 11.6
            set rotational vector x				= 1
            set rotational vector y				= 0
            set rotational vector z				= 0
        end
    end


Simulation control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection simulation control
      set time step                 		 = 1e-6
      set time end       					 = 15
      set log frequency				         = 1000
      set output frequency            		 = 1000
    end


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


Running the simulation
----------------------
This simulation can be launched (in parallel mode on 64 processes) by:

.. code-block:: text

  mpirun -np 64 dem_3d rotating_drum_with_post_processing.prm

