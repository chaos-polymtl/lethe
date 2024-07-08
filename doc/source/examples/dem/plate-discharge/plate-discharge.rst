==================================
Plate Discharge
==================================

This example compares the angles of repose and performance results of a plate discharging particles with performance enhancement methods.

----------------------------------
Features
----------------------------------

- Solvers: ``lethe-particles``
- Three-dimensional problem
- Uses the `adaptive sparse contacts <../../../parameters/dem/model_parameters.html#adaptive-sparse-contacts-asc>`_
- Uses the `dynamic load balancing <../../../parameters/dem/model_parameters.html#load-balancing>`_
- Post-processes results and compares them to the literature


---------------------------
Files Used in this Example
---------------------------

All the files mentioned below are located in the example folder ``examples/dem/plate-discharge``.

- There are 4 parameters files: a baseline case and three other cases with different features using one or a combination of enhanced performance methods. The parameters files are:

  .. list-table::
     :width: 100%
     :widths: 30 30 30
     :header-rows: 1
     :align: center

     * - Name of the .prm file
       - `Adaptive Sparse Contacts <../../../parameters/dem/model_parameters.html#adaptive-sparse-contacts-asc>`_
       - `Load Balancing <../../../parameters/dem/model_parameters.html#load-balancing>`_
     * - ``plate-discharge_base.prm``
       -
       -
     * - ``plate-discharge_asc.prm``
       - ×
       -
     * - ``plate-discharge_lb.prm``
       -
       - ×
     * - ``plate-discharge_asc-lb.prm``
       - ×
       - ×

- Those parameters files are ready for the simulations. Since we want the physical data from the cases and the computational performance information of the simulation not affected by the writing of data files, we run 2 sets of simulation (performance and data). The performance analysis parameter files are in the folder ``performance/``, and the ones for the data analysis are in the folder ``data/``.

-----------------------
Description of the Case
-----------------------

This example simulates the discharge of particles at the side of a plate in a rectangular container in order to get the angle of repose of the granular material as done by Zhou *et al*. [#zhou2002]_. The example compares the angles of repose and the performance of the simulations with the use of the adaptive sparse contacts and the load balancing methods. The angles are also compared to the literature.

.. figure:: images/plate-discharge-diagram.png
    :width: 50%
    :alt: Plate discharge diagram
    :align: center

    Diagram of the container (black) with the plate (green) and the particles (blue).


-------------------
DEM Parameter files
-------------------

Baseline case simulation
~~~~~~~~~~~~~~~~~~~~~~~~

In this section we introduce the different sections of the parameter file ``plate-discharge_base.prm`` with do not use any performance enhancement methods.


Simulation Control
------------------

The simulation lasts 15 seconds and the DEM time step is 0.0001 seconds. The output are generated every 0.01 seconds for the simulation for data analysis.

.. code-block:: text

   subsection simulation control
     set time step        = 1e-4
     set time end         = 15
     set log frequency    = 500
     set output frequency = 100
     set output path      = ./output_base/
   end

Mesh
----

The rectangular container is a 1 m x 1 m x 0.2 m box with a 0.9 m x 0.2 m plate placed at a height of 0.4 m.

.. code-block:: text

   subsection mesh
     set type               = dealii
     set grid type          = subdivided_hyper_rectangle
     set grid arguments     = 5,5,1 : -0.5, 0.0, 0.0 : 0.5, 1.0, 0.2 : true
     set initial refinement = 3
   end

Lagrangian Physical Properties
------------------------------

The lagrangian properties are relatively arbitrary. The simulation contains 52000 particles with a diameter of 0.01 m, a density of 2400 kg/m³. Both properties of particle-particles and particle-wall interactions are the same.

.. code-block:: text

   subsection lagrangian physical properties
     set g                        = 0, -9.81, 0.0
     set number of particle types = 1
     subsection particle type 0
       set size distribution type            = uniform
       set diameter                          = 0.01
       set number of particles               = 52000
       set density particles                 = 2400
       set young modulus particles           = 1e6
       set poisson ratio particles           = 0.3
       set restitution coefficient particles = 0.9
       set friction coefficient particles    = 0.2
       set rolling friction particles        = 0.1
     end
     set young modulus wall           = 1e6
     set poisson ratio wall           = 0.3
     set restitution coefficient wall = 0.9
     set friction coefficient wall    = 0.2
     set rolling friction wall        = 0.1
   end

Insertion Info
--------------

The particles are inserted above the plate with the volume insertion method.

.. code-block:: text

   subsection insertion info
     set insertion method                               = volume
     set inserted number of particles at each time step = 52000
     set insertion frequency                            = 20000
     set insertion box points coordinates               = -0.45, 0.4, 0 : 0.45, 1.0, 0.2
     set insertion distance threshold                   = 1.25
     set insertion maximum offset                       = 0.1
     set insertion prn seed                             = 20
     set insertion direction sequence                   = 0, 2, 1
   end

Floating Walls
--------------

The particles all stay on the plate with floating walls that are placed vertically at both extremities of the plate. The walls are removed after 0.75 seconds of simulation, starting the discharge.

.. code-block:: text

   subsection floating walls
     set number of floating walls = 2
     subsection wall 0
       subsection point on wall
         set x = -0.45
         set y = 0
         set z = 0
       end
       subsection normal vector
         set nx = 1
         set ny = 0
         set nz = 0
       end
       set start time = 0
       set end time   = 0.75
     end
     subsection wall 1
       subsection point on wall
         set x = 0.45
         set y = 0
         set z = 0
       end
       subsection normal vector
         set nx = 1
         set ny = 0
         set nz = 0
       end
       set start time = 0
       set end time   = 0.75
     end
   end


Solid Objects
-------------

The plate is a solid object with a simple mesh of 2 triangles placed at 0.4 m in the container.

.. code-block:: text

   subsection solid objects
     subsection solid surfaces
       set number of solids = 1
       subsection solid object 0
         subsection mesh
           set type                = gmsh
           set file name           = plate.msh
           set simplex             = true
           set initial translation = 0, 0.4, 0
         end
       end
     end
   end

Model Parameters
----------------

The model parameters are quite standard for a DEM simulation with the non-linear Hertz-Mindlin contact force model, a constant rolling resistance torque, and the velocity Verlet integration method. For the baseline case, we do not use any performance enhancement methods.

.. code-block:: text

   subsection model parameters
     subsection contact detection
       set contact detection method                = dynamic
       set dynamic contact search size coefficient = 0.9
       set neighborhood threshold                  = 1.3
     end
     subsection load balancing
       set load balance method = none
     end
     set particle particle contact force method = hertz_mindlin_limit_overlap
     set rolling resistance torque method       = constant_resistance
     set particle wall contact force method     = nonlinear
     set integration method                     = velocity_verlet
     subsection adaptive sparse contacts
       set enable adaptive sparse contacts = false
     end
   end


Timer
-------

The timer is enabled since we want to profile the computational performance of the simulations. The timer prints the total wallclock time elapsed since the start at every log frequency iterations.

.. code-block:: text

   subsection timer
     set type = iteration
   end


ASC Simulation
~~~~~~~~~~~~~~~~~~

The only differences between ``plate-discharge_base.prm`` and ``plate-discharge_asc.prm`` are the enabling of the adaptive sparse contacts and the name of the folder for outputs.

Model Parameters
----------------

Here the ASC is enabled with a granular temperature threshold of 0.0001 m²/s² and a solid fraction threshold of 0.4.

.. code-block:: text

   subsection model parameters
     subsection contact detection
       set contact detection method                = dynamic
       set dynamic contact search size coefficient = 0.9
       set neighborhood threshold                  = 1.3
     end
     subsection load balancing
       set load balance method = none
     end
     set particle particle contact force method = hertz_mindlin_limit_overlap
     set rolling resistance torque method       = constant_resistance
     set particle wall contact force method     = nonlinear
     set integration method                     = velocity_verlet
     subsection adaptive sparse contacts
       set enable adaptive sparse contacts = true
       set granular temperature threshold  = 1e-4
       set solid fraction threshold        = 0.4
     end
   end


Load Balancing Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~

The only differences between ``plate-discharge_base.prm`` and ``plate-discharge_lb.prm`` are the usage of the load balancing and the name of the folder for outputs.

Model Parameters
----------------

Here, the dynamic load balancing checks if a load balancing is needed every 2500 iterations with a load threshold of 0.5.

.. code-block:: text

   subsection model parameters
     subsection contact detection
       set contact detection method                = dynamic
       set dynamic contact search size coefficient = 0.9
       set neighborhood threshold                  = 1.3
     end
     subsection load balancing
       set load balance method     = dynamic
       set threshold               = 0.5
       set dynamic check frequency = 2500
     end
     set particle particle contact force method = hertz_mindlin_limit_overlap
     set rolling resistance torque method       = constant_resistance
     set particle wall contact force method     = nonlinear
     set integration method                     = velocity_verlet
     subsection adaptive sparse contacts
       set enable adaptive sparse contacts = false
     end
   end


ASC with Load Balancing Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The only differences between ``plate-discharge_base.prm`` and ``plate-discharge_asc-lb.prm`` are the usage of the ASC method with the load balancing, and the name of the folder for outputs.

Model Parameters
----------------

Here, we use the ASC with the dynamic load balancing still checking if a load balancing is needed every 2500 iterations with a load threshold of 0.5. In with case, the mobility status of the cells from the adaptive sparse contacts will influence the weight, or the computational contribution, of the cell in the load balancing evaluation. The weight factor of the active cells is 0.7, and the weight factor of the inactive cells is 0.5, while the mobile cells always have a fixed weight factor of 1.

.. code-block:: text

   subsection model parameters
     subsection contact detection
       set contact detection method                = dynamic
       set dynamic contact search size coefficient = 0.9
       set neighborhood threshold                  = 1.3
     end
     subsection load balancing
       set load balance method     = dynamic_with_sparse_contacts
       set threshold               = 0.5
       set dynamic check frequency = 2500
       set active weight factor    = 0.7
       set inactive weight factor  = 0.5
     end
     set particle particle contact force method = hertz_mindlin_limit_overlap
     set rolling resistance torque method       = constant_resistance
     set particle wall contact force method     = nonlinear
     set integration method                     = velocity_verlet
     subsection adaptive sparse contacts
       set enable adaptive sparse contacts = true
       set granular temperature threshold  = 1e-4
       set solid fraction threshold        = 0.4
     end
   end


-----------------------
Running the Simulations
-----------------------

Simulations can be launched individually with the executable ``lethe-particles`` and the parameter files with the logging of the display in the terminal.
To make things easier a script is provided to run all the simulations in sequence from the ``dem/3d-plate-discharge/`` folder.

In order to run the simulations for the performance analysis, you can use the following command:

.. code-block:: text
  :class: copy-button

  bash run-performance-simulation.sh

Which corresponds to:

.. code-block:: bash

  simulations=("base" "asc" "lb" "asc-lb")

  cd performance/

  for sim in "${simulations[@]}"
  do
     echo "Running the $sim simulation"
     time mpirun -np 8 lethe-particles plate-discharge_$sim.prm | tee log_$sim.out
  done

Or you can run the simulations in the ``performance/`` folder with the following commands:

.. code-block:: text
  :class: copy-button

  time mpirun -np 8 lethe-particles plate-discharge_base.prm | tee log_base.out
  time mpirun -np 8 lethe-particles plate-discharge_asc.prm | tee log_asc.out
  time mpirun -np 8 lethe-particles plate-discharge_lb.prm | tee log_lb.out
  time mpirun -np 8 lethe-particles plate-discharge_asc-lb.prm | tee log_asc-lb.out

In order to run the simulations for the data analysis, you can use the following script:

.. code-block:: text
  :class: copy-button

  bash data-performance-simulation.sh

.. note::
   Running the simulations for the performance analysis using 8 cores takes between 25 and 45 minutes per simulation, for a total of around 2 hours. Running the simulations for data analysis takes a few minutes longer per simulation.

-------
Results
-------

The simulations should look like the following video:

[WIP: Add video here]


Post-Processing Code
~~~~~~~~~~~~~~~~~~~~

The data is extracted with the Lethe PyVista tool and post-processed with custom functions in the files ``pyvista_utilities.py`` and ``log_utilities.py``.
Extraction, post-processing and plotting are automated in the script ``plate-discharge_post-processing.py``:

.. code-block:: text
  :class: copy-button

  python3 plate-discharge_post-processing.py

The script will generate the figures. If you want to modify the path or the filenames, you have to modify the script.

Performance Analysis
~~~~~~~~~~~~~~~~~~~~

The log files (outputs displayed in the terminal) are read to extract the simulation and wall times.

The speedup is calculated with the baseline case as the reference. The results are plotted in the following figure, where the solid lines show the walltime during the simulation, the dashed lines show the speedup, and the points show to total speedup.


.. figure:: images/performance.png
   :alt: Performance results
   :align: center
   :name: plate-discharge-performance-graph

   The walltime during the simulations (solid line) and the speedup (dashed line) for the performance enhancement methods with the Adaptive Sparse Contacts (ASC) and the Load Balancing (LB) compared to the baseline case.

.. note::
   The slight oscillations of the speedup are caused by the scientific notation format of the walltime by the timer feature after 1000 seconds. The walltimes are attenuated by the moving average, but the division operation for the speedup accentuates the lack of time precision.

Without going into much details, the load balancing method will help the performance of the simulation from the start, since the particles move within the domain during the discharge. The load balancing allows to distribute the particles, therefore all their related computation, more evenly between the cores. Once the discharge of the particles is mostly done, and only a few particles are still falling from the top part, the performance gain brought by the load balancing stays constant since the load across the cores is already balanced.
The adaptive sparse contacts method will help the performance of the simulation mostly when there are large areas of motionless particles. As it was showed in the video, those areas are located in the core of the pile at the top and at the corners of settled particles below the plate. This explains why the ASC gives a limited performance gain at the start of the simulation (only from the core of the pile) and an increasing gain through the simulation (accumulation of motionless particles at the bottom part). Given that both methods help the computation performance at different times, the combination of both methods gives the best performance as observed.


Angle of Repose
~~~~~~~~~~~~~~~

The angle of repose is calculated from the data extracted from the VTU output files. Angles of repose are calculated from the pile of particles on the plate for comparison with the literature, and from the piles formed by the discharge for curiosity.

The configuration of the case gives a symmetrical formation of the piles, meaning that there are 2 angles of repose to calculate over and below the plate.
In order to show how the results may fluctuated in regards of that, we show the angle obtained from the particle positions from the left and the right sides.


The angles of repose are calculated by linear regressions from the highest particle positions in y-axis from -0.35 m to -0.15 m for the left angles and from 0.15 m to 0.35 m for the right angles in x-axis. The following figure shows the areas where the angles are calculated. The given angles of reposes are the linear regressions from the positions with absolute x coordinates.

.. figure:: images/angle-areas.png
   :alt: Angle of repose areas
   :align: center
   :name: plate-discharge-angle-areas

   The areas where the angle of repose is calculated for the left (blue) and right (red) sides of the piles.

.. note::
   The calculated angles ignore the wall effects since all the particles in the depth of the container are taken into account.

According to Zhou *et al.* [#zhou2002]_, the angle of repose from this type of configuration is calculated with the following formula:

.. math::
   \phi = 68.61 \mu_{\text{f,pp}}^{0.27} \mu_{\text{f,pw}}^{0.22} \mu_{\text{r,pp}}^{0.06} \mu_{\text{r,pw}}^{0.12} d_p^{-0.2}


where :math:`\mu_{\text{f,pp}}` and :math:`\mu_{\text{f,pw}}` are the friction coefficients of the particle-particle and particle-wall interactions, respectively, :math:`\mu_{\text{r,pp}}` and :math:`\mu_{\text{r,pw}}` are the rolling friction coefficients, and :math:`d_p` is the particle diameter.


The meaning of the rolling friction coefficient by the authors [#zhou2002]_ is different than found in Lethe. They express the coefficient as aa length in the `rolling friction model <../../../theory/multiphase/cfd_dem/dem.html#rolling-friction-models>`_. However, they also use the constant torque, therefore the rolling friction coefficient in Lethe as to be multiplied by the effective radius of the particle for results comparison:

.. math::
   \mu_{\text{r}}^{\text{eqt}} = \mu_{\text{r}}^{\text{lethe}}D_p

.. figure:: images/angle-of-repose.png
   :alt: Angle of repose results
   :align: center
   :name: plate-discharge-angle-graph

   The angles of repose calculated from the simulation data. The solid lines are the angles computed from the highest particles on both side, while the shaded areas represent the angles for the left and the right.

The theoretical angle of repose is 19.7°. We did not compute the mean of the angled of repose in order to compare the results with the literature since, even after 15 seconds of simulation, some particles are still falling from the top and the angles are still not oscillating around the same value. We can however state that the angles are close to the literature.

Here we can see that the top angles from all simulations are in a range of around ±1.5 from the baseline case, which can be considered as a good agreement. We can clearly see a trend in the bottom angles using the ASC. The angles of repose are about 2° below the baseline and load balancing cases. It seems to be caused by the accumulation of particles at the bottom of the piles.


----------
References
----------

.. [#zhou2002] \Y.C. Zhou, B.H. Xu, A.B. Yu, P. Zulli, “An experimental and numerical study of the angle of repose of coarse spheres,” *Powder Technology*, vol. 125, pp. 45-54, 2002. doi: `10.1016/S0032-5910(01)00520-4 <https://doi.org/10.1016/S0032-5910(01)00520-4>`_\.
