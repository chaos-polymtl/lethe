===============================
3D Dam Break With an Obstacle
===============================

This example simulates a dam break experiment from the Maritime Research Institute Netherlands (MARIN) `[1] <https://www.spheric-sph.org/tests/test-02>`_.

.. warning::
    This example displays the need for improvement of low-viscosity fluid flow simulation of the current numeric model. Further work will be done to improve this aspect of the model.

----------------------------------
Features
----------------------------------
- Solver: ``gls_navier_stokes_3d`` (Q1-Q1)
- Two phase flow handled by the Volume-of-Fluids (VOF) approach with interface sharpening
- Mesh adaptation using phase fraction
- Unsteady problem handled by an adaptive BDF1 time-stepping scheme
- The use of a python script for post-processing data


Files used in this example
~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Parameter file: ``examples/multiphysics/3d-dam-break/3d-dam-break.prm``
- Geometry file: ``examples/multiphysics/3d-dam-break/tank_with_obstacle.geo``
- Python script for post-processing data: ``examples/multiphysics/3d-dam-break/3d-dam-break_postprocess.py``
- Experimental data file: ``examples/multiphysics/3d-dam-break/experimental_data.txt``


.. _Description of the case:

-------------------------
Description of the case
-------------------------

For this example, the simulated fluids are water and air. Initially, the water is at rest on the right side of the tank (represented in dashed blue lines in the figures below). At :math:`t = 0` s, the gate opens instantaneously and the water starts flowing under the action of gravity, :math:`\mathbf{g} = (-9.81 \  \mathbf{j}) \frac{\text{m}}{\text{s}^2}`. The tank in which this experiment happens has the following dimensions: :math:`3.22 \times 1.00 \times 1.00` m. On all boundaries, ``slip`` conditions were applied. On the left side of the tank, a rectangular box-shaped obstacle is presented (colored in grey in the figures).


Along the x-axis, the water height is measured at 4 different positions. These positions are represented by red crosses in the figure below.

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/3d-dam-break-figure.png                                                                       |
|     :alt: View from above the tank. Initially, the water is at rest on the right side of the tank. The            |
|      water column is 1.228 m long. 0.6635 m from the left edge, there is a rectangular box-shaped obstacle.       |
|      This obstacle is 1.1675 m away from the resting water. The water height is measured along the x-axis         |
|      at respectively 0.496 m, 0.992 m, 1.488 m, and 2.632 m from the left side of the tank.                       |
|     :align: center                                                                                                |
|     :name: Initial state top view                                                                                 |
|                                                                                                                   |
|     Initial state top view                                                                                        |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

|

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/3d-dam-break-figure_side.png                                                                  |
|     :alt: View from the side of the tank. Initially, the water is at rest on the right side of the tank.          |
|      The water column is 1.228 m long, and 0.55 m in height. 0.6635 m from the left edge, there is a rectangular  |
|      box-shaped obstacle. This obstacle is 1.1675 m away from the resting water. The water height is measured     |
|      along the x-axis at respectively 0.496 m, 0.992 m, 1.488 m, and 2.632 m from the left side of the tank.      |
|     :align: center                                                                                                |
|     :name: Initial state side view                                                                                |
|                                                                                                                   |
|     Initial state side view                                                                                       |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

|

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/geo.png                                                                                       |
|     :alt: 3D view of the initial system                                                                           |
|     :align: center                                                                                                |
|     :name: Initial state in 3D                                                                                    |
|                                                                                                                   |
|     Initial state in 3D                                                                                           |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

-----------------
Parameter file
-----------------

Simulation control
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Time integration is handled by a 1st order backward differentiation scheme (`bdf1`), for a :math:`6` s simulation time with an initial time step of :math:`0.001` s. Time-step adaptation is enabled using ``adapt=true``
and the max CFL is :math:`0.5`.

.. note::
    This example uses an adaptive time-stepping method, where the time-steps are modified during the simulation to keep the maximum value of the CFL condition below the given threshold (0.5).

.. code-block:: text

    #---------------------------------------------------
    # Simulation Control
    #---------------------------------------------------

    subsection simulation control
      set method                       = bdf1
      set time end                     = 6
      set time step                    = 0.001
      set adapt                        = true
      set max cfl                      = 0.5
      set output name                  = 3d-dam-break
      set output frequency             = 5
      set output path                  = ./output/
      set adaptative time step scaling = 1.05
      set output boundaries            = true
    end

Multiphysics
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``multiphysics`` subsection enables to turn on `(true)`
and off `(false)` the physics of interest. Here ``VOF`` is chosen.
Note that the fluid dynamics are solved by default.

.. code-block:: text

    #---------------------------------------------------
    # Multiphysics
    #---------------------------------------------------

    subsection multiphysics
      set VOF = true
    end

Physical properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``physical properties`` subsection defines the physical properties of the fluids. In this example, we need two fluids with densities of :math:`1.204 \ \frac{kg}{m^3}` (air) and :math:`1000 \ \frac{kg}{m^3}` (water). However, the current numerical model was not able to solve with the real dynamic viscosities of the fluids. Therefore, they were altered in order to run the simulation.

.. warning::
    Altering the dynamic viscosities of the fluids will surely have an impact on the results. We will show this impact in the `<Results_>`_ section.

.. code-block:: text

    #---------------------------------------------------
    # Physical Properties
    #---------------------------------------------------

    subsection physical properties
      set number of fluids = 2
      subsection fluid 0
        set density             = 1.204
        set kinematic viscosity = 0.01516
      end
      subsection fluid 1
        set density             = 1000
        set kinematic viscosity = 0.001
      end
    end

Initial condition
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the ``initial condition`` subsection, we need to define the interface between the two fluids. We define this interface by using a function expression in the ``VOF`` subsection of the ``initial condition``.

.. code-block:: text

    #---------------------------------------------------
    # Initial Condition
    #---------------------------------------------------

    subsection initial conditions
      set type = nodal
      subsection uvwp
        set Function expression = 0; 0; 0; 0
      end

      subsection VOF
        set Function expression = if (x>1.992 & z<0.55 & y>=-0.5, 1, 0)
      end
    end

Source term
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the ``source term`` subsection, we define the gravitational acceleration.

.. code-block:: text

    #---------------------------------------------------
    # Source term
    #---------------------------------------------------

    subsection source term
      set enable = true
      subsection xyz
        set Function expression = 0;0;-9.81;0
      end
    end

VOF
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the ``VOF`` subsection, we enable ``interface sharpening`` to reconstruct the interface and keep it sharp during the simulation. Here we use the default ``constant`` method for interface sharpening.

.. code-block:: text

    #---------------------------------------------------
    # VOF
    #---------------------------------------------------

    subsection VOF
      subsection interface sharpening
        set enable                  = true
        set threshold               = 0.5
        set interface sharpness     = 1.5
        set frequency               = 10
        set type                    = constant
      end
    end

Mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the ``mesh`` subsection, we specify the mesh used in this example. The structured mesh used in this example can be generated from the ``tank.geo`` file using `Gmsh <https://gmsh.info/#Download>`_. The initial refinement is set to :math:`3`.

.. code-block:: text

    #---------------------------------------------------
    # Mesh
    #---------------------------------------------------
    subsection mesh
        set type                 = gmsh
        set file name            = tank.msh
        set initial refinement   = 3
    end


Mesh Adaptation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``mesh adaptation`` section controls the dynamic mesh adaptation. Here, we choose ``phase`` and ``pressure`` as the ``refinement variables``. The maximum and minimum refinement levels are respectively set to :math:`4` and :math:`2`.

.. code-block:: text

    #---------------------------------------------------
    # Mesh Adaptation
    #---------------------------------------------------

    subsection mesh adaptation
      set type                     = kelly
      set variable                 = phase, pressure
      set fraction type            = fraction
      set max refinement level     = 4
      set min refinement level     = 2
      set frequency                = 2
      set fraction refinement      = 0.999, 0.4
      set fraction coarsening      = 0.001, 0.05
      set initial refinement steps = 5
    end


-----------------------
Running the simulation
-----------------------

We call the gls_navier_stokes_3d by invoking:

``mpirun -np $number_of_CPU_cores gls_navier_stokes_3d 3d-dam-break.prm``

.. warning::
    Make sure to compile Lethe in `Release` mode and run in parallel using ``mpirun``. This simulation took :math:`\approx` 15.5 hours on 40 processes (runned on the `Béluga <https://docs.alliancecan.ca/wiki/B%C3%A9luga/en>`_ cluster).

.. _Results:

-----------------
Results
-----------------

The following video shows the results of the simulation:

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/XGp7pxBQvWY" frameborder="0" allowfullscreen></iframe>


In the following figure, we compare the water height evolution at 4 the positions mentioned in the `<Description of the case_>`_ section with the experimental results obtained from MARIN (available `here <https://www.spheric-sph.org/tests/test-02>`_):

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/H1_to_H4_evolution.png                                                                        |
|     :alt: In this figure, the water height evolutions are compared with the experimental results of MARIN.        |
|      In the top left corner, we compare the evolution of the water height at 0.496 m away from the left side of   |
|      the tank. In the top right corner, we compare the evolution of the water height at 0.992 m away from the     |
|      left side of the tank. In the bottom left corner, we compare the evolution of the water height at 1.488 m    |
|      away from the left side of the tank. In the bottom right corner, we compare the evolution of the water       |
|      height at 1.638 m away from the left side of the tank.                                                       |
|     :align: center                                                                                                |
|     :name: Comparison of the water height at different position in the tank with the experimental data of MARIN   |
|                                                                                                                   |
|     Comparison of the water height evolution                                                                      |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

As we can see, the simulated general evolution of the height seems to follow the experimentation results. However, on all 4 subplots, we notice that the height is overestimated. We also notice a slight shift to the right for :math:`H2`,  :math:`H3`, and :math:`H4` evolutions. These observations may be explained by the "highly viscous air" (fluid 0) that acts as an obstacle to the free flow of the water. Additionally, fluid 1 representing the water is 1000 times more viscous than regular water. With these results, we can see that the model needs to be improved to be able to accurately simulate low-viscosity fluids such as air.


-----------
References
-----------


`[1] <https://www.spheric-sph.org/tests/test-02>`_ Issa, R., & Violeau, D. (2006). Test-case 2, 3D dambreaking, Release 1.1. ERCOFTAC, SPH European Research Interest Community SIG, Électricité de France, Laboratoire National d’Hydraulique et Environnement. 
