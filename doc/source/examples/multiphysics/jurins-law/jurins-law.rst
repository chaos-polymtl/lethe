===============================
Jurin's Law
===============================

This example simulates the capillary rise of a fluid between two planes caused by a constraint on the angle of contact between the wall and the fluid.


----------------------------------
Features
----------------------------------

- Solver: ``lethe-fluid`` (Q1-Q1)
- Two phase flow handled by the Cahn-Hilliard-Navier-Stokes (CHNS) approach
- Angle of contact boundary condition
- Dimensionality of the length
- Parametric sweep generation on the value of the angle of contact
- Post-processing of height difference between the outside fluid and the meniscus


--------------------------
Files Used in This Example
--------------------------

All files mentioned below are located in the example's folder (``examples/multiphysics/jurins-law``).

- Pointwise mesh file: ``jurins-law-2d-axisymmetric-dimensioned.pw``
- Mesh file: ``jurins-law-2d-mesh-dimensioned.msh``
- Parameter file: ``jurins-law-2d.prm``
- Postprocessing Python script: ``jurins_law_multiple_folders.py``


.. _Description of the case:

-------------------------
Description of the Case
-------------------------

Thanks to the symmetry of the problem, only one side is considered in the example. At :math:`t = 0`, a denser fluid (fluid 1) occupies the lower half of the domain with a lighter fluid (fluid 0) on the top. Because of the boundary condition imposing an angle of contact between the wall and the denser fluid, the surface is curved and a pressure gradient appears. Depending on the value of the angle, the height of the fluid will increase (or decrease) and reach an equilibrium height. The parameter file is ``jurins-law-2d.prm``
The computational domain is described in the following figure (not to scale):

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/jurins-law.svg                                                                                |
|     :alt: Computational domain definition and example of a capillary rise for an arbitrary (<90°) angle of contact|
|     :align: center                                                                                                |
|     :name: Computational domain and example of a capillary rise for an arbitrary (<90°) angle of contact          |
|                                                                                                                   |
|     Representation of the initial (left) and final (right) states of the capillary rise                           |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

The quantity of interest of this problem is the difference in height between the tip of the meniscus and the height of the fluid outside of the central part, denoted by :math:`\Delta H` in the previous figure. The Jurin's law [#liu2018]_ gives an asymptotical value of :math:`\Delta H` :

.. math::
    \Delta H = \frac{\sigma\cos{\alpha_c}}{\rho_0gR}

with :math:`\sigma` the surface tension coefficient, :math:`\alpha_c` the angle of contact imposed on the wall, :math:`\rho_1` the density of the denser fluid, :math:`g` the gravitational acceleration and :math:`R` the radius of the central part.

:math:`\Omega`
:math:`\Omega`
:math:`\Omega`
:math:`\Omega`
:math:`\Omega`


-----------------
Parameter File
-----------------

Simulation Control
~~~~~~~~~~~~~~~~~~

Time integration is handled by a 1st order backward differentiation scheme (`bdf1`), for a :math:`6 \ \text{s}` simulation time with an initial time step of :math:`0.001 \ \text{s}`. Time-step adaptation is enabled using ``adapt=true``
and the max CFL is :math:`0.8`.

.. note::
    This example uses an adaptive time-stepping method, where the time-steps are modified during the simulation to keep the maximum value of the CFL condition below the given threshold (0.5).

.. code-block:: text

    subsection simulation control
      set method                       = bdf1
      set time end                     = 6
      set time step                    = 0.001
      set adapt                        = true
      set max cfl                      = 0.8
      set output name                  = 3d-dam-break
      set output frequency             = 5
      set output path                  = ./output/
      set adaptative time step scaling = 1.05
      set output boundaries            = true
    end

Multiphysics
~~~~~~~~~~~~

The ``multiphysics`` subsection enables to turn on `(true)`
and off `(false)` the physics of interest. Here ``VOF`` is chosen.
Note that the fluid dynamics are solved by default.

.. code-block:: text

    subsection multiphysics
      set VOF = true
    end

Physical Properties
~~~~~~~~~~~~~~~~~~~

The ``physical properties`` subsection defines the physical properties of the fluids. In this example, we need two fluids with densities of :math:`1.204 \ \frac{\text{kg}}{\text{m}^3}` (air) and :math:`1000 \ \frac{\text{kg}}{\text{m}^3}` (water). However, the current numerical model was not able to solve with the real dynamic viscosities of the fluids. Therefore, they were altered in order to run the simulation.

.. warning::
    Altering the dynamic viscosities of the fluids will surely have an impact on the results. We will show this impact in the `<Results_>`_ section.

.. code-block:: text

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

Initial Conditions
~~~~~~~~~~~~~~~~~~

In the ``initial conditions`` subsection, we need to define the interface between the two fluids. We define this interface by using a function expression in the ``VOF`` subsection of ``initial conditions``. A projection step is applied to ensure a smooth definition of the initial condition.

.. code-block:: text

    subsection initial conditions
      set type = nodal
      subsection uvwp
        set Function expression = 0; 0; 0; 0
      end

      subsection VOF
        set Function expression = if (x>1.992 & z<0.55 & y>=-0.5, 1, 0)
        subsection projection step
          set enable           = true
          set diffusion factor = 1
        end
      end
    end

Source Term
~~~~~~~~~~~

In the ``source term`` subsection, we define the gravitational acceleration.

.. code-block:: text

    subsection source term
      subsection fluid dynamics
        set Function expression = 0;0;-9.81;0
      end
    end

VOF
~~~

In the ``VOF`` subsection, we select the ``tanh`` filter to filter the phase fraction and get a more defined interface. We set the value of beta to 10.

.. code-block:: text

    subsection VOF
      subsection phase filtration
        set type   = tanh
        set beta   = 10
      end
    end

Mesh
~~~~

In the ``mesh`` subsection, we specify the mesh used in this example. The structured mesh used in this example can be generated from the ``tank.geo`` file using `Gmsh <https://gmsh.info/#Download>`_. The initial refinement is set to :math:`3`.

.. code-block:: text

    subsection mesh
        set type                 = gmsh
        set file name            = tank.msh
        set initial refinement   = 3
    end


Mesh Adaptation
~~~~~~~~~~~~~~~

The ``mesh adaptation`` section controls the dynamic mesh adaptation. Here, we choose ``phase`` and ``pressure`` as the ``refinement variables``. The maximum and minimum refinement levels are respectively set to :math:`4` and :math:`2`.

.. code-block:: text

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
Running the Simulation
-----------------------

We call the lethe-fluid by invoking:

``mpirun -np $number_of_CPU_cores lethe-fluid 3d-dam-break.prm``

.. warning::
    Make sure to compile Lethe in `Release` mode and run in parallel using ``mpirun``. This simulation took :math:`\approx` 17 hours on 64 processes (runned on the `Narval <https://docs.alliancecan.ca/wiki/Narval/en>`_ cluster).

.. _Results:

-----------------
Results
-----------------

The following video shows the results of the simulation:

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/gaz4PiqhOzg"  frameborder="0" allowfullscreen></iframe>


As we can see, the simulated general evolution of the height seems to follow the experimentation results. However, on all 4 subplots, we notice that the height is overestimated. We also notice a slight shift to the right for :math:`H2`,  :math:`H3`, and :math:`H4` evolutions. These observations may be explained by the "highly viscous air" (fluid 0) that acts as an obstacle to the free flow of the water. Additionally, fluid 1 representing the water is 1000 times more viscous than regular water. With these results, we can see that the model needs to be improved to be able to accurately simulate low-viscosity fluids such as air. Furthermore, we observe that the wave formed at the impact with the obstacle doesn't collapse the right way due to the lack of compressibility of the air being simulated.


-----------
References
-----------


.. [#liu2018] \S. Liu, S. Li, and J. Liu, ‘Jurin’s law revisited: Exact meniscus shape and column height’, Eur. Phys. J. E, vol. 41, no. 3, p. 46, Mar. 2018, doi: 10.1140/epje/i2018-11648-1.

