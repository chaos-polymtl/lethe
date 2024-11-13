======================================
Tracer in Static Mixer
======================================

In this example, we extend the :doc:`../../sharp-immersed-boundary/3d-rbf-static-mixer/3d-rbf-static-mixer` example by using a passive tracer to study mixing effects more formally.

We inject two phases through a static mixer to verify the mixing effects of a given geometry.

+-----------------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/static_mixer_stl_casing_arrows.png                                                                      |
|     :align: center                                                                                                          |
|     :width: 800                                                                                                             |
|     :name: Surface grid representation of a helix static mixer with its casing.                                             |
|                                                                                                                             |
+-----------------------------------------------------------------------------------------------------------------------------+

----------------------------------
Features
----------------------------------

- Solver: ``lethe-fluid-sharp``
- Transient problem
- Displays the use of the tracer physics
- Usage of the ``average_velocity_profile`` initial condition to accelerate transient auxiliary physics resolution


---------------------------
Files Used in this Example
---------------------------

All files mentioned below are located in the example's folder (``examples/multiphysics/tracer-in-static-mixer``).

* Parameter file for solving the flow: ``flow_in_long_mixer.prm``;
* Parameter file for solving the tracer transport: ``flow_in_long_mixer_only_tracer.prm``;
* Composite geometry file: ``mixer_long.composite``;
* RBF geometry file: ``RBF_helix.output``. The extension is ``.output`` because it was named from a `bitpit <https://github.com/optimad/bitpit>`_ perspective. RBF preparation was done in the previous :doc:`../../sharp-immersed-boundary/3d-rbf-static-mixer/3d-rbf-static-mixer` example;
* Python post-processing script: ``postprocess_tracer.py``.

-----------------------
Description of the Case
-----------------------

In this example, we use the auxiliary tracer physics to  study the mixing effects of a static mixer wherein the flow is steady. To achieve this, two simulations are run:

#. Solving a pseudo-steady problem using a ``bdf1`` integration scheme with a small ``type step``, until a steady solution is reached;
#. Solving a transient tracer simulation with a large ``time step`` while reusing the computed velocity field from the previous simulation, which dramatically reduces computing time.

--------------
Parameter File
--------------

Flow resolution
~~~~~~~~~~~~~~~~~~~~~~~~~

Simulation control
******************

.. code-block:: text

    subsection simulation control
      set method    = bdf1
      set time end  = 40e-4
      set time step = 1e-4

      set adapt             = false
      set number mesh adapt = 0

      set output path      = ./output/
      set output name      = output
      set output control   = iteration
      set output frequency = 5
    end

#. We use a ``bdf1`` time integration scheme with a short ``time step = 1e-4``. After ``time end = 40e-4``, we consider that the velocity field has reached a steady-state.
#. We do not allow time step adaptation, since it is a parameter that we control ourselves between simulation: ``adapt = false``.

Restart
******************

We use the checkpoint/restart mechanism of Lethe to achieve these simulations.

.. code-block:: text

    subsection restart
      set checkpoint = true
      set frequency  = 5
      set filename   = restart
      set restart    = false
    end

#. We ``checkpoint`` the simulation at each ``5`` time steps, both for safety (each time step takes a long time to complete) and for reuse in the second simulation.
#. We do not ``restart`` the simulation, unless it is stopped before it reaches ``time end``.

Multiphysics
******************

Both ``fluid dynamics`` and ``tracer`` are enabed for the first simulation, although nothing happens tracer-wise.

.. code-block:: text

    subsection multiphysics
      set fluid dynamics = true
      set tracer         = true
    end

Physical Properties
*******************

The physical properties here are defined using centimeter length-unit and seconds time-unit, assuming that we have a passive tracer in water.

.. code-block:: text

    subsection physical properties
      subsection fluid 0
        set kinematic viscosity      = 0.01
        set tracer diffusivity model = immersed solid tanh
        subsection immersed solid tanh
          set tracer diffusivity inside  = 1e-10
          set tracer diffusivity outside = 1e-5
          set thickness                  = 5e-1
        end
      end
    end

#. The ``tracer diffusivity model`` is ``immersed solid tanh``. This model is used for tracer simulations involving immersed solids when using ``lethe-fluid-sharp``.
#. The ``tracer diffusivity outside`` is ``1e-5``, as this is a typical value for a passive tracer inside a liquid.
#. The ``tracer diffusivity inside`` is set to ``1e-10``. We desire no diffusivity inside the solid, and this value is low enough to achieve this effect while providing numerical stability.
#. The ``thickness`` is ``5e-1``. At the scale of the problem, this provides a smooth transition without generating oscillations between liquid and solid phases.

Tracer Boundary Conditions
***************************

.. code-block:: text

    subsection boundary conditions tracer
      set number         = 6
      set time dependent = true
      subsection bc 0
        set id   = 0
        set type = dirichlet
        subsection dirichlet
          set Function expression = if(y<0,t<0.1?0:(t<1.1?1:0),0)
        end
      end
      subsection bc 1
        set id   = 1
        set type = outlet
      end
      subsection bc 2
        set id   = 2
        set type = outlet
      end
      subsection bc 3
        set id   = 3
        set type = outlet
      end
      subsection bc 4
        set id   = 4
        set type = outlet
      end
      subsection bc 5
        set id   = 5
        set type = outlet
      end
    end

#. The boundary conditions are ``time dependent`` because of the inlet.
#. All other boundary conditions are ``outlet``. This condition is natural for the outlet of the problem. For lateral walls, this condition is also appropriate since no flow is occurring through these; it is then equivalent to a impermeable wall.

Post-processing
*******************

.. code-block:: text

    subsection post-processing
      set verbosity = verbose

      set calculate average velocities = true
      set initial time                 = 30e-4

      # Flow post-processing
      set calculate pressure drop = true
      set calculate flow rate     = true
      set inlet boundary id       = 0
      set outlet boundary id      = 1

      # Tracer post-processing
      set calculate tracer statistics = true
      set tracer statistics name      = tracer_statistics
      set calculate tracer flow rate  = true
      set tracer flow rate name       = tracer_flow_rate
    end

#. ``calculate average velocities`` is enabled, and its ``initial time = 30e-4``. This means that the last :math:`25\%` of the first simulation will be used to calculate the average velocity profile.
#. ``calculate tracer physics`` and ``calculate tracer flow rate`` are enabled to provide data for understanding mixing effects.

Tracer transport
~~~~~~~~~~~~~~~~~~~~~~~~~

In the second simulation, we reuse the velocity profile from the first simulation and simply transport the passive tracer. This allows to increase the time step and duration of the simulation while keeping numerical stability.

Simulation control
******************

.. code-block:: text

    subsection simulation control
      set method    = bdf1
      set time end  = 50
      set time step = 0.1
      ...
    end

#. We use a longer ``time step`` of ``0.1``.
#. We let the simulation run until ``time end = 50`` to allow the tracer to pass through.

Restart
******************

We enable the ``restart`` from the previous simulation.

.. code-block:: text

    subsection restart
      ...
      set restart    = true
    end

Multiphysics
******************

We disable ``fluid dynamics``, as we simply reuse the time-averaged velocity profile from the previous simulation.

.. code-block:: text

    subsection multiphysics
      set fluid dynamics = false
      set tracer         = true
    end

Initial conditions
******************

We use ``average_velocity_profile`` from the previous simulation as the initial (and constant) condition for the fluid dynamics. Since this physics is disable, the velocity and pressure profiles will remain constant throughout the simulation.

.. code-block:: text

    subsection initial conditions
      set type = average_velocity_profile
    end


----------------------
Running the Simulation
----------------------

As previously mentionned, the case is run in two steps:
#. Simulation to reach a pseudo steady-state for the flow field;
#. Transient simulation to transport tracer through the domain.

The simulation can be launched on multiple cores using ``mpirun`` and the ``lethe-fluid-sharp`` executable. Using 6 CPU cores and assuming that the ``lethe-fluid-sharp`` executable is within your path, the simulation can be launched by typing:

.. code-block:: text
  :class: copy-button

  mpirun -np 6 lethe-fluid-sharp flow_in_long_mixer.prm
  mpirun -np 6 lethe-fluid-sharp flow_in_long_mixer_only_tracer.prm


-------
Results
-------

The results in ``.pvd`` format can then be viewed using visualisation software such as Paraview. 

.. image:: images/paraview-tracer.png
    :alt: Simulation results in Meshgrid format
    :align: center

The higher presence of tracer in the outlet on the same side as the tracer inlet may indicate poor mixing.
As the tracer diffusivity is low, the mixing between the streams comes mainly from advection.
However, since the kinematic viscosity is high, the flow is laminar (i.e. dominated by viscous forces) and
the streamlines do not cross. 