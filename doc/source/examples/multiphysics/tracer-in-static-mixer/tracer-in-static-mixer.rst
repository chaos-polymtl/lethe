======================================
Tracer in Static Mixer
======================================

In this example, we extend the :doc:`../../sharp-immersed-boundary/3d-rbf-static-mixer/3d-rbf-static-mixer` example using a passive tracer to quantify the mixing profile.

We inject two phases through a static mixer to verify how geometry affects their mixing.

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
- Tracer physics
- Accelerated transient auxiliary physics resolution with ``average_velocity_profile`` initial condition. 


---------------------------
Files Used in this Example
---------------------------

All mentioned files are found in the example's folder (``examples/multiphysics/tracer-in-static-mixer``).

* Parameter file for solving the flow: ``flow_in_long_mixer.prm``;
* Parameter file for solving the tracer transport: ``flow_in_long_mixer_only_tracer.prm``;
* Composite geometry file: ``mixer_long.composite``;
* RBF geometry file: ``RBF_helix.output``. The extension is ``.output`` because it was named from a `bitpit <https://github.com/optimad/bitpit>`_ perspective. RBF preparation was done in the previous :doc:`../../sharp-immersed-boundary/3d-rbf-static-mixer/3d-rbf-static-mixer` example, based on this surface grid from Thingiverse [#thingiverse]_ under CC BY 4.0;
* Python post-processing script: ``postprocess_tracer.py``.

-----------------------
Description of the Case
-----------------------

In this example, we use the auxiliary tracer physics to quantify the mixing of two fluids in a static mixer. We assume the flow is steady. Hence, we use a two-step simulation, consisting of:

#. Simulating the transient flow with a small ``time step`` until the flow has reached the steady state;
#. Simulate the flow of the tracer with a larger ``time step`` while reusing the computed velocity field from the previous simulation.

This strategy reduces the computational cost of the simulation as the fluid flow solution is reused in the second step.

--------------
Parameter File
--------------

Simulation 1: Flow resolution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simulation control
******************

.. code-block:: text

    subsection simulation control
      set method    = bdf1
      set time end  = 40e-4
      set time step = 1e-4

      set adapt = false

      set output path      = ./output/
      set output name      = output
      set output control   = iteration
      set output frequency = 5
    end

We use a ``bdf1`` time integration scheme with a short ``time step = 1e-4`` :math:`\text{s}`. After ``time end = 40e-4`` :math:`\text{s}`, we consider that the velocity field has reached a steady state. According to :doc:`../../sharp-immersed-boundary/3d-rbf-static-mixer/3d-rbf-static-mixer`, the pressure drop varies by less than :math:`5` % in the range [30, 40] :math:`\text{s}`.

Restart
******************

We use Lethe's checkpoint/restart mechanism to feed the steady-state flow information to the second step of the simulation.

.. code-block:: text

    subsection restart
      set checkpoint = true
      set frequency  = 5
      set filename   = restart
      set restart    = false
    end

We ``checkpoint`` the simulation at every ``5`` time steps, both for safety (each time step takes a long time to complete) and for reuse in the second simulation.

Multiphysics
******************

Both ``fluid dynamics`` and ``tracer`` are enabled for the first simulation to ensure initialization of the ``tracer`` field. However, the tracer injection only begins at the second step.

.. code-block:: text

    subsection multiphysics
      set fluid dynamics = true
      set tracer         = true
    end

Physical Properties
*******************

In this case we consider that we have a passive tracer in water. The units used for the presented case are :math:`\text{cm}`, :math:`\text{s}`, and :math:`\text{g}`.

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

#. The ``tracer diffusivity model`` is ``immersed solid tanh``. This model is used in ``lethe-fluid-sharp`` for tracer flow percolating immersed solids.
#. The ``tracer diffusivity outside`` is ``1e-5`` :math:`\text{cm²/s}`, as this is a typical value for a passive tracer in a liquid.
#. The ``tracer diffusivity inside`` is set to ``1e-10`` :math:`\text{cm²/s}`. The low value prevents diffusivity inside the solid while providing numerical stability (:math:`> 0`).
#. The ``thickness`` is ``5e-1`` :math:`\text{cm}`. At the scale of the problem, this provides a smooth transition without generating oscillations between liquid and solid phases. The thickness is of the order of magnitude of the smallest cell length to restrict the transition to one cell thickness.

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
          set Function expression = if(y<0,t<11?0:(t<61?1:0),0)
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

#.  We use ``time dependent`` boundary conditions, defined at :math:`y<0` and :math:`11 < t < 61` to inject a pulse on the lower half of the inlet.
#. All other boundary conditions are ``outlet``. This condition is natural for the outlet of the problem. For lateral walls, this condition represents an impermeable wall since velocity is perpendicular.

Post-processing
*******************

.. code-block:: text

    subsection post-processing
      set verbosity = verbose

      set calculate average velocities = true
      set initial time                 = 30e-4

      # Tracer post-processing
      set calculate tracer flow rate  = true
      set tracer flow rate name       = tracer_flow_rate
    end

#. ``calculate average velocities`` is enabled, beginning at ``initial time = 30e-4`` :math:`\text{s}`. This means that the last :math:`25\%` of the first simulation will be used to calculate the time-averaged fluid velocity profile.
#. ``calculate tracer flow rate`` is enabled to provide data for mixing quantification.

Simulation 2: Tracer transport
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the second simulation, we reuse the velocity profile from the first simulation and simply transport the passive tracer. This allows us to increase the time step and duration of the simulation while keeping realistic velocity and pressure solutions.

Simulation control
******************

.. code-block:: text

    subsection simulation control
      set method    = bdf1
      set time end  = 500
      set time step = 10
      ...
    end

#. We use a longer ``time step`` of ``10`` :math:`\text{s}`.
#. We simulate until ``time end = 500`` :math:`\text{s}` to allow the tracer to flow through the entire length of the domain.

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

We disable ``fluid dynamics``, as we use the time-averaged velocity profile from the previous simulation.

.. code-block:: text

    subsection multiphysics
      set fluid dynamics = false
      set tracer         = true
    end

Initial conditions
******************

We use ``average_velocity_profile`` from the previous simulation as the initial (and persistent) condition for the fluid dynamics. Since this physics is disabled, the velocity and pressure profiles will remain fixed throughout the simulation.

.. code-block:: text

    subsection initial conditions
      set type = average_velocity_profile
    end


----------------------
Running the Simulation
----------------------

As previously mentionned, the case is run in two steps:

#. Simulation to reach a pseudo steady-state of the flow field;
#. Transient simulation to transport tracer through the domain.

The simulation can be launched on multiple cores using ``mpirun`` and the ``lethe-fluid-sharp`` executable. Using 6 CPU cores and assuming that the ``lethe-fluid-sharp`` executable is within your path, the simulation can be launched by typing:

.. code-block:: text
  :class: copy-button

  mpirun -np 6 lethe-fluid-sharp flow_in_long_mixer.prm
  mpirun -np 6 lethe-fluid-sharp flow_in_long_mixer_only_tracer.prm


-------
Results
-------

The following movie shows the flow of the tracer through the static mixer, both as a colored slice and colored streamlines:

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/dSBRXLukX8E" frameborder="0" allowfullscreen></iframe>


The tracer evolution through the inlet and outlet can be monitored by plotting their values in time, extracted from ``/output/tracer_flow_rate.dat``. The Python script ``postprocess_tracer.py`` generates the following plot:

.. code-block:: text
  :class: copy-button

    python3 postprocess_tracer.py


+-----------------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/tracer_flow_rates.svg                                                                                   |
|     :align: center                                                                                                          |
|     :width: 600                                                                                                             |
|     :name: Tracer flow rates                                                                                                |
|                                                                                                                             |
+-----------------------------------------------------------------------------------------------------------------------------+

As the Figure shows, the concentration of the tracer flattens as it flows. The gap between between the inlet and outlet peaks is of :math:`170` :math:`\text{s}`. When compared to the theoretical time of :math:`150` :math:`\text{s}` (:math:`d_x/u_x`, with :math:`d_x` the domain length and :math:`u_x` the inlet velocity), this difference can be explained by retention effects and the tortuous paths that the tracer travels through.

---------
Reference
---------

.. [#thingiverse] "Group 9., Helix Static Mixer," *Thingiverse* Available: https://www.thingiverse.com/thing:3915237
