Multiphysics
--------------
This subsection defines the multiphysics interface of Lethe and enables the solution of Auxiliary Physics in addition to traditional fluid dynamics simulations.

.. code-block:: text

  subsection multiphysics
    # Thermal physics
    set heat transfer 		= false
    set viscous dissipation 	= false
    set buoyancy force 		= false

    # Tracer
    set tracer 			= false

    # Multiphase flow
    set VOF 			= false
    set interface sharpening 	= false
    set conservation monitoring = false
    set fluid monitored		= 1
  end


* ``heat transfer``: controls if the heat transfer auxiliary physics are solved. This is an advection-diffusion equation. 

  When ``set heat transfer = true``, these optional parameters can be used:
   * ``viscous dissipation``: controls if the viscous dissipation is taken into account in the heat transfer equation.
   * ``buoyancy force``: controls if the buoyancy force is taken into account in the Navier-Stokes equations. The buoyancy force is calculated using the Boussinesq approximation.

.. seealso::

  The heat transfer solver is used in the example :doc:`../../examples/multiphysics/warming-up-a-viscous-fluid/warming-up-a-viscous-fluid`.

* The ``tracer`` parameter adds a passive tracer auxiliary physics. This is an advection-diffusion equation.

.. seealso::

  The tracer solver is used in the example :doc:`../../examples/multiphysics/tracer-through-cad-junction/tracer-through-cad-junction`.

* ``VOF``: enables multiphase flow simulations, with two fluids separated by a free surface. Volume-of-Fluid method is used, with the phase parameter following an advection equation. 

  See :doc:`initial_conditions` for the definition of the VOF conditions and `Physical properties - two phase simulations <https://lethe-cfd.github.io/lethe/parameters/cfd/physical_properties.html#two-phase-simulations>`_ for the definition of the physical properties of both fluids.

  When ``set VOF = true``, these optional parameters can be used:
    * ``interface sharpening``: controls if the interface sharpening method is used. Additional parameters can be found at :doc:`interface_sharpening`.
    * ``conservation monitoring``: controls if conservation is monitored at each iteration, through the volume computation of the ``fluid monitored``. Results are outputted in a data table (`VOF_monitoring_fluid_0.dat` or `VOF_monitoring_fluid_1.dat`).
    * ``fluid monitored``: index of the monitored fluid, if ``set conservation monitoring = true``. 

      .. important::

        Must match with the index set in `Physical properties - two phase simulations <https://lethe-cfd.github.io/lethe/parameters/cfd/physical_properties.html#two-phase-simulations>`_. For instance, if the physical properties for the fluid of interest is defined in ``subsection fluid 0``, use ``set fluid monitored = 0``.


.. warning::

  At the moment, a maximum of two fluids is supported, with respective index of ``0`` and ``1``. By convention, air is usually the ``fluid 0`` and the other fluid of interest is the ``fluid 1``.

.. seealso::

  The VOF solver is used in the example :doc:`../../examples/multiphysics/dam-break-VOF/dam-break-VOF`.


