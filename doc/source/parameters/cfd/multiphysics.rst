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
    set fluid index 		= 0
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
    * ``conservation monitoring``: controls if conservation is monitored at each iteration, through the volume computation of the fluid with the given ``fluid index``. Results are outputted in a data table (`free_surface_monitoring.dat`).


.. warning::

  At the moment, a maximum of two fluids is supported, with respective ``fluid index`` of ``0`` and ``1``.

.. seealso::

  The VOF solver is used in the example :doc:`../../examples/multiphysics/dam-break-VOF/dam-break-VOF`.


