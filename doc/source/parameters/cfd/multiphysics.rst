Multiphysics
--------------
This subsection defines the multiphysics interface of Lethe and enables the solution of Auxiliary Physics in addition to traditional fluid dynamics simulations.

.. code-block:: text

  subsection multiphysics
    # Thermal physics
    set fluid mechanics 	= true
    set heat transfer 		= false
    set viscous dissipation 	= false
    set buoyancy force 		= false
    set use average velocity field = false

    # Tracer
    set tracer 			= false

    # Multiphase flow
    set VOF 			= false
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

* ``VOF``: enables multiphase flow simulations, with two fluids separated by a free surface, using the Volume-of-Fluid method. 

  See :doc:`volume_of_fluid` for advanced VOF parameters, :doc:`initial_conditions` for the definition of the VOF conditions and `Physical properties - two phase simulations <https://lethe-cfd.github.io/lethe/parameters/cfd/physical_properties.html#two-phase-simulations>`_ for the definition of the physical properties of both fluids.
  
* ``use average velocity field``:  controls if the advection of the subphysics is done using the average velocity field instead of the current velocity field. This is useful when the flow dynamics and the subphysics reach a pseudo-steady state at different time scales. To use this feature, the user should launch a simulation with the fluid mechanics solver and use the time averaging feature. Once the time average of the velocity field is sufficiently established, the simulation should be stopped and restarted without the fluid mechanics solver. The subphysics can then be solved using a larger time step.

.. note::
	Currently, this mechanism is only implemented in the heat transfer physic.
	
.. important::
   The subphysic must also be enabled (for example, ``set heat transfer = true``)
   
.. seealso::

  The VOF solver is used in the example :doc:`../../examples/multiphysics/dam-break-VOF/dam-break-VOF`.


