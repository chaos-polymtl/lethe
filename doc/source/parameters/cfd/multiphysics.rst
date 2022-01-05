Multiphysics
--------------
This subsection defines the multiphysics interface of Lethe and enables the solution of Auxiliary Physics in addition to traditional fluid dynamics simulations.

.. code-block:: text

  subsection multiphysics
    set heat transfer = false
    set buoyancy force = false
    set tracer = false
    set VOF = false
    set interface sharpening = false
    set viscous dissipation = false
  end


* The ``heat transfer`` parameter adds the heat transfer auxiliary physics. This is an advection-diffusion equation. If  ``viscous dissipation`` is set to true, viscous dissipation is taken into account in the heat transfer equation.

* The ``buoyancy force`` parameter adds the buoyancy force in the Navier-Stokes equations. The buoyancy force is calculated using the Boussinesq approximation. Note that ``heat transfer`` is a prerequisite for ``buoyancy force``, and it must be set to true when we set ``buoyancy force`` to true.

* The ``tracer`` parameter adds a passive tracer auxiliary physics. This is an advection-diffusion equation.

* The ``VOF`` parameter enables defining a multiphasic simulation, with 2 fluids separated by a free surface. Volume-of-Fluid method is used, with the phase parameter following an advection equation. See :doc:`initial_conditions` for the definition of the VOF conditions and :doc:`physical_properties` for the definition of the physical properties of both fluids. At the moment, a maximum of two fluids is supported.

* The ``interface sharpening`` parameter activates the interface sharpening method in VOF simulations.
