============
Multiphysics
============

This subsection defines the multiphysics interface of Lethe and enables the solution of Auxiliary Physics in addition to traditional fluid dynamics simulations.

.. code-block:: text

  subsection multiphysics
    set fluid dynamics                  = true

    # Thermal physics
    set heat transfer                   = false
    set viscous dissipation             = false
    set thermal buoyancy force          = false

    # Tracer
    set tracer                          = false

    # Multiphase flow
    # Conservative Level-Set method
    set cls                             = false

    # Cahn-Hilliard equations
    set cahn hilliard                   = false

    # Electromagnetics
    set electromagnetics                = false
    set microwave heating               = false
  end



* ``fluid dynamics``: controls if the fluid dynamics are solved. This is ``true`` by default and can be turned to ``false`` to enable calculation of an auxiliary physic only. When appropriate, this can decrease drastically the computation time. 

.. tip::

  ``fluid dynamics = false`` and ``heat transfer = true`` enables to solve the heat transfer equations with a large time step, with the fluid velocity being defined from the :doc:`initial_conditions`.

* ``heat transfer``: controls if the heat transfer auxiliary physics are solved. This is an advection-diffusion equation. 

  When ``set heat transfer = true``, these optional parameters can be used:
   * ``viscous dissipation``: controls if the viscous dissipation is taken into account in the heat transfer equation.

   * ``thermal buoyancy force``: controls if the thermal buoyancy force is taken into account in the Navier-Stokes equations. The thermal buoyancy force is calculated using the Boussinesq approximation.

.. seealso::

  The heat transfer solver is used in the example :doc:`../../examples/multiphysics/warming-up-a-viscous-fluid/warming-up-a-viscous-fluid`.

* The ``tracer`` parameter adds a passive tracer auxiliary physics. This is an advection-diffusion equation.

.. seealso::

  The tracer solver is used in the example :doc:`../../examples/multiphysics/tracer-in-static-mixer/tracer-in-static-mixer`.

* ``CLS``: enables multiphase flow simulations, with two fluids separated by a free surface, using the Conservative Level-Set method. 

  See :doc:`conservative_level_set` for advanced CLS parameters, :doc:`initial_conditions` for the definition of the CLS conditions and `Physical properties - two phase simulations <https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/physical_properties.html#two-phase-simulations>`_ for the definition of the physical properties of both fluids.

.. seealso::

  The CLS solver is used in the example :doc:`../../examples/multiphysics/dam-break/dam-break`.

* The ``cahn hilliard`` parameter enables multiphase flow simulations, with two fluids separated by a free surface, using the Cahn-Hilliard equations. 

  See :doc:`cahn_hilliard` for advanced Cahn-Hilliard parameters, :doc:`initial_conditions` for the definition of the Cahn-Hilliard conditions and `Physical properties - two phase simulations <https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/physical_properties.html#two-phase-simulations>`_ for the definition of the physical properties of both fluids. 

* The ``electromagnetics`` parameter enables the solution of the time-harmonic Maxwell equations. 

  See :doc:`time_harmonic_maxwell` for advanced time-harmonic Maxwell parameters, :doc:`boundary_conditions_multiphysics` for the definition of the electromagnetic boundary conditions and :doc:`physical_properties` for the definition of the physical properties of the medium.

  When ``set heat transfer = true``, in addition to the electromagnetic solver, the optional parameter ``microwave heating`` can be used to enable the calculation of the heat source due to the electromagnetic fields. If ``set microwave heating = true``, the following heat source is calculated and added to the right-hand side of the heat transfer equation:

  .. math::

      Q_\text{em}=-\nabla \cdot \overline{\mathbf{S}} = \frac{1}{2}\sigma|\mathbf{E}|^2 + \frac{1}{2}\omega\varepsilon_0\varepsilon_\mathrm{im}|\mathbf{E}|^2 + \frac{1}{2}\omega\mu_0\mu_\mathrm{im}|\mathbf{H}|^2,

  where :math:`\overline{\mathbf{S}}` is the time-averaged Poynting vector, :math:`\sigma` is the conductivity, :math:`\varepsilon_0` is the vacuum permittivity, :math:`\varepsilon_\mathrm{im}` is the imaginary part of the relative permittivity, :math:`\mu_0` is the vacuum permeability, :math:`\mu_\mathrm{im}` is the imaginary part of the relative permeability, and :math:`\mathbf{E}` and :math:`\mathbf{H}` are the electric and magnetic fields, respectively.

.. seealso::

  The electromagnetic solver can be used on its own, an example is available in :doc:`../../examples/multiphysics/waveguide/waveguide`.