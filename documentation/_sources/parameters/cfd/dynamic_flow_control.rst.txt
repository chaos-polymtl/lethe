====================
Dynamic Flow Control
====================

The purpose of this subsection is to enable dynamic flow control. It is important when we want to simulate an average
velocity on a specific boundary (CFD) or the whole domain (for CFD-DEM). To control the average velocity of the flow, the code
calculates a :math:`\beta` coefficient at each time step that is used as a source term in the momentum equation to keep the average velocity at a targeted value.

The main controller of the average velocity is the following equation and is based on approach of Wang [#wang2023]_:

.. math::
    \beta^{n+1} = \beta^n + \frac{\alpha}{\Delta t} \left[ \bar{U}^{0} - 2\bar{U}^{n} + \bar{U}^{n-1} \right]

The default parameters are:

.. code-block:: text

  subsection flow control
    set enable               = false
    set enable beta particle = false
    set average velocity     = 0.0
    set inlet boundary id    = 0
    set flow direction       = 0
    set initial beta         = 0.0
    set alpha                = 1.0
    set beta threshold       = 0.0
    set verbosity            = quiet
  end

* The ``enable`` parameter is set to ``true`` to enable an imposed average velocity.

* The ``enable beta particle`` parameter is set to ``true`` to apply a proportional beta force on particles:

.. math::
    \beta^{n+1}_{p} = \frac{\rho}{\rho_{p}}\beta^{n+1}

* The ``average velocity`` parameter specifies the targeted average velocity (:math:`m/s`). The value is compared to the calculated value at a boundary (for CFD solvers) or the whole domain (for CFD-DEM solver).

.. note::

  If you are simulating a problem with periodic boundary conditions, it is important to control the average velocity by specifying the boundary id that corresponds to the inlet only.

* The ``inlet boundary id`` parameter is the boundary where the average velocity intended is controlled (for CFD solvers only).

* The ``flow direction`` parameter sets the absolute direction. It should be ``0`` if the flow is in the :math:`x` direction, ``1`` in the :math:`y` direction and ``2`` in the :math:`z` direction.

* The ``initial beta`` parameter sets an initial :math:`\beta` coefficient that may speed up the convergence of the average velocity target. Some tests should be done to find a good initial :math:`\beta` coefficient.

.. tip:: 

  A good method to find a reasonable initial beta is to test two or three different initial beta parameters, write down the given flow rate at the first time step in the simulation and do a regression. The correlation is linear and giving a proper value will reduce the time it takes to reach the target velocity.

* The ``alpha`` parameter is a relaxation coefficient that is used to control the convergence speed or stability. The higher the value, the faster the convergence. However, if the value is too high, there might be oscillations and the simulation could be unstable.

* The ``beta threshold`` parameter is the threshold on beta calculated at the previous time step that prevents the new calculated beta force to be updated. If ``abs(beta_n - beta_n+1) < abs(beta_n * beta_threshold)``, the previous beta value is kept. This could avoid the reassembly of the matrix because of the updated force term when reuse matrix option is enabled for the non-linear solver.

Reference
---------
.. [#wang2023] \W. Wang, “A non-body conformal grid method for simulations of laminar and turbulent flows with a compressible large eddy simulation solver,” Ph.D., Ann Arbor, United States, 2009. Accessed: May 4, 2023. [Online]. Available: https://www.proquest.com/docview/304905306\.
