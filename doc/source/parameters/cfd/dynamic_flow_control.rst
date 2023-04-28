Dynamic flow control
~~~~~~~~~~~~~~~~~~~~

The purpose of this subsection is to enable dynamic flow control. It is important when we want to simulate an average
velocity on a specific boundary (CFD) or the whole domain (for CFD-DEM). To control the average velocity of the flow, the code
calculates a :math:`\beta` coefficient at each time step that is used to keep the average velocity in an acceptable
threshold. If the boundary on which an average velocity is imposed is an inlet, the average velocity
targeted must be negative since the outward normal vector is used for the calculation. If the 
boundary is an outlet, average velocity must be positive. In other to keep the sign consistency, the targeted
velocity should have an negative sign when used with CFD-DEM even if no normal vector is used since the
whole domain is taken into account.

The main controller of the average velocity is the following equation and is based on approach of Wang [1]:

.. math::
    \beta^{n+1} = \beta^n - \frac{\alpha}{\Delta t} \left[ (\bar{U})^{0} - 2(\bar{U})^{n} + (\bar{U})^{n-1} \right]

The default parameters are:

.. code-block:: text

  subsection flow control
    set enable                  = false
    set enable beta particle    = false
    set average velocity target = 0
    set boundary id             = 0
    set flow direction          = 0
    set initial beta            = 0
    set alpha                   = 1.0
    set verbosity               = quiet
  end

* The ``enable`` parameter is set to ``true`` to enable an imposed average velocity.

* The ``enable beta particle`` parameter is set to ``true`` to apply a proportional beta force on particles.

* The ``average velocity target`` parameter specifies the average velocity (:math:`m/s`). The value is compared to the calculated value at a boundary (for CFD solvers) or the whole domain (for CFD-DEM solver).

.. warning::

  This value should be negative if imposed on an the inlet boundary (or CFD-DEM), and positive if imposed on an outlet boundary.

.. note::

  If you are simulating a problem with periodic boundary conditions, it is important to control the average velocity by specifying the boundary id that corresponds to the inlet only (where the flow would have a negative value).

* The ``boundary id`` parameter is the boundary where the average velocity intended is controlled (for CFD solvers only).

* The ``flow direction`` parameter sets the absolute direction. It should be ``0`` if the flow is in the :math:`x` direction, ``1`` in the :math:`y` direction and ``2`` in the :math:`z` direction.

* The ``initial beta`` parameter sets an initial :math:`\beta` coefficient that may speed up the convergence of the average velocity target. Some tests should be done to find a good initial :math:`\beta` coefficient.

.. tip:: 

  A good method to find a reasonable initial beta is to test two or three different initial beta parameters, write down the given flow rate at the first time step in the simulation and do a regression. The correlation is linear and giving a proper value will reduce the time it takes to reach the target velocity.

* The ``alpha`` parameter is a relaxation coefficient that is used to control the convergence speed or stability. The higher the value, the faster the convergence. However, if the value is too high, the convergence may be unstable.

Reference
---------
`[1] <https://www.proquest.com/dissertations-theses/non-body-conformal-grid-method-simulations/docview/304905306/se-2>`_ Wang, W. (2009). A non-body conformal grid method for simulations of laminar and turbulent flows with a compressible large eddy simulation solver. ProQuest Dissertations & Theses Global (Order No. 3355753), p.41.
