Volume of Fluid
--------------------

In this subsection, the parameters for multiphase flow simulation using the volume of fluid method (VOF) are specified. 

In this method, the two fluids considered are given index of :math:`0` and :math:`1` respectively. The amount of fluid at any given quadrature point is represented by a phase fraction between :math:`0` and :math:`1`. The interface is therefore considered located where the phase fraction :math:`= 0.5`. The interface between the two fluids is moved by a transport equation on the phase fraction.

.. note::

  At the moment, a maximum of two fluids is supported. By convention, air is usually the ``fluid 0`` and the other fluid of interest is the ``fluid 1``.

The default values of the VOF parameters are given in the text box below.

.. code-block:: text

	subsection VOF	
		subsection interface sharpening
			set enable 	= false
			set verbosity 	= quiet
			set sharpening threshold   = 0.5
			set interface sharpness    = 2
			set sharpening frequency   = 10
		end

		subsection peeling wetting
			set enable 	= false
			set verbosity 	= quiet
			set peeling pressure value 	= -0.05
			set peeling pressure gradient 	= -1e-3
			set wetting pressure value 	= 0.05
			set wetting phase distance 	= 0
			set diffusivity = 0
		end

		subsection mass conservation
			set skip mass conservation in fluid 0 = false
			set skip mass conservation in fluid 1 = false
			set monitoring 		= false
			set fluid monitored 	= 1
		end

		subsection surface tension force
			set enable 	= false
			set verbosity 	= quiet
			set output auxiliary fields 	= false
			set surface tension coefficient = 0.0
			set phase fraction gradient filter = 0.5
			set curvature filter 		= 0.5			
		end

	end

.. warning::
  Do not forget to ``set VOF = true`` in the :doc:`multiphysics` subsection of the ``.prm``.


.. seealso::
  See :doc:`initial_conditions` for the definition of the VOF initial conditions and `Physical properties - two phase simulations <https://lethe-cfd.github.io/lethe/parameters/cfd/physical_properties.html#two-phase-simulations>`_ for the definition of the physical properties of both fluids.


* ``subsection interface sharpening``: defines parameters to counter numerical diffusion of the VOF method and to avoid the interface between the two fluids becoming more and more blurry after each time step.

  * ``enable``: controls if interface sharpening is enabled.
  * ``verbosity``: enables the display of the residual at each non-linear iteration, to monitor the progress of the linear iterations, similarly to the ``verbosity`` option in :doc:`linear_solver_control`. Choices are: ``quiet`` (default, no output) and ``verbose``.
  * ``sharpening threshold``: phase fraction threshold, between 0 and 1. It is by default equal to :math:`0.5`, the phase fraction at which the interphase is considered located, but its value can be changed to :ref:`counter fluid loss or creation <tip on sharpening threshold>`.
  * ``interface sharpness``: sharpness of the moving interface (parameter :math:`a` in the `interface sharpening model <https://www.researchgate.net/publication/287118331_Development_of_efficient_interface_sharpening_procedure_for_viscous_incompressible_flows>`_).
  
  .. tip::
    This parameter must be larger than 1 for interface sharpening. Choosing values less than 1 leads to interface smoothing instead of sharpening. A good value would be between 1 and 2.
  
  * ``sharpening frequency``: sets the frequency (in number of iterations) for the interface sharpening computation.

.. seealso::

  The :doc:`../../examples/multiphysics/dam-break-VOF/dam-break-VOF` example using VOF represents well the interface sharpening issue.

* ``subsection peeling wetting``: Peeling and wetting mechanisms are very important to consider when there are solid boundaries in the domain, like a wall. If the fluid is already on the wall and its velocity drives it away from it, the fluid should be able to detach from the wall, meaning to `peel` from it. If the fluid is not already on the wall and its velocity drives it toward it, the fluid should be able to attach to the wall, meaning to `wet` it. This subsection defines the parameters for peeling and wetting mechanisms at the VOF boundaries, as defined in :doc:`boundary_conditions_multiphysics`. 

  * ``enable``: controls if peeling/wetting mechanism is enabled.
  * ``verbosity``: enables the display of the number of peeled and wet cells at each time-step. Choices are: ``quiet`` (default, no output) and ``verbose``.

  .. admonition:: Example of a ``set verbosity = verbose`` output:
  
    .. code-block:: text

      Peeling/wetting correction at step 2
        -number of wet cells: 24
        -number of peeled cells: 1

  * Peeling of the higher density fluid occurs where those conditions are met:

    * the cell is in the domain of the higher density fluid,
    * the cell pressure value is below ``peeling pressure value``, and
    * more than half of the quadrature points in the cell have a pressure gradient below ``peeling pressure gradient``.

    The cell is then filled with the lower density fluid by changing its phase value.

  * Wetting of the lower density fluid occurs where those conditions are met: 

    * the cell is in the domain of the lower density fluid,
    * the cell pressure value is above ``wetting pressure value``, and
    * the distance (on the phase value) to the interface is above ``wetting phase distance``.

    The cell is then filled with the higher density fluid by changing its phase value.

    .. tip::

      For ``set wetting phase distance = 0``, the wetting can only occur at the interface (considered at ``phase value = 0.5``).

      For ``set wetting phase distance`` :math:`> 0`, the wetting can occur in the area where is larger than the area occupied by the higher density fluid. For example:

      * if the ``fluid 1`` has a higher density than ``fluid 0``, and ``set wetting phase distance = 0.1``, the wetting can occur where the phase value is below :math:`= 0.4`.
      * if the ``fluid 0`` has a higher density than ``fluid 1``, and ``set wetting phase distance = 0.1``, the wetting can occur where the phase value is above :math:`= 0.6`.

    * ``diffusivity``: value of the diffusivity (diffusion coefficient) in the transport equation of the phase fraction. Default value is 0 to have pure advection. This can be used to :ref:`improve wetting`.

* ``subsection mass conservation``: By default, mass conservation (continuity) equations are solved on the whole domain, i.e. on both fluids. However, replacing the mass conservation by a zero-pressure condition on one of the fluid (typically, the air), so that it can get in and out of the domain, can be useful to :ref:`improve wetting`. This subsection defines parameters that can be used to skip mass conservation in one of the fluid, and to monitor the surface/volume (2D/3D) occupied by the other fluid of interest.

  * mass conservation can be skipped on the fluid with index 0 or 1, as defined in the subsection `Physical properties - two phase simulations <https://lethe-cfd.github.io/lethe/parameters/cfd/physical_properties.html#two-phase-simulations>`_, with ``skip mass conservation in fluid 0`` and ``skip mass conservation in fluid 1`` respectively.
  * ``monitoring``: controls if conservation is monitored at each iteration, through the volume computation of the fluid with index ``fluid monitored``. Results are outputted in a data table (`VOF_monitoring_fluid_0.dat` or `VOF_monitoring_fluid_1.dat`).

  .. admonition:: Example of file output, `VOF_monitoring_fluid_1.dat`:

    The ``volume_fluid_1`` column gives the surface/volume (2D/3D) occupied by the fluid with index 1.
  
    .. code-block:: text

	 time   volume_fluid_1 
	 0.0000     4.9067e-01 
	 0.0100     4.8425e-01 
	 0.0200     4.8564e-01 
	 0.0300     4.7858e-01 
	 0.0400     4.8245e-01 
	 0.0500     4.7693e-01 
	 0.0600     4.8154e-01 
	 0.0700     4.7590e-01 
	 0.0800     4.8133e-01 
	 0.0900     4.7604e-01 
	 0.1000     4.8198e-01 

.. _tip on sharpening threshold:

.. tip::
  Due to numerical diffusion of the interface, the ``peeling wetting`` mechanism or an added ``diffusivity``, the method used is not strictly conservative at every iteration. The ``sharpening threshold`` can then be adapted to counter fluid loss (e.g. ``set sharpening threshold = 0.45``) or creation (e.g. ``set sharpening threshold = 0.55``).


* ``subsection surface tension force``: Surface tension is the tendency of a liquid to maintain the minimum possible surface area. This subsection defines parameters to ensure an accurate interface between the two phases, used when at least one phase is liquid. 

  * ``enable``: controls if ``surface tension force`` is considered.
  * ``verbosity``: enables the display of the output from the surface tension force calculations. Choices are: ``quiet`` (default, no output) and ``verbose``.
  * ``output auxiliary fields``: enables the display of the filtered ``phase fraction gradient`` and filtered ``curvature``. Used for debugging purposes.
  * ``surface tension coefficient``: surface tension coefficient in :math:`Nm^{-1}`, as used to define the Weber number (:math:`We`):

    .. math::
        We = Re \cdot \frac{\mu_\text{ref} \; u_\text{ref}}{\sigma} 

    where :math:`Re` is the Reynolds number, :math:`\mu_\text{ref}` and :math:`u_\text{ref}` are some reference viscosity and velocity characterizing the flow problem, and :math:`\sigma` is the surface tension coefficient.

  * ``phase fraction gradient filter``: value used to apply a `projection step <https://onlinelibrary.wiley.com/doi/full/10.1002/fld.2643>`_ to damp high frequency errors, that are magnified by differentiation, in the phase fraction gradient (:math:`\bf{\psi}`), following the equation:

    .. math::
        \int_\Omega \left( {\bf{v}} \cdot {\bf{\psi}} + \eta_n \nabla {\bf{v}} \cdot \nabla {\bf{\psi}} \right) d\Omega = \int_\Omega \left( {\bf{v}} \cdot \nabla {\phi} \right) d\Omega

    where :math:`\bf{v}` is a piecewise continuous vector-valued test function, :math:`\bf{\psi}` is the filtered phase fraction gradient, :math:`\eta_n \geq 0` is the ``phase fraction gradient filter`` value, and :math:`\phi` is the phase fraction.

  .. tip::

    The ``phase fraction gradient filter`` must be a small value larger than 0. Use the procedure suggested in: :ref:`choosing values for the surface tension force filters`.

  * ``curvature filter``: value used to apply a `projection step <https://onlinelibrary.wiley.com/doi/full/10.1002/fld.2643>`_ to damp high frequency errors, that are magnified by differentiation, in the curvature (:math:`k`), following the equation:

    .. math:: 
        \int_\Omega \left( v k + \eta_k \nabla v \cdot \nabla k \right) d\Omega = \int_\Omega \left( \nabla v \cdot \frac{\bf{\psi}}{|\bf{\psi}|} \right) d\Omega

    where :math:`v` is a test function, :math:`k` is the filtered curvature, :math:`\eta_k` is the ``curvature filter`` value, and :math:`\bf{\psi}` is the filtered phase fraction gradient. 

  .. tip::

    Use the procedure suggested in: :ref:`choosing values for the surface tension force filters`.




.. _improve wetting:

Improving the Wetting mechanism
+++++++++++++++++++++++++++++++++++

In the framework of incompressible fluids, a layer of the lowest density fluid (e.g. air) can form between the highest density fluid (e.g. water) and the boundary, preventing its wetting. Two strategies can be used to improve the wetting mechanism:

1. Add a small ``diffusivity`` to the transport equation (e.g. ``set diffusivity = 1e-3``), so that the higher density fluid spreads to the boundary location. 

.. tip::
  It is strongly advised to sharpen the interface more often (e.g. ``set sharpening frequency = 2``) to limit interface blurriness due the added diffusivity. As peeling-wetting is handled after the transport equation is solved, but before interface sharpening, this will not prevent the wetting from occuring.

2. Remove the conservation condition on the lowest density fluid (e.g. ``set skip mass conservation in fluid 0 = false``). The mass conservation equation in the cells of interest is replaced by a zero-pressure condition, to allow the fluid to get out of the domain. 

.. tip::
  This can give more precise results as the interface remains sharp, but the time step (in :doc:`simulation_control`) must be low enough to prevent numerical instabilities.


.. _choosing values for the surface tension force filters:

Choosing values for the surface tension force filters
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

The following procedure is recommended to choose proper values for the ``phase fraction gradient filter`` and ``curvature filter``: 

1. Use ``set output auxiliary fields = true``.
2. Choose a small value, still larger than :math:`0`, for example :math:`h/10` with :math:`h` the smallest mesh size.
3. Run the simulation and check whether the filtered phase fraction gradient field is smooth and without oscillation.
4. If the filtered field (``phase fraction gradient`` or ``curvature``) shows oscillations, increase the value, for example :math:`h/5`, and repeat this process until reaching a smooth filtered field without oscillations.

