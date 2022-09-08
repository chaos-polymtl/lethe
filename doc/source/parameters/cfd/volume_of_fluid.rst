Multiphase Flow - Volume of Fluid
----------------------------------

In this subsection, the parameters for multiphase flow simulation using the volume of fluid method (VOF) are specified. 

In this method, the two fluids considered are given index of :math:`0` and :math:`1` respectively. The amount of fluid at any given quadrature point is represented by a phase fraction between :math:`0` and :math:`1`. The interface is therefore considered located where the phase fraction :math:`= 0.5`. The interface between the two fluids is moved by a transport equation on the phase fraction.

.. note::

  At the moment, a maximum of two fluids is supported. By convention, air is usually the ``fluid 0`` and the other fluid of interest is the ``fluid 1``.

The default values of the VOF parameters are given in the text box below.

.. code-block:: text

	subsection VOF	
		set viscous dissipative fluid = fluid 1
		set diffusivity = 0

		subsection interface sharpening
			set enable 	= false
			set frequency   = 10			
			set interface sharpness    = 2

			set type 	= constant

			# parameter for constant sharpening
			set threshold   = 0.5

			# parameters for adaptative sharpening
			set threshold max deviation = 0.20
			set max iterations = 20

			set verbosity 	= quiet
		end

		subsection peeling wetting
			set enable 	= false
			set verbosity 	= quiet
		end

		subsection mass conservation
			set conservative fluid  = both
			set monitoring 		= false
			set monitored fluid 	= fluid 1

			# parameters used with adaptative sharpening
			set tolerance		= 1e-6
			set verbosity 		= quiet
		end

		subsection surface tension force
			set enable 	= false
			set verbosity 	= quiet
			set output auxiliary fields 	   = false
			set surface tension coefficient    = 0.0
			set phase fraction gradient filter = 0.5
			set curvature filter 		   = 0.5	
            
			subsection marangoni effect
				set enable = false
				set surface tension gradient = 0.0
			end
		end

	end

.. warning::
  Do not forget to ``set VOF = true`` in the :doc:`multiphysics` subsection of the ``.prm``.


.. seealso::
  See :doc:`initial_conditions` for the definition of the VOF initial conditions and `Physical properties - two phase simulations <https://lethe-cfd.github.io/lethe/parameters/cfd/physical_properties.html#two-phase-simulations>`_ for the definition of the physical properties of both fluids.

* ``viscous dissipative fluid``: defines fluid(s) to which viscous dissipation is applied. 

  Choices are: ``fluid 0``, ``fluid 1`` (default) or ``both``, with the fluids defined in the :doc:`./physical_properties` for two phase simulations.

  .. warning::

	This parameter is used only if ``set heat transfer = true`` and ``set viscous dissipation = true`` in :doc:`./multiphysics`. 

  .. tip::

	Applying viscous dissipation in one of the fluids instead of both is particularly useful when one of the fluids is air. For numerical stability, the ``kinematic viscosity`` of the air is usually increased. However, but we do not want to have viscous dissipation in the air, because it would result in an unrealistic increase in its temperature.


* ``diffusivity``: value of the diffusivity (diffusion coefficient) in the transport equation of the phase fraction. Default value is ``0`` to have pure advection. Increase ``diffusivity`` to :ref:`improve wetting`.


  * ``subsection interface sharpening``: defines parameters to counter numerical diffusion of the VOF method and to avoid the interface between the two fluids becoming more and more blurry after each time step.

  * ``enable``: controls if interface sharpening is enabled.
  * ``frequency``: sets the frequency (in number of iterations) for the interface sharpening computation.
  * ``interface sharpness``: sharpness of the moving interface (parameter :math:`a` in the `interface sharpening model <https://www.researchgate.net/publication/287118331_Development_of_efficient_interface_sharpening_procedure_for_viscous_incompressible_flows>`_).
  
  .. tip::
    This parameter must be larger than 1 for interface sharpening. Choosing values less than 1 leads to interface smoothing instead of sharpening. A good value would be between 1 and 2.

  * ``type``: defines the interface sharpening type, either ``constant`` or ``adaptative``

    * ``set type = constant``: the sharpening ``threshold`` is the same throughout the simulation. This ``threshold``, between ``0`` and ``1`` (``0.5`` by default), corresponds to the phase fraction at which the interphase is considered located.
    * ``set type = adaptative``: the sharpening threshold is searched in the range :math:`\left[0.5-c_\text{dev} \; ; 0.5+c_\text{dev}\right]`, with :math:`c_\text{dev}` the ``threshold max deviation`` (``0.2`` by default), to ensure mass conservation. The search algorithm will stop either if the mass conservation ``tolerance`` is reached (see ``subsection mass conservation``), or if the number of search steps reach the number of ``max iterations``. If the ``tolerance`` is not reached, a warning message will be printed.

    .. warning::

      In case of adaptative interface sharpening (``set type = adaptative``), mass conservation must be monitored (``set monitoring = true`` in ``mass conservation`` subsection).

    .. admonition:: Example of a warning message if sharpening is adaptative but the mass conservation tolerance is not reached:
  
      .. code-block:: text

	  WARNING: Maximum number of iterations (5) reached in the 
	  adaptative sharpening threshold algorithm, remaining error
	  on mass conservation is: 0.02
	  Consider increasing the sharpening threshold range or the 
	  number of iterations to reach the mass conservation tolerance.

    .. tip::

      Usually the first iterations with sharpening are the most at risk to reach the ``max iterations`` without the ``tolerance`` being met, particularly if the mesh is quite coarse. 

      As most of the other iterations converge in only one step (corresponding to a final threshold of :math:`0.5`), increasing the sharpening search range through a higher ``threshold max deviation`` will relax the condition on the first iterations with a limited impact on the computational cost.

  * ``verbosity``: enables the display of the residual at each non-linear iteration, to monitor the progress of the linear iterations, similarly to the ``verbosity`` option in :doc:`linear_solver_control`. Choices are: ``quiet`` (default, no output), ``verbose`` (indicates sharpening steps) and ``extra verbose`` (details of the linear iterations).

    .. tip::
      
      The ``adaptive`` sharpening algorithm calls for the sharpening method multiple times to test different values of sharpening threshold. It is therefore advised to avoid using ``set verbosity = extra verbose`` in the ``subsection interface sharpening``.

.. seealso::

  The :doc:`../../examples/multiphysics/dam-break-VOF/dam-break-VOF` example using VOF represents well the interface sharpening issue.

* ``subsection peeling wetting``: Peeling and wetting mechanisms are very important to consider when there are solid boundaries in the domain, like a wall. If the fluid is already on the wall and its velocity drives it away from it, the fluid should be able to detach from the wall, meaning to `peel` from it. If the fluid is not already on the wall and its velocity drives it toward it, the fluid should be able to attach to the wall, meaning to `wet` it. This subsection defines the parameters for peeling and wetting mechanisms at the VOF boundaries, as defined in :doc:`boundary_conditions_multiphysics`. 

  .. important::
    This peeling/wetting mechanism implementation is an heuristic. It has been developed to meet the need of specific projects and gave satisfactory results as is, but it has not been broadly tested nor demonstrated, so its results should be considered with cautions. Do not hesitate to write to the team through the `Lethe github page <https://github.com/lethe-cfd/lethe>`_ would you need assistance.


  * ``enable``: controls if peeling/wetting mechanism is enabled.
  * ``verbosity``: enables the display of the number of peeled and wet cells at each time-step. Choices are: ``quiet`` (default, no output) and ``verbose``.

    .. admonition:: Example of a ``set verbosity = verbose`` output:
  
      .. code-block:: text

        Peeling/wetting correction at step 2
          -number of wet cells: 24
          -number of peeled cells: 1

  * Peeling occurs in a cell where the following conditions are met:

    * the cell is in the domain of the higher density fluid,
    * the cell pressure value is below the average pressure of the ``monitored fluid`` (``fluid 1`` by default, see ``subsection mass conservation``), and
    * the pressure gradient is negative for more than half of the quadrature points.

    The cell is then filled with the lower density fluid by changing its phase value progressively.

  * Wetting occurs in a cell where those conditions are met: 

    * the cell is in the domain of the lower density fluid,
    * the cell pressure value is above the average pressure of the ``monitored fluid`` (``fluid 1`` by default, see ``subsection mass conservation``), and
    * the pressure gradient is positive for more than half of the quadrature points.

    The cell is then filled with the higher density fluid by changing its phase value.

    .. tip::
      Even if ``monitoring`` is not enabled, the ``monitored fluid`` (``fluid 1`` by default) will be considered as the fluid of interest for the average pressure calculation in the peeling/wetting mechanism.

.. warning::

  As peeling/wetting mechanisms result in fluid generation and loss, is it highly advised to monitor the mass conservation of the fluid of interest (``subsection mass conservation``) and to change the type of sharpening threshold to adaptative (``subsection sharpening``).

* ``subsection mass conservation``: By default, mass conservation (continuity) equations are solved on the whole domain, i.e. on both fluids (``set conservative fluid = both``). However, replacing the mass conservation by a zero-pressure condition on one of the fluid (typically, the air), so that it can get in and out of the domain, can be useful to :ref:`improve wetting`. This subsection defines parameters that can be used to solve mass conservation in one fluid instead of both, and to monitor the surface/volume (2D/3D) occupied by the other fluid of interest.

  * ``conservative fluid``: defines fluid(s) to which conservation is solved. 

    Choices are: ``fluid 0``, ``fluid 1`` or ``both`` (default), with the fluids defined in the :doc:`./physical_properties` for two phase simulations.

  * ``monitoring``: controls if conservation is monitored at each iteration, through the volume computation of the fluid given as ``monitored fluid`` (``fluid 0`` or ``fluid 1`` (default)). Results are outputted in a data table (`VOF_monitoring_fluid_0.dat` or `VOF_monitoring_fluid_1.dat`).

    .. admonition:: Example of file output, `VOF_monitoring_fluid_1.dat`:

      The ``volume_fluid_1`` column gives the surface/volume (2D/3D) occupied by the fluid with index 1, its total mass, and the sharpening threshold used for this iteration.
  
      .. code-block:: text

	 time  volume_fluid_1 mass_fluid_1 sharpening_threshold 
	0.0000     4.9067e-01   3.8125e+02               0.5000 
	0.0050     4.9297e-01   3.8304e+02               0.5000 
	0.0100     4.9150e-01   3.8189e+02               0.5000 
	0.0150     4.9001e-01   3.8074e+02               0.5000 
	0.0200     4.8844e-01   3.7952e+02               0.5000 
	0.0250     4.9762e-01   3.8665e+02               0.5000 
	0.0300     4.9588e-01   3.8530e+02               0.5000 
	0.0350     4.9437e-01   3.8413e+02               0.5000 
	0.0400     4.9294e-01   3.8302e+02               0.5000 
	0.0450     4.9144e-01   3.8185e+02               0.5000 
	0.0500     5.0639e-01   3.9346e+02               0.5000 

  * ``tolerance``: value for the tolerance on the mass conservation of the monitored fluid, used with adaptative sharpening (see the ``subsection sharpening``). 
  
    For instance, with ``set tolerance = 0.02`` the sharpening threshold will be adapted so that the mass of the ``monitored fluid`` varies less than :math:`\pm 2\%` from the initial mass (at :math:`t = 0.0` sec).

  * ``verbosity``: states whether from the mass conservation data should be printed. Choices are quiet (no output), verbose (output information from the ``adaptive`` sharpening threshold) and extra verbose (output of the monitoring table in the terminal at the end of the simulation).

    .. admonition:: Example of mass conservation verbosity output (``verbose`` or ``extra verbose``):

      .. code-block:: text

	Sharpening interface at step 2
	   Adapting the sharpening threshold
	   ... step 1 of the search algorithm, min, avg, max mass deviation is : -0.1 -0.05 0.05
	   ... step 1 of the search algorithm, min, avg, max mass deviation is : -0.05 -0.025 0.04
	   ... search algorithm took : 2 step(s) 
	   ... error on mass conservation reached: -0.03
	   ... final sharpening is : 0.458224


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

  * ``subsection marangoni effect``: Marangoni effect is a thermocapillary effect, considered in simulations if ``set enable = true`` and if the ``surface tension gradient`` is not zero :math:`\left(\frac{\partial \sigma}{\partial T} \neq 0\right)`.

.. seealso::

  The surface tension force is used in the :doc:`../../examples/multiphysics/rising-bubble-VOF/rising-bubble-VOF` example.


.. _improve wetting:

Improving the Wetting mechanism
+++++++++++++++++++++++++++++++++++

In the framework of incompressible fluids, a layer of the lowest density fluid (e.g. air) can form between the highest density fluid (e.g. water) and the boundary, preventing its wetting. Two strategies can be used to improve the wetting mechanism:

1. Increase the ``diffusivity`` to the transport equation (e.g. ``set diffusivity = 1e-2``), so that the higher density fluid spreads even more to the boundary location. 

.. tip::
  It is strongly advised to sharpen the interface more often (e.g. ``set frequency = 2`` or even ``1``) to limit interface blurriness due the added diffusivity. As peeling-wetting is handled after the transport equation is solved, but before interface sharpening, sharpening will not prevent the wetting from occuring.

2. Remove the conservation condition on the lowest density fluid (e.g. ``set conservative fluid = fluid 1``). The mass conservation equation in the cells of interest is replaced by a zero-pressure condition, to allow the fluid to get out of the domain. 

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

