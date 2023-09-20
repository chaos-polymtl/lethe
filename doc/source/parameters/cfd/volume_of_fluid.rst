=================================
Volume of Fluid (Multiphase Flow)
=================================

In this subsection, the parameters for multiphase flow simulation using the volume of fluid method (VOF) are specified. 

In this method, the two fluids considered are given index of :math:`0` and :math:`1` respectively. The amount of fluid at any given quadrature point is represented by a phase fraction between :math:`0` and :math:`1`. The interface is therefore considered located where the phase fraction :math:`= 0.5`. The interface between the two fluids is moved by a transport equation on the phase fraction.

.. note::

  At the moment, a maximum of two fluids is supported. By convention, air is usually the ``fluid 0`` and the other fluid of interest is the ``fluid 1``.    See :doc:`initial_conditions` for the definition of the VOF initial conditions and :ref:`Physical properties - Two Phase Simulations<two phase simulations>` for the definition of the physical properties of both fluids.  Do not forget to ``set VOF = true`` in the :doc:`multiphysics` subsection of the ``.prm``.


The default values of the VOF parameters are given in the text box below.

.. code-block:: text

  subsection VOF

    set viscous dissipative fluid = fluid 1
    set diffusivity               = 0
    set compressible              = false

    subsection interface sharpening
      set enable                  = false
      set frequency               = 10
      set interface sharpness     = 2
      set type                    = constant

      # parameter for constant sharpening
      set threshold               = 0.5

      # parameters for adaptive sharpening
      set threshold max deviation = 0.20
      set max iterations          = 20

      set verbosity               = quiet
    end

    subsection phase filtration
      set type      = none
      set verbosity = quiet

      # parameter for the tanh filter
      set beta      = 20
    end

    subsection surface tension force
      set enable                                = false
      set verbosity                             = quiet
      set output auxiliary fields               = false
      set phase fraction gradient filter factor = 4
      set curvature filter factor               = 1
      set enable marangoni effect               = false
    end

    subsection mass conservation
      set monitoring         = false
      set monitored fluid    = fluid 1

      # parameters used with adaptive sharpening
      set tolerance          = 1e-6
      set verbosity          = quiet
    end
  end

* ``viscous dissipative fluid``: defines fluid(s) to which viscous dissipation is applied.

  Choices are: ``fluid 0``, ``fluid 1`` (default) or ``both``, with the fluid IDs defined in Physical properties - :ref:`two phase simulations`.

  .. tip::
    Applying viscous dissipation in one of the fluids instead of both is particularly useful when one of the fluids is air. For numerical stability, the ``kinematic viscosity`` of the air is usually increased. However, we do not want to have viscous dissipation in the air, because it would result in an unrealistic increase in its temperature. This parameter is used only if ``set heat transfer = true`` and ``set viscous dissipation = true`` in :doc:`./multiphysics`.

* ``diffusivity``: value of the diffusivity (diffusion coefficient) in the transport equation of the phase fraction. Default value is ``0`` to have pure advection. 
* ``compressible``: enables interface compression (:math:`\phi \nabla \cdot \mathbf{u}`) in the VOF equation.  This term should be kept to its default value of ``false`` except when compressible equations of state are used.


Interface Sharpening
~~~~~~~~~~~~~~~~~~~~~

* ``subsection interface sharpening``: defines parameters to counter numerical diffusion of the VOF method and to avoid the interface between the two fluids becoming more and more blurry after each time step. The reader is refered to the Interface sharpening section of :doc:`../../../theory/multiphysics/vof` theory guide for additional details on this sharpening method.

  * ``enable``: controls if interface sharpening is enabled.
  * ``frequency``: sets the frequency (in number of iterations) for the interface sharpening computation.
  * ``interface sharpness``: sharpness of the moving interface (parameter :math:`a` in the `interface sharpening model <https://www.researchgate.net/publication/287118331_Development_of_efficient_interface_sharpening_procedure_for_viscous_incompressible_flows>`_). This parameter must be larger than 1 for interface sharpening. Choosing values less than 1 leads to interface smoothing instead of sharpening. A good value would be around 1.5.

  * ``type``: defines the interface sharpening type, either ``constant`` or ``adaptive``

    * ``set type = constant``: the sharpening ``threshold`` is the same throughout the simulation. This ``threshold``, between ``0`` and ``1`` (``0.5`` by default), corresponds to the phase fraction at which the interface is located.
    * ``set type = adaptive``: the sharpening threshold is searched in the range :math:`\left[0.5-c_\text{dev} \; ; 0.5+c_\text{dev}\right]`, with :math:`c_\text{dev}` the ``threshold max deviation`` (``0.2`` by default), to ensure mass conservation. The search algorithm will stop either if the mass conservation ``tolerance`` is reached (see ``subsection mass conservation``), or if the number of search steps reach the number of ``max iterations``. If the ``tolerance`` is not reached, a warning message will be printed.

    .. warning::

      In case of adaptive interface sharpening (``set type = adaptive``), mass conservation must be monitored (``set monitoring = true`` in ``mass conservation`` subsection).

    .. admonition:: Example of a warning message if sharpening is adaptive but the mass conservation tolerance is not reached:

      .. code-block:: text

        WARNING: Maximum number of iterations (5) reached in the
        adaptive sharpening threshold algorithm, remaining error
        on mass conservation is: 0.02
        Consider increasing the sharpening threshold range or the
        number of iterations to reach the mass conservation tolerance.

    .. tip::

      Usually the first iterations with sharpening are the most at risk to reach the ``max iterations`` without the ``tolerance`` being met, particularly if the mesh is quite coarse.

      As most of the other iterations converge in only one step (corresponding to a final threshold of :math:`0.5`), increasing the sharpening search range through a higher ``threshold max deviation`` will relax the condition on the first iterations with a limited impact on the computational cost.

  * ``verbosity``: enables the display of the residual at each non-linear iteration, to monitor the progress of the linear iterations, similarly to the ``verbosity`` option in :doc:`linear_solver_control`. Choices are: ``quiet`` (default, no output), ``verbose`` (indicates sharpening steps) and ``extra verbose`` (details of the linear iterations).

  .. seealso::

    The :doc:`../../examples/multiphysics/dam-break/dam-break` example discussed the interface sharperning mechanism.


Phase Filtration
~~~~~~~~~~~~~~~~~~

* ``subsection phase filtration``: This subsection defines the filter applied to the phase fraction. This affects the definition of the interface.

* ``type``: defines the filter type, either ``none`` or ``tanh``

  * ``set type = none``: the phase fraction is not filtered
  * ``set type = tanh``: the filter function described in the Interface filtration section of :doc:`../../../theory/multiphysics/vof` theory guide is applied.
* ``beta``: value of the :math:`\beta` parameter of the ``tanh`` filter
* ``verbosity``: enables the display of filtered phase fraction values. Choices are ``quiet`` (no output) and ``verbose`` (displays values)


Surface Tension Force
~~~~~~~~~~~~~~~~~~~~~~

* ``subsection surface tension force``: Surface tension is the tendency of a liquid to maintain the minimum possible surface area. This subsection defines parameters to ensure an accurate interface between the two phases, used when at least one phase is liquid. 

  * ``enable``: controls if ``surface tension force`` is considered.

    .. attention::

      When the surface tension force is enabled, a ``fluid-fluid`` material interaction and a ``surface tension model`` with its parameters must be specified in the :doc:`physical_properties` subsection.

  * ``verbosity``: enables the display of the output from the surface tension force calculations. Choices are: ``quiet`` (default, no output) and ``verbose``.
  * ``output auxiliary fields``: enables the display of the filtered ``phase fraction gradient`` and filtered ``curvature``. Used for debugging purposes.

  * ``phase fraction gradient filter factor``: value of the factor :math:`\alpha` applied in the filter :math:`\eta_n = \alpha h^2`, where :math:`h` is the cell size. This filter is used to apply a `projection step <https://onlinelibrary.wiley.com/doi/full/10.1002/fld.2643>`_ to damp high frequency errors, that are magnified by differentiation, in the phase fraction gradient (:math:`\bf{\psi}`), following the equation:

    .. math::
        \int_\Omega \left( {\bf{v}} \cdot {\bf{\psi}} + \eta_n \nabla {\bf{v}} \cdot \nabla {\bf{\psi}} \right) d\Omega = \int_\Omega \left( {\bf{v}} \cdot \nabla {\phi} \right) d\Omega

    where :math:`\bf{v}` is a piecewise continuous vector-valued test function, :math:`\bf{\psi}` is the filtered phase fraction gradient, and :math:`\phi` is the phase fraction.


  * ``curvature filter factor``: value of the factor :math:`\beta` applied in the filter :math:`\eta_\kappa = \beta h^2`, where :math:`h` is the cell size. This filter is used to apply a `projection step <https://onlinelibrary.wiley.com/doi/full/10.1002/fld.2643>`_ to damp high frequency errors, that are magnified by differentiation, in the curvature :math:`\kappa`, following the equation:

    .. math:: 
        \int_\Omega \left( v \kappa + \eta_\kappa \nabla v \cdot \nabla \kappa \right) d\Omega = \int_\Omega \left( \nabla v \cdot \frac{\bf{\psi}}{|\bf{\psi}|} \right) d\Omega

    where :math:`v` is a test function, :math:`\kappa` is the filtered curvature, and :math:`\bf{\psi}` is the filtered phase fraction gradient.

  .. tip::

    Use the procedure suggested in: :ref:`choosing values for the surface tension force filters`.

  * ``enable marangoni effect``: Marangoni effect is a thermocapillary effect. It is considered in simulations if this parameter is set to ``true``. Additionally, the ``heat transfer`` auxiliary physics must be enabled (see: :doc:`./multiphysics`) and a non constant ``surface tension model`` with its parameters must be specified in the ``physical properties`` subsection (see: :doc:`./physical_properties`).

.. seealso::

  The surface tension force is used in the :doc:`../../examples/multiphysics/rising-bubble/rising-bubble` example.

.. _choosing values for the surface tension force filters:

Choosing Values for the Surface Tension Force Filters
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

The following procedure is recommended to choose proper values for the ``phase fraction gradient filter factor`` and ``curvature filter factor``:

1. Use ``set output auxiliary fields = true`` to write filtered phase fraction gradient and filtered curvature fields.
2. Choose a value close to 1, for example, :math:`\alpha = 4` and :math:`\beta = 1`.
3. Run the simulation and check whether the filtered phase fraction gradient field is smooth and without oscillation.
4.  If the filtered phase fraction gradient and filtered curvature fields show oscillations, increase the value :math:`\alpha` and :math:`\beta` to larger values, and repeat this process until reaching smooth filtered phase fraction gradient and filtered curvature fields without oscillations.


.. _mass conservation:

Mass Conservation
~~~~~~~~~~~~~~~~~~~~~

* ``subsection mass conservation``: By default, mass conservation (continuity) equations are solved on the whole domain, i.e. on both fluids. This subsection defines parameters that are used to solve mass conservation in both fluids, and to monitor the surface/volume (2D/3D) occupied by the fluid of interest (``monitored fluid``).

  * ``monitoring``: controls if conservation is monitored at each iteration, through the volume (3D) or surface (2D) computation of the fluid given as ``monitored fluid`` (``fluid 1`` (default) or ``fluid 0``). Results are outputted in a data table (`VOF_monitoring_fluid_0.dat` or `VOF_monitoring_fluid_1.dat`).

    .. tip::
      In 2D, the mass returned is in the dimension of :math:`\left[\text{mass}/\text{length}\right]`: multiply this value by the depth of your system to get the ``monitored fluid`` mass (in :math:`\text{kg}` if using SI units).

    .. admonition:: Example of file output, `VOF_monitoring_fluid_1.dat`:

      The ``volume_fluid_1`` column gives the volume occupied by the fluid with index 1, its total mass, and the sharpening threshold used for this iteration.
  
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

  * ``tolerance``: value for the tolerance on the mass conservation of the monitored fluid, used with adaptive sharpening (see the ``subsection sharpening``).
  
    For instance, with ``set tolerance = 0.02`` the sharpening threshold will be adapted so that the mass of the ``monitored fluid`` varies less than :math:`\pm 2\%` from the initial mass (at :math:`t = 0.0` sec).

  * ``verbosity``: states whether from the mass conservation data should be printed. Choices are quiet (no output), verbose (output information from the ``adaptive`` sharpening threshold) and extra verbose (output of the monitoring table in the terminal at the end of the simulation).
