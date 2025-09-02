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
    
    subsection interface regularization method
      set type       = none
      set frequency  = 10
      set verbosity  = quiet
      
      subsection projection-based interface sharpening
        set interface sharpness     = 2
        set type                    = constant

        # parameter for constant projection-based interface sharpening
        set threshold               = 0.5

        # parameters for adaptive projection-based interface sharpening
        set threshold max deviation = 0.20
        set max iterations          = 20
        set tolerance               = 1e-6
        set monitored fluid         = fluid 1
      end
      
      subsection geometric interface reinitialization
        set max reinitialization distance = 1.0
        set transformation type           = tanh
        set tanh thickness                = 1.0
      end

      subsection algebraic interface reinitialization
        set output reinitialization steps = false
        set steady-state criterion        = 1e-2
        set max steps number              = 5
        set diffusivity multiplier        = 1.0
        set diffusivity power             = 1.0
        set reinitialization CFL          = 1.0
      end
    end

    subsection phase filtration
      set type      = none
      set verbosity = quiet

      # parameter for the tanh filter
      set beta      = 20
    end

    subsection surface tension force
      set enable                                   = false
      set verbosity                                = quiet
      set output auxiliary fields                  = false
      set phase fraction gradient diffusion factor = 4
      set curvature diffusion factor               = 1
      set enable marangoni effect                  = false
    end

  end

* ``viscous dissipative fluid``: defines fluid(s) to which viscous dissipation is applied.

  Choices are: ``fluid 0``, ``fluid 1`` (default) or ``both``, with the fluid IDs defined in Physical properties - :ref:`two phase simulations`.

  .. tip::
    Applying viscous dissipation in one of the fluids instead of both is particularly useful when one of the fluids is air. For numerical stability, the ``kinematic viscosity`` of the air is usually increased. However, we do not want to have viscous dissipation in the air, because it would result in an unrealistic increase in its temperature. This parameter is used only if ``set heat transfer = true`` and ``set viscous dissipation = true`` in :doc:`./multiphysics`.

* ``diffusivity``: value of the diffusivity (diffusion coefficient) in the transport equation of the phase fraction. Default value is ``0`` to have pure advection. 
* ``compressible``: enables interface compression (:math:`\phi \nabla \cdot \mathbf{u}`) in the VOF equation.  This term should be kept to its default value of ``false`` except when compressible equations of state are used.

Interface Regularization Method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``subsection interface regularization method`` defines parameters to counter numerical diffusion of the VOF method and to avoid the interface between the two fluids becoming more and more blurry after each time-step. 

* ``type``: sets the method of regularization. There are four methods available:``none``, ``projection-based interface sharpening``, ``geometric interface reinitialization``, and ``algebraic interface reinitialization``. If ``none`` is selected, the interface is not regularized. The three other types are described bellow along with their corresponding subsection.
* ``frequency``: indicates the frequency at which the regularization process is applied to the VOF phase fraction field. For instance, if the user specifies ``frequency = 2``, the interface will be regularized once every :math:`2` time-steps.

* ``verbosity``: displays the solution process of the regularization method. The different levels of verbosity are:

  * ``quiet``: default verbosity level; no information on the process is displayed.

    .. warning::
      The verbosity of the algebraic interface reinitialization (``type = algebraic``) depends also on the verbosity level of the non-linear and linear solvers. If they are set to ``verbose``, the console outputs of the iteration progress (e.g., norms of the residual and Newton update) may remain.

  * ``verbose``: displays regularization steps progression. For the ``algebraic interface reinitialization``, it only indicates the details of the non-linear and linear iterations if the corresponding solvers are also set to ``verbose``.

  * ``extra verbose``: for the ``projection-based interface sharpening``, indicates the details of the linear iterations. For the ``algebraic interface reinitialization``, in addition to what is displayed at the ``verbose`` level, it displays the steady-state criterion progression through reinitialization steps. This may be used for debugging purposes.
  
Projection-Based Interface Sharpening
+++++++++++++++++++++++++++++++++++++

The ``type = projection-based interface sharpening`` corresponds to a projection-based regularization method in which the phase indicator is projected into a sharper space. The reader is referred to the Projection-Based Interface Sharpening section of :doc:`../../../theory/multiphase/cfd/vof` theory guide for additional details on this regularization method. The ``subsection projection-based interface sharpening`` defines the relevant parameters.

* ``interface sharpness``: sharpness of the moving interface, denoted :math:`\alpha` in the Interface Sharpening section of :doc:`../../../theory/multiphase/cfd/vof` and :math:`a` in the `interface sharpening model <https://www.researchgate.net/publication/287118331_Development_of_efficient_interface_sharpening_procedure_for_viscous_incompressible_flows>`_ paper. This parameter must be larger than 1 for interface sharpening. Choosing values less than 1 leads to interface smoothing instead of sharpening. A good value would be around 1.5.

* ``type``: defines the projection-based interface sharpening type, either ``constant`` or ``adaptive``

  * ``set type = constant``: the sharpening ``threshold`` is the same throughout the simulation. This ``threshold``, between ``0`` and ``1`` (``0.5`` by default), corresponds to the phase fraction at which the interface is located.
  * ``set type = adaptive``: the sharpening threshold is searched in the range :math:`\left[0.5-c_\text{dev} \; ; 0.5+c_\text{dev}\right]`, with :math:`c_\text{dev}` the ``threshold max deviation`` (``0.2`` by default), to ensure mass conservation. The search algorithm will stop either if the mass conservation ``tolerance`` is reached, or if the number of search steps reaches the number of ``max iterations``. If the ``tolerance`` is not reached, a warning message will be printed.

  .. admonition:: Example of a warning message if the sharpening is adaptive but the mass conservation tolerance is not reached:

    .. code-block:: text

      WARNING: Maximum number of iterations (5) reached in the
      adaptive sharpening threshold algorithm, remaining error
      on mass conservation is: 0.02
      Consider increasing the sharpening threshold range or the
      number of iterations to reach the mass conservation tolerance.

  .. tip::

    Usually the first iterations with sharpening are the most at risk to reach the ``max iterations`` without the ``tolerance`` being met, particularly if the mesh is quite coarse.

    As most of the other iterations converge in only one step (corresponding to a final threshold of :math:`0.5`), increasing the sharpening search range through a higher ``threshold max deviation`` will relax the condition on the first iterations with a limited impact on the computational cost.
    
* ``monitored fluid``: Fluid in which the mass conservation is monitored to find the adaptive sharpening threshold. The choices are ``fluid 1`` (default) or ``fluid 0``.

* ``tolerance``: Value of the tolerance on the mass conservation of the monitored fluid.

  For instance, with ``set tolerance = 0.02`` the sharpening threshold will be adapted so that the mass of the ``monitored fluid`` varies less than :math:`\pm 2\%` from the initial mass (at :math:`t = 0.0` sec).

.. seealso::

  The :doc:`../../examples/multiphysics/dam-break/dam-break` example discussed the interface sharperning mechanism.

Geometric Interface Reinitialization
++++++++++++++++++++++++++++++++++++

The ``type = geometric interface reinitialization`` reinitializes the phase fraction field by computing the signed distance from the interface. The latter is then converted back to a phase fraction using a transformation function. The reader is referred to the *Geometric Interface Reinitialization* section of the :doc:`Volume of Fluid method theory guide<../../../theory/multiphase/cfd/vof>` for additional details on this method. The ``geometric interface reinitialization`` sunsection defines the relevant parameters.

* ``max reinitialization distance``: the maximum distance to the interface up to which the signed distance is computed. Above this value, the signed distance is set to the ``max reinitialization distance``.

* ``transformation type``: type of the transformation function used to convert the signed distance to a phase fraction. The choices are: ``tanh`` and ``piecewise polynomial``.
  
  * ``tanh``: the regularized phase fraction is given by :math:`\phi = 0.5-0.5\tanh(d/\varepsilon)`, where :math:`d` is the signed distance to the interface and :math:`\varepsilon` is a measure of the interface thickness and is set by the parameter ``tanh thickness``.
  
  * ``piecewise polynomial``: this transformation uses a piecewise polynomial function of degree 4. It takes the form:
  
    .. math::
      \phi =
      \begin{cases}
        0.5 - 0.5(4d' + 6d'^2 + 4d'^3 + d'^4) \text{ if } d' < 0.0 \\
        0.5 - 0.5(4d' - 6d'^2 + 4d'^3 - d'^4) \text{ if } d' > 0.0
      \end{cases}
    
    where :math:`d' = d/d_\mathrm{max}` is the dimensionless distance to the interface and :math:`d_\mathrm{max}` is the ``max reinitialization distance``.
  
Algebraic Interface Reinitialization
++++++++++++++++++++++++++++++++++++

The ``type = algebraic interface reinitialization`` corresponds to a PDE-based reinitialization method. Alike the projection-based interface sharpening, this aims to reduce numerical diffusion of the phase fraction and redefine the interface sharply by resolving a PDE.  The reader is referred to the *Algebraic Interface Reinitialization* section of the :doc:`Volume of Fluid method theory guide<../../../theory/multiphase/cfd/vof>` for additional details on this method. The ``subsection algebraic interface reinitialization`` defines parameters used to reinitialize the interface in VOF simulations. 

* ``output reinitialization steps``: when set to ``true``, it enables outputs in parallel vtu format of the algebraic reinitialization steps. The files are stored in a folder named ``algebraic-reinitialization-steps-output`` located inside the ``output path`` directory specified in the :doc:`simulation control<./simulation_control>` subsection.

  Outputted quantities of interest are:
    * Reinitialized phase fraction scalar-field (``reinit_phase_fraction``);
    * VOF phase fraction scalar-field (``vof_phase_fraction``);
    * VOF projected phase gradient vector-field (``vof_phase_gradient``) and;
    * VOF projected curvature scalar-field (``vof_curvature``).

  .. tip::
    This feature can be used for debugging purposes by observing how the reinitialization steps affect the phase fraction field.

The interface reinitialization process ends either when steady-state (``steady-state criterion``) is reached or when an imposed maximum number of steps (``max steps number``) is reached.

* ``steady-state criterion``: one of the two stop criteria of the interface reinitialization process. This parameter :math:`(\alpha_\text{ss})` acts as a tolerance for reaching steady-state when solving the algebraic interface reinitialization partial differential equation (PDE).

  .. math::
   \alpha_\text{ss} \geq \frac{ \lVert \phi_\text{reinit}^{\tau + 1} - \phi_\text{reinit}^{\tau} \rVert_2}{\Delta \tau}


  where :math:`\tau` is the pseudo-time used to solve the reinitialization PDE and :math:`\Delta \tau` is the associated pseudo-time-step.

* ``max steps number``: indicates the maximum number of interface reinitialization steps that can be applied before the process ends.

The algebraic interface reinitialization PDE contains a diffusion term. This term contains a diffusion coefficient :math:`(\varepsilon)` given by:

.. math::
  \varepsilon = C h_\text{min}^d

* ``diffusivity multiplier``: factor :math:`(C)` multiplying the smallest cell-size value :math:`(h_\text{min})` in the evaluation of the diffusion coefficient of the PDE.

* ``diffusivity power``: power :math:`(d)` to which the smallest cell-size value :math:`(h_\text{min})` is elevated in the evaluation of the diffusion coefficient of the PDE.

* ``reinitialization CFL``: CFL condition of the interface reinitialization process. This is used to evaluate the pseudo-time-step :math:`(\Delta\tau)`.

  .. math::
    \Delta \tau = C_\text{CFL} \, h_\text{min}

  where :math:`C_\text{CFL}` is the imposed CFL condition and :math:`h_\text{min}` is the size of the smallest cell.

Phase Filtration
~~~~~~~~~~~~~~~~~~

* ``subsection phase filtration``: defines the filter applied to the phase fraction. This affects the definition of the interface.

* ``type``: defines the filter type, either ``none`` or ``tanh``

  * ``set type = none``: the phase fraction is not filtered
  * ``set type = tanh``: the filter function described in the Interface filtration section of :doc:`../../../theory/multiphase/cfd/vof` theory guide is applied.
* ``beta``: value of the :math:`\beta` parameter of the ``tanh`` filter
* ``verbosity``: enables the display of filtered phase fraction values. Choices are ``quiet`` (no output) and ``verbose`` (displays values)


Surface Tension Force
~~~~~~~~~~~~~~~~~~~~~~

* ``subsection surface tension force``: Surface tension is the tendency of a liquid to maintain the minimum possible surface area. This subsection defines parameters to ensure an accurate interface between the two phases, used when at least one phase is liquid. 

  * ``enable``: controls if ``surface tension force`` is considered.

    .. attention::

      When the surface tension force is enabled, a ``fluid-fluid`` material interaction and a ``surface tension model`` with its parameters must be specified in the :doc:`physical_properties` subsection.

  * ``verbosity``: enables the display of the output from the surface tension force calculations. Choices are: ``quiet`` (default, no output) and ``verbose``.
  * ``output auxiliary fields``: enables the display of the projected ``phase fraction gradient`` and projected ``curvature``. Used for debugging purposes.

  * ``phase fraction gradient diffusion factor``: value of the factor :math:`\alpha` in :math:`\eta_n = \alpha h^2`, where :math:`h` is the cell size. This diffusion coefficient (:math:`\eta_n`) is used in a `projection step <https://onlinelibrary.wiley.com/doi/full/10.1002/fld.2643>`_ to damp high frequency errors, that are magnified by differentiation, in the phase fraction gradient (:math:`\bf{\psi}`), following the equation:

    .. math::
        \int_\Omega \left( {\bf{v}} \cdot {\bf{\psi}} + \eta_n \nabla {\bf{v}} \cdot \nabla {\bf{\psi}} \right) d\Omega = \int_\Omega \left( {\bf{v}} \cdot \nabla {\phi} \right) d\Omega

    where :math:`\bf{v}` is a piecewise continuous vector-valued test function, :math:`\bf{\psi}` is the projected phase fraction gradient, and :math:`\phi` is the phase fraction.


  * ``curvature diffusion factor``: value of the factor :math:`\beta` in  :math:`\eta_\kappa = \beta h^2`, where :math:`h` is the cell size. This diffusion coefficient (:math:`\eta_\kappa`) is used in a `projection step <https://onlinelibrary.wiley.com/doi/full/10.1002/fld.2643>`_ to damp high frequency errors, that are magnified by differentiation, in the curvature :math:`\kappa`, following the equation:

    .. math:: 
        \int_\Omega \left( v \kappa + \eta_\kappa \nabla v \cdot \nabla \kappa \right) d\Omega = \int_\Omega \left( \nabla v \cdot \frac{\bf{\psi}}{|\bf{\psi}|} \right) d\Omega

    where :math:`v` is a test function, :math:`\kappa` is the projected curvature, and :math:`\bf{\psi}` is the projected phase fraction gradient.

  .. tip::

    Use the procedure suggested in: :ref:`choosing values for the surface tension force diffusion factors`.

  * ``enable marangoni effect``: Marangoni effect is a thermocapillary effect. It is considered in simulations if this parameter is set to ``true``. Additionally, the ``heat transfer`` auxiliary physics must be enabled (see: :doc:`./multiphysics`) and a non constant ``surface tension model`` with its parameters must be specified in the ``physical properties`` subsection (see: :doc:`./physical_properties`).

.. seealso::

  The surface tension force is used in the :doc:`../../examples/multiphysics/rising-bubble/rising-bubble` example.

.. _choosing values for the surface tension force diffusion factors:

Choosing Values for the Surface Tension Force Diffusion Factors
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The following procedure is recommended to choose proper values for the ``phase fraction gradient diffusion factor`` and ``curvature diffusion factor``:

1. Use ``set output auxiliary fields = true`` to write projected phase fraction gradient and projected curvature fields.
2. Choose a value close to 1, for example, :math:`\alpha = 4` and :math:`\beta = 1`.
3. Run the simulation and check whether the projected phase fraction gradient field is smooth and without oscillation.
4.  If the projected phase fraction gradient and projected curvature fields show oscillations, increase the value :math:`\alpha` and :math:`\beta`, and repeat this process until reaching smooth projected phase fraction gradient and projected curvature fields without oscillations.
