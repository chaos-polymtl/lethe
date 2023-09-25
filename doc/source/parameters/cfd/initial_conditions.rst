==================
Initial Conditions
==================

It is often necessary to set-up complex initial conditions when simulating transient problems. Furthermore, setting a coherent initial condition can greatly help the convergence of steady-state problems. Lethe offers three different strategies to set initial conditions : Nodal (through a function), L2 projection (through a function) and viscous. The viscous strategy solves the same CFD problem but with a user-defined kinematic viscosity. This can allow the user to initialize the velocity field with a Stokes-like flow for which convergence can be obtained more easily.

.. code-block:: text

  subsection initial conditions
    # Type of initial conditions. Choices are L2projection, viscous or nodal
    set type                = nodal

    # Kinematic viscosity for viscous initial conditions
    set kinematic viscosity = 1

    subsection uvwp
      set Function expression = 0; 0; 0; 0
    end

    subsection VOF
      set Function expression = if (x<0.5 & y<1, 1, 0)

      subsection projection step
        set enable           = false
        set diffusion factor = 1
      end
    end

    subsection temperature
      set Function expression = 0
    end

    subsection cahn hilliard
     set Function expression = if (x<0.5 & y<1, 1, -1); 0
    end

    subsection ramp
      subsection kinematic viscosity
        set initial kinematic viscosity = 1.0
        set iterations                  = 0
        set alpha                       = 0.5
      end
      subsection n
        set initial n  = 1.0
        set iterations = 0
        set alpha      = 0.5
      end
    end
  end


* The ``type`` parameter indicates which strategy is used to impose the initial conditions. The choices are : ``L2projection``, ``viscous``, ``nodal`` or ``ramp``.

* The ``kinematic viscosity`` parameter controls the kinematic viscosity that is  used to solve the Navier-Stokes equations when the viscous initial condition is chosen.

* The ``subsection uvwp`` allows the user to select a function (velocity-pressure) to set a nodal or L2 projection initial condition.

* The ``subsection VOF`` defines the areas where both fluids lay at the initial state (see section :doc:`multiphysics`). In this example, the ``Function expression`` is used with a boolean condition to establish the region where the fluid indicator is :math:`0` or :math:`1`: ``if (condition, value if true, value if false)``. ``if (x<0.5 & y<1, 1, 0)`` means that ``fluid 1`` initially fills the surface where ``x<0.5`` and ``y<1``, the rest being filled with ``fluid 0``.

  .. note::
    The ``Function expression`` can be used to establish an even more complex free surface initial geometry. For example, one can create a circle of fluid : ``if ( (x^2+y^2)<=(r)^2 ,1,0)``

  * The ``subsection projection step`` allows to smooth the VOF initial condition using a projection step and avoid a staircase definition of the free surface.

    * When the parameter ``enable`` is set to ``true``, the initial condition is projected following :

    .. math::
      \psi(\Omega_K) = \int_{\Omega_K} \phi d\Omega

    .. math::
      \int_\Omega \left( \psi^* v + \eta_\psi \nabla \psi^* \cdot \nabla v  \right) d\Omega = \int_\Omega \psi v  d\Omega

    where :math:`\psi(\Omega_K)` corresponds to a color function value on the Kth element, :math:`\phi` is the phase fraction, :math:`\psi^*` is the smoothed phase fraction, :math:`\eta_\psi = \alpha h^2` with :math:`\alpha` corresponding to the ``diffusion factor`` and :math:`h` to the cell size, and :math:`v` is a test function.

* The ``subsection cahn hilliard`` defines the areas where both fluids lay at the initial state (see section :doc:`multiphysics`). It works similarly to the ``subsection VOF`` for the first component, which corresponds to the phase order parameter. The user also has the choice to specify initial conditions for the chemical potential, although it is often more suitable to leave it at :math:`0`.
* The ``subsection temperature`` allows the user to define an initial temperature for the fluid domain (if ``set heat tranfer = true`` in :doc:`multiphysics`).

* The ``subsection ramp`` holds the parameters to operate a ramp on either or both the kinematic viscosity and the ``n`` parameter in rheological models (see :doc:`physical_properties` for more information on this parameter). When ramping on the kinematic viscosity value,

  * The ``initial kinematic viscosity`` is the kinematic viscosity with which the initial condition starts off. An initial kinematic viscosity of :math:`1.0` is suggested.
  * The ``iterations`` parameter sets the number of kinematic viscosity iterations before reaching the simulation kinematic viscosity.
  * The ``alpha`` parameter sets the stepping length between kinematic viscosity iterations, as seen in the following equation, where :math:`\eta` is the kinematic viscosity and :math:`i` stands for the iteration number.

.. math::
  \eta_{i+1} = \eta_i + \alpha (\eta_{\text{end}} - \eta_i)

.. note::
  The ramped up kinematic viscosity in the Carreau model in :math:`\eta_0`, and :math:`\eta_{\infty}` stays unchanged. See :doc:`physical_properties` for more details.


Likewise, in the ``subection n``, the parameters for ramping on the ``n`` value are the following.
  * The ``initial n`` is the :math:`n` value with which the initial condition starts off. An initial :math:`n` of :math:`1.0` is suggested.
  * The ``iterations`` parameter sets the number of :math:`n` iterations before reaching the simulation :math:`n`.
  * The ``alpha`` parameter sets the stepping length between :math:`n` iterations, as seen in the following equation, :math:`i` stands for the iteration number.

.. math::
  n_{i+1} = n_i + \alpha (n_{\text{end}} - n_i)
