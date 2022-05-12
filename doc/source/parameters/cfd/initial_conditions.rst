Initial Conditions
-------------------
It is often necessary to set-up complex initial conditions when simulating transient problems. Furthermore, setting a coherent initial condition can greatly help the convergence of steady-state problems. Lethe offers three different strategies to set initial conditions : Nodal (through a function), L2 projection (through a function) and viscous. The viscous strategy solves the same CFD problem but with a user-defined viscosity. This can allow the user to initialize the velocity field with a Stokes-like flow for which convergence can be obtained more easily.

.. code-block:: text

 subsection initial conditions
   # Type of initial conditions. Choices are L2projection, viscous or nodal
   set type      = nodal

   # viscosity for viscous initial conditions
   set viscosity = 1

   subsection uvwp
     set Function expression = 0; 0; 0; 0 
   end

   subsection VOF
     set Function expression = if (x<0.5 & y<1, 1, 0)
   end

   subsection ramp
    subsection viscosity
      set initial viscosity = 1.0
      set iterations = 0
      set alpha = 0.5
    end
    subsection n
      set initial n = 1.0
      set iterations = 0
      set alpha = 0.5
    end
  end


* The ``type`` parameter indicates which strategy is used to impose the initial conditions. The choices are : ``L2projection``, ``viscous``, ``nodal`` or ``ramp``.

* The ``viscosity`` parameter controls the viscosity that is  used to solve the Navier-Stokes equations when the viscous initial condition is chosen.

* The ``subsection uvwp`` allows the user to select a function (velocity-pressure) to set a nodal or L2 projection initial condition.

* The ``subsection VOF`` defines the areas where both fluids lay at the initial state (see section :doc:`multiphysics`). In this example, the ``Function expression`` is used with a boolean condition to establish the region where the fluid indicator is 0 or 1: ``if (condition, value if true, value if false)``. ``if (x<0.5 & y<1, 1, 0)`` means that ``fluid 1`` initially fills the surface where ``x<0.5`` and ``y<1``, the rest being filled with ``fluid 0``.

.. note::
   The ``Function expression`` can be used to establish an even more complex free surface initial geometry. For example, one can create a circle of fluid : ``if ( (x^2+y^2)<=(r)^2 ,1,0)``


* The ``subsection ramp`` holds the parameters to operate a ramp on either or both the viscosity and the ``n`` parameter in rheological models (see :doc:`physical_properties` for more information on this parameter). When ramping on the viscosity value, 

  * The ``initial viscosity`` is the viscosity with which the initial condition starts off. An initial viscosity of 1.0 is suggested.
  * The ``iterations`` parameter sets the number of viscosity iterations before reaching the simulation viscosity. 
  * The ``alpha`` parameter sets the stepping length between viscosity iterations, as seen in the following equation, where :math:`\eta` is the viscosity and :math:`i` stands for the iteration number.

.. math::
  \eta_{i+1} = \eta_i + \alpha (\eta_{\text{end}} - \eta_i)

.. note::
  The ramped up viscosity in the Carreau model in :math:`\eta_0`, and :math:`\eta_{\infty}` stays unchanged. See :doc:`physical_properties` for more details.


Likewise, in the ``subection n``, the parameters for ramping on the ``n`` value are the following.
  * The ``initial n`` is the :math:`n` value with which the initial condition starts off. An initial math:`n` of 1.0 is suggested.
  * The ``iterations`` parameter sets the number of math:`n` iterations before reaching the simulation math:`n`. 
  * The ``alpha`` parameter sets the stepping length between math:`n` iterations, as seen in the following equation, :math:`i` stands for the iteration number.

.. math::
  n_{i+1} = n_i + \alpha (n_{\text{end}} - n_i)
