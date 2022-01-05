Physical Properties
---------------------
.. note:: 
    Lethe supports both single phase and two phase (through VOF) simulations. The same subsection is used to manage both types of simulation using the fluid subsections.

.. code-block:: text

  subsection physical properties
  set number of fluids	= 1
     subsection fluid 0
       set density 		= 1
       set kinematic viscosity 	= 1
       set specific heat 	= 1
       set thermal conductivity = 1
       set tracer diffusivity   = 0
     end
  end

* The ``kinematic viscosity`` parameter is the kinematic viscosity of the fluid in units of :math:`Length^{2} \cdot Time^{-1}`. In SI this is :math:`m^{2} \cdot s^{-1}`.

* The ``density`` parameter is the density of the fluid in units of :math:`Mass \cdot Length^{-3}`

* The ``specific heat`` parameter is the specific heat of the fluid in units of :math:`Energy \cdot Temperature^{-1} \cdot Mass^{-1}` .

* The ``thermal expansion`` parameter is the thermal expansion coefficient of the fluid with dimension of :math:`Temperature^{-1}`.

* The ``tracer diffusivity`` parameter is the diffusivity coefficient of the tracer in units of :math:`Length^{2} \cdot Time^{-1}` . In SI this is :math:`m^{2} \cdot s^{-1}` .

Two phase simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. note:: 
  Two phase simulations require that ``set VOF = true`` in the Multiphysics subsection. By convention, air is usually the ``fluid 0`` and the other fluid of interest is the ``fluid 1``.

For two phases, the properties are defined for each fluid. Default values are:

.. code-block:: text

  subsection physical properties
  set number of fluids		= 2
      subsection fluid 0
         set density 		= 1
         set kinematic viscosity 	= 1
         set specific heat 	= 1
         set thermal conductivity = 1
         set tracer diffusivity   = 0
      end
      subsection fluid 1
         set density 		= 1
         set kinematic viscosity 	= 1
         set specific heat 	= 1
         set thermal conductivity = 1
         set tracer diffusivity   = 0
      end
  end

* ``number of fluids = 2`` is required for a free surface simulation, otherwise an error will be thrown in the terminal.
* ``subsection fluid 0`` indicates the properties of fluid where the phase indicator = 0 (Volume of Fluid method), as defined when initializing the free surface (see the initial conditions subsection), and correspondingly ``fluid 1`` is located where the phase indicator = 1.
