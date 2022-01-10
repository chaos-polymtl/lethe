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

* The ``kinematic viscosity`` parameter is the kinematic viscosity of the fluid in units of :math:`\text{Length}^{2} \cdot \text{Time}^{-1}`. In SI this is :math:`\text{m}^{2} \cdot \text{s}^{-1}`.

* The ``density`` parameter is the density of the fluid in units of :math:`\text{Mass} \cdot \text{Length}^{-3}`

* The ``specific heat`` parameter is the specific heat of the fluid in units of :math:`\text{Energy} \cdot \text{Temperature}^{-1} \cdot \text{Mass}^{-1}` .

* The ``thermal expansion`` parameter is the thermal expansion coefficient of the fluid with dimension of :math:`\text{Temperature}^{-1}`.

* The ``tracer diffusivity`` parameter is the diffusivity coefficient of the tracer in units of :math:`\text{Length}^{2} \cdot \text{Time}^{-1}` . In SI this is :math:`\text{m}^{2} \cdot \text{s}^{-1}` .

Two phase simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. note:: 
  Two phase simulations require that ``set VOF = true`` in the :doc:`multiphysics` subsection. By convention, air is usually the ``fluid 0`` and the other fluid of interest is the ``fluid 1``.

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
* ``subsection fluid 0`` indicates the properties of fluid where the phase indicator = 0 (Volume of Fluid method), as defined when initializing the free surface (see the :doc:`initial_conditions` subsection), and correspondingly ``fluid 1`` is located where the phase indicator = 1.

Rheological models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generalized non Newtonian rheologies (for shear thinning and shear thickening flows) are supported in Lethe. 

.. note:: 
  Currently, non Newtonian flow simulations are only supported using one single fluid. A possibility for multiple Newtonian/non Newtonian fluids is upcoming!
  
Default values for a non Newtonian fluid are

.. code-block:: text

    subsection physical properties
    set non newtonian flow	= false
        subsection non newtonian
        set model 		= carreau
        end
    end
    
* The ``non newtonian flow`` parameter has to be set to ``true`` to use a rheological model.

* The ``model`` parameter sets which rheological model you are using. The available options are:
    * ``carreau``
    * ``power-law`` 

The Carreau model is in reality the five parameter Carreau model :

.. math::

  \eta(\dot{\gamma}) =\eta_{\infty} + (\eta_0 - \eta_{\infty}) \left[ 1 + (\dot{\gamma}\lambda)^a\right]^{\frac{n-1}{a}}
 
where :math:`\eta` is the **kinematic viscosity** and :math:`\dot{\gamma}` is the shear rate.

The parameters for the Carreau model are defined by the ``carreau`` subsection. The default values are:

.. code-block:: text

  subsection physical properties
    set non newtonian flow	= true
      subsection non newtonian
        set model 		= carreau
        subsection carreau
          set viscosity_0	= 1.0
          set viscosity_inf = 1.0
          set a = 2.0
          set lambda = 1.0
          set n = 0.5
        end
    end
  end

* The ``viscosity_0`` parameter represents the viscosity when the shear rate on the fluid tends to 0.

* The ``viscosity_inf`` parameter represents the viscosity when the shear rate on the fluid becomes large.

* The ``a`` is the Carreau parameter, generally set to 2.

* The ``lambda`` is the relaxation time associated to the fluid.

* The ``n`` is a power parameter. It sets the slope in the log-log :math:`\eta = f(\dot{\gamma})` graph.


The power-law model is a simple rheological model:

.. math::

  \eta(\dot{\gamma}) = K \dot{\gamma}^{n-1}


where :math:`\eta` is the **kinematic viscosity** and :math:`\dot{\gamma}` is the shear rate.
When using the Power-Law model, the default values are:

.. code-block:: text

    subsection physical properties
    set non newtonian flow	= true
        subsection non newtonian
        set model 		= power-law
            subsection power-law
              set K = 1.0
              set n = 0.5
              set shear rate min = 1e-3
            end
        end
    end

* The ``K`` parameter is a fluid consistency index. It represents the fluid viscosity is it were Newtonian.

* The ``n`` parameter is the flow behavior index. low  It sets the slope in the log-log :math:`\eta = f(\dot{\gamma})` graph.

* The ``shear rate min`` parameter yields the magnitude of the shear rate tensor for which the viscosity is calculated. Since the model uses a power operation, a nul shear rate magnitude leads to an invalid viscosity. To ensure numerical stability, the shear rate cannot go below this threshold when the viscosity  calculated.


