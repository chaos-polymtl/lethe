===================
Physical Properties
===================

.. note:: 
    Lethe supports single phase, two phase (using VOF) and conjugate simulations. This is managed using the fluid and solid subsections.

.. code-block:: text

  subsection physical properties
    set number of fluids = 1
    subsection fluid 0
      # Rheology
      set rheological model          = newtonian
      set kinematic viscosity        = 1
      
      # Density
      set density model              = constant
      set density                    = 1
      
      # Specific heat
      set specific heat model        = constant
      set specific heat              = 1
      
      # Thermal conductivity
      set thermal conductivity model = constant
      set thermal conductivity       = 1
      
      # Thermal expansion
      set thermal expansion model    = constant
      set thermal expansion          = 0
      
      # Tracer diffusivity
      set tracer diffusivity model   = constant
      set tracer diffusivity         = 0
    end

    set number of solids = 0

    set number of material interactions = 1
    subsection material interaction 0
      set type = fluid-fluid
      subsection fluid-fluid interaction
        set first fluid id              = 0
        set second fluid id             = 1

        # Surface tension
        set surface tension model       = constant
        set surface tension coefficient = 0
      end

      # if fluid-solid interaction
      subsection fluid-solid interaction
        set fluid id                    = 0
        set solid id                    = 0

        # Surface tension
        set surface tension model       = constant
        set surface tension coefficient = 0
      end
    end
  end
 
* The ``number of fluids`` parameter controls the number of fluids in the simulation. This parameter is set to ``1`` except in `Two Phase Simulations`_ .

  * The ``rheological model`` parameter sets the choice of rheological model. The choices are between ``newtonian``, ``power-law``, ``carreau`` and ``phase_change``. For more details on the rheological models, see  `Rheological Models`_ .

  * The ``kinematic viscosity`` parameter is the kinematic viscosity of the newtonian fluid in units of :math:`\text{Length}^{2} \cdot \text{Time}^{-1}`. In SI, this is :math:`\text{m}^{2} \cdot \text{s}^{-1}`. This viscosity is only used when ``rheological model = newtonian``.

  * The ``density model`` specifies the model used to calculate the density. At the moment, a ``constant`` density and an ``isothermal_ideal_gas`` model are supported. For more details on the density models, see `Density Models`_.

  * The ``density`` parameter is the constant density of the fluid in units of :math:`\text{Mass} \cdot \text{Length}^{-3}`

  * The ``specific heat model`` specifies the model used to calculate the specific heat. At the moment, only a constant specific heat is supported.

  * The ``specific heat`` parameter is the constant specific heat of the fluid in units of :math:`\text{Energy} \cdot \text{Temperature}^{-1} \cdot \text{Mass}^{-1}` .

  * The ``thermal conductivity model`` specifies the model used to calculate the thermal conductivity. At the moment, ``constant`` and ``linear`` thermal conductivity are available. For more details on the thermal conductivity models, see `Thermal Conductivity Models`_.

  * The ``thermal conductivity`` parameter is the thermal conductivity coefficient of the fluid with units of :math:`\text{Power} \cdot \text{Temperature}^{-1} \cdot \text{Length}^{-1}`.

  * The ``thermal expansion model`` specifies the model used to calculate the thermal expansion coefficient. At the moment, only a constant thermal expansion is supported.

  * The ``thermal expansion`` parameter is the thermal expansion coefficient of the fluid with dimension of :math:`\text{Temperature}^{-1}`. It is used to define the buoyancy-driven flow (natural convection) using the Boussinesq approximation, which leads to the definition of the following source term that is added to the Navier-Stokes equation:

    .. math::

      {\bf{F_{B}}} = -\beta {\bf{g}} (T-T_0)

    where :math:`F_B` denotes the buoyant force source term, :math:`\beta` is the thermal expansion coefficient, :math:`T` is temperature, and :math:`T_0` is the base temperature.

  * The ``tracer diffusivity model`` specifies the model used to calculate the tracer diffusivity. At the moment, only a constant tracer diffusivity is supported.

  * The ``tracer diffusivity`` parameter is the diffusivity coefficient of the tracer in units of :math:`\text{Length}^{2} \cdot \text{Time}^{-1}` . In SI, this is :math:`\text{m}^{2} \cdot \text{s}^{-1}`.

* The ``number of solids`` parameter controls the number of solid regions. Solid regions are currently only implemented for `Conjugate Heat Transfer`_.

* The ``number of material interactions`` parameter controls the number of physical properties that are due to the interaction between two materials. At the moment, only the surface tension between two fluids is implemented in `Two Phase Simulations`_.

  * The material interaction ``type`` can either be ``fluid-fluid`` (default) or ``fluid-solid``.

  * In the ``fluid-fluid`` subsection we define the pair of fluids and their physical properties.

    * The ``first fluid id`` is the id of the first fluid.

    * The ``second fluid id`` is the id of the second fluid.

      .. attention::
          The ``second fluid id`` should be greater than the ``first fluid id``.

    * The ``surface tension model`` specifies the model used to calculate the surface tension coefficient of the fluid-fluid pair. At the moment, only the ``constant`` model is supported.

    * The ``surface tension coefficient`` parameter is the constant surface tension coefficient of the two interacting fluids in units of :math:`\text{Mass} \cdot \text{Time}^{-2}`. In SI, this is :math:`\text{N} \cdot \text{m}^{-1}`. The surface tension coefficient is used to define the Weber number (:math:`We`):

      .. math::
          We = Re \cdot \frac{\mu_\text{ref} \; u_\text{ref}}{\sigma}

      where :math:`Re` is the Reynolds number, :math:`\mu_\text{ref}` and :math:`u_\text{ref}` are some reference viscosity and velocity characterizing the flow problem, and :math:`\sigma` is the surface tension coefficient.

  * In the ``fluid-solid`` subsection we define the fluid-solid pair and their physical properties.

    * The ``fluid id`` is the id of the fluid.

    * The ``solid id`` is the id of the solid.

    * The ``surface tension model``  and ``surface tension coefficient`` are the same as described in the ``fluid-fluid`` subsection above.

.. note:: 
  The default values for all physical properties models in Lethe is ``constant``. Consequently, it is not necessary (and not recommended) to specify the physical property model unless this model is not constant. This generates parameter files that are easier to read.

.. _two phase simulations:

Two Phase Simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. note:: 
  Two phase simulations require that either ``set VOF = true`` or ``set cahn hilliard = true`` in the :doc:`multiphysics` subsection. By convention, air is usually the ``fluid 0`` and the other fluid of interest is the ``fluid 1``.

For two phases, the properties are defined for each fluid. Default values are:

.. code-block:: text

  subsection physical properties
  set number of fluids = 2
      subsection fluid 0
         set density              = 1
         set kinematic viscosity  = 1
         set specific heat        = 1
         set thermal conductivity = 1
         set tracer diffusivity   = 0
      end
      subsection fluid 1
         set density              = 1
         set kinematic viscosity  = 1
         set specific heat        = 1
         set thermal conductivity = 1
         set tracer diffusivity   = 0
      end
  end

* ``number of fluids = 2`` is required for a free surface simulation, otherwise an error will be thrown in the terminal.
* ``subsection fluid 0`` indicates the properties of fluid where the phase indicator = 0 (Volume of Fluid method), as defined when initializing the free surface (see the :doc:`initial_conditions` subsection), and correspondingly ``fluid 1`` is located where the phase indicator = 1.

.. warning:: 
  Lethe now supports the use of physical properties models that are different for both phases. For example, the liquid could have a carreau rheological model and the air could have a newtonian rheological model. However, this feature has not been fully tested and could lead to unpredictable results. Use with caution.


.. _conjugate heat transfer:

Conjugate Heat Transfer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Conjugate heat transfer enables the addition of solid regions in which the fluid dynamics is not solved for. To enable the presence of a solid region, ``number of solids`` must be set to 1. By default, the region with the ``material_id=0`` will be the fluid region whereas the region with ``material_id=1`` will be the solid region. The physical properties of the solid region are set in an identical fashion as those of the fluid.

.. warning::
  This is an experimental feature. It has not been tested on a large range of application cases. 

.. code-block:: text

  subsection physical properties
    set number of fluids = 1
    subsection fluid 0
      ...
    end
    set number of solids = 1
    subsection solid 0
      # Density
      set density model              = constant
      set density                    = 1
      
      # Specific heat
      set specific heat model        = constant
      set specific heat              = 1
      
      # Thermal conductivity
      set thermal conductivity model = constant
      set thermal conductivity       = 1
    end
  end

.. _rheological_models:

Rheological Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two families of rheological models are supported in Lethe. The first one are generalized non Newtonian rheologies (for shear thinning and shear thickening flows). In these models, the viscosity depends on the shear rate. The second family of rheological models possess a viscosity that is independent of the shear rate, but that may depend on other fields such as the temperature.

The ``rheological model`` parameter sets which rheological model you are using. The default ``rheological model`` is ``newtonian``, which uses a constant ``kinematic viscosity``.

.. code-block:: text

    subsection physical properties
      set number of fluids = 1
      subsection fluid 0
        set rheological model   = newtonian
        set kinematic viscosity = 1.0
      end
    end

The rheological model available options are:
    * ``newtonian``
    * ``power-law`` 
    * ``carreau``
    * ``phase_change``

Power-Law Model
^^^^^^^^^^^^^^^

The power-law model is the simplest rheological model, using only 2 parameters 

.. math::

  \eta(\dot{\gamma}) = K \dot{\gamma}^{n-1}


where :math:`\eta` is the **kinematic viscosity** and :math:`\dot{\gamma}` is the local shear rate magnitude.

.. image:: images/physical_properties_powerlaw.png
    :width: 600
    :align: center

When using the power-law model, the default values are:

.. code-block:: text

  subsection physical properties
    set number of fluids = 1
    subsection fluid 0
      set rheological model   = power-law
      subsection non newtonian
        subsection power-law
          set K               = 1.0
          set n               = 0.5
          set shear rate min  = 1e-3
        end
      end
    end
  end

* The ``K`` parameter is a fluid consistency index. It represents the fluid viscosity for a local shear rate of :math:`1.0`.

* The ``n`` parameter is the flow behavior index. It sets the slope in the log-log :math:`\eta = f(\dot{\gamma})` graph.

* The ``shear rate min`` parameter yields the magnitude of the shear rate tensor for which the viscosity is calculated. Since the model uses a power operation, a null shear rate magnitude leads to an invalid viscosity. To ensure numerical stability, the shear rate cannot go below this threshold when the viscosity  calculated.

Carreau Model
^^^^^^^^^^^^^^^

The Carreau model is in reality the five parameter Carreau model:

.. math::

  \eta(\dot{\gamma}) =\eta_{\infty} + (\eta_0 - \eta_{\infty}) \left[ 1 + (\dot{\gamma}\lambda)^a\right]^{\frac{n-1}{a}}
 
where :math:`\eta` is the **kinematic viscosity** and :math:`\dot{\gamma}` is the shear rate.

.. image:: images/physical_properties_carreau.png
    :width: 600
    :align: center

The parameters for the Carreau model are defined by the ``carreau`` subsection. The default values are:

.. code-block:: text

  subsection physical properties
    set number of fluids = 1
    subsection fluid 0
      set rheological model   = carreau
      subsection non newtonian
        subsection carreau
          set viscosity_0     = 1.0
          set viscosity_inf   = 1.0
          set a               = 2.0
          set lambda          = 1.0
          set n               = 0.5
        end
      end
    end
  end

* The ``viscosity_0`` parameter represents the viscosity when the shear rate on the fluid tends to 0.

* The ``viscosity_inf`` parameter represents the viscosity when the shear rate on the fluid becomes large.

* The ``a`` is the Carreau parameter, generally set to 2.

* The ``lambda`` is the relaxation time associated to the fluid.

* The ``n`` is a power parameter. It sets the slope in the log-log :math:`\eta = f(\dot{\gamma})` graph just like in the power-law model.

.. note::
    The Carreau model is only suitable for Newtonian and shear-thinning flows.

Phase-Change Model
^^^^^^^^^^^^^^^^^^^ 

The phase change model is a simple rheological model in which the viscosity depends on the temperature. This model is used to model melting and freezing of components. The kinematic viscosity :math:`\nu` is given by :

.. math::

  \nu =   c^{*}_p  = \begin{cases} \nu_s \; \text{if} \; T<T_{s} \\
              \frac{T-T_s}{T_l-T_s} \nu_l + (1-\frac{T-T_s}{T_l-T_s}) \nu_s \; \text{if} \; T_{l}>T>T_{s}\\
              \nu_l \; \text{if} \; T>T_{l}
              \end{cases}

where :math:`T_l` and :math:`T_s` are the liquidus and solidus temperature. The underlying hypothesis of this model is that the melting and the solidification occurs over a phase change interval. Melting will occur between :math:`T_s` and :math:`T_l` and solidification will occur between :math:`T_l` and :math:`T_s`.

This model is parameterized using the ``phase change`` subsection

.. code-block:: text

  subsection phase change
    # Temperature of the liquidus
    set liquidus temperature = 1
  
    # Temperature of the solidus
    set solidus temperature  = 0

    # Specific heat of the liquid phase
    set viscosity liquid     = 1
  
    # viscosity of the solid phase
    set viscosity solid      = 1
  end


* The ``liquidus temperature`` is :math:`T_l`

* The ``solidus temperature`` is :math:`T_s`

* The ``specific heat liquid`` is :math:`\nu_{l}`

* The ``specific heat solid`` is :math:`\nu_{s}`

.. note::
  The phase change subsection is used to parametrize *both* ``rheological model = phase_change`` *and* ``specific heat model = phase_change``. This prevents parameter duplication.

.. _density_models:

Density Models
~~~~~~~~~~~~~~

Lethe supports both ``constant`` and ``isothermal_ideal_gas`` density models. Constant density assumes a constant density value. Isothermal ideal gas density assumes that the fluid's density varies according the following state equation:

.. math::
  \rho = \rho_{ref} + \psi p = \rho_{ref} + \frac{1}{R T} \ p

where :math:`\rho_{ref}` is the density of the fluid at the reference state, :math:`\psi = \frac{1}{R T}` is the compressibility factor derived from the ideal gas law with :math:`R= \frac{R_u}{M}` the specific gas constant (universal gas constant (:math:`R_u`) divided by the molar mass of the gas (:math:`M`)) and :math:`T` the temperature of the gas, finally, :math:`p` is the differential pressure between the reference state and the current state.

This model is parametrized using the ``isothermal_ideal_gas`` subsection:

.. code-block:: text

  subsection physical properties
    set number of fluids = 1
    subsection fluid 0
      set density model = isothermal_ideal_gas
      subsection isothermal_ideal_gas
        set density_ref = 1.2
        set R           = 287.05
        set T           = 293.15
      end
    end
  end

where:

* ``density_ref`` corresponds to :math:`\rho_{ref}`

* ``R`` corresponds to :math:`R`

* ``T`` corresponds to :math:`T`

By default, parameters are set to the values of dry air evaluated under normal temperature and pressure conditions :math:`(20 \ \text{°C}`, :math:`1 \ \text{atm})`.

.. caution::
  When defining the initial pressure condition in the ``initial conditions`` subsection (see :doc:`initial_conditions`), make sure to set it to :math:`0`, as it represents the reference state for the calculated pressure. In solving the Navier-Stokes equations, the pressure is defined to within a constant. Therefore, it is more appropriate to interpret it as a differential pressure.

  .. attention::
    Currently, in two phase flow simulations, the ``isothermal_ideal_gas`` density model can only be used in conjunction with the incompressible Navier-Stokes equations. However, it is meant to be used with the isothermal formulation of compressible Navier-Stokes equations to account for weakly compressible flows. In a future update, this change will be implemented.

.. _thermal_conductivity_models:

Thermal Conductivity Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Constant, linear and phase_change thermal conductivities are supported in Lethe. Constant thermal conductivity assumes a constant value of the thermal conductivity. Linear thermal conductivity assumes that that the thermal conductivity :math:`k` varies linearly with the temperature, taking the following form:

.. math::
  k = k_{A,0}+ k_{A,1} T 

where :math:`k_{A,0}` and :math:`k_{A,1}` are constants and :math:`T` is the temperature. This enables a linear variation of the thermal conductivity as a function of the temperature.

In the ``phase_change`` thermal conductivity model, two different values (``thermal conductivity liquid``, and ``thermal conductivity solid``) are required for calculating the thermal conductivities of the liquid and solid phases, respectively. For the liquid phase (T>T_liquidus), the ``thermal conductivity liquid`` is applied, while for the solid phase (T<T_solidus), the model uses the ``thermal conductivity solid``. In the mushy zone between T_solidus and T_liquidus, the thermal conductivity is equal to:

.. math::

  k = \alpha_l k_l + (1 - \alpha_l) k_s


where :math:`k_l`, :math:`k_s` and  :math:`\alpha_l` denote thermal conductivities of the liquid and solid phases and the liquid fraction.

Specific Heat Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lethe supports two types of specific heat models. Setting ``specific heat=constant`` sets a constant specific heat. Lethe also supports a ``phase_change`` specific heat model. This model can simulate the melting and solidification of a material. The model follows the work of Blais & Ilinca `[1] <https://doi.org/10.1016/j.compfluid.2018.03.037>`_. This approach defines the specific heat :math:`C_p` as:

.. math::

  C_p = \frac{H(T)-H(T_0)}{T-T_0}


where :math:`T` is the temperature, :math:`T_0` is the temperature at the previous time and :math:`H(T)` is the enthalpy, as a function of the temperature, to be:

.. math::
  H(T) = H_0 + \int_{T_0}^{T} c^{*}_p (T^*) dT


where :math:`H_0` is a reference enthalpy, taken to be 0, and :math:`c^{*}_p` is:

.. math::
  c^{*}_p  = \begin{cases} C_{p,s}\\
              \frac{C_{p,s}+C_{p,l}}{2}+\frac{h_l}{T_l-T_s}\\
              C_{p,l}
              \end{cases}

where :math:`C_{p,s}` and :math:`C_{p,l}` are the solid and liquid specific heat, respectively. :math:`h_l` is the latent enthalpy (enthalpy related to the phase change), :math:`T_l` and :math:`T_s` are the liquidus and solidus temperature. The underlying hypothesis of this model is that the melting and the solidification occurs over a phase change interval. Melting will occur between :math:`T_s` and :math:`T_l` and solidification will occur between :math:`T_l` and :math:`T_s`.

This model is parameterized using the following section:

.. code-block:: text

  subsection phase change
    # Enthalpy of the phase change
    set latent enthalpy      = 1
  
    # Temperature of the liquidus
    set liquidus temperature = 1
  
    # Temperature of the solidus
    set solidus temperature  = 0
  
    # Specific heat of the liquid phase
    set specific heat liquid = 1
  
    # Specific heat of the solid phase
    set specific heat solid  = 1
  end

* The ``latent enthalpy`` is the latent enthalpy of the phase change: :math:`h_l`

* The ``liquidus temperature`` is :math:`T_l`

* The ``solidus temperature`` is :math:`T_s`

* The ``specific heat liquid`` is :math:`C_{p,l}`

* The ``specific heat solid`` is :math:`C_{p,s}`


`[1] <https://doi.org/10.1016/j.compfluid.2018.03.037>`_ B. Blais and F. Ilinca, “Development and validation of a stabilized immersed boundary CFD model for freezing and melting with natural convection,” *Comput. Fluids*, vol. 172, pp. 564–581, Aug. 2018, doi: 10.1016/j.compfluid.2018.03.037.
