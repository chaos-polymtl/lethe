============================
Quadrature Method of Moments
============================

This subsection defines the parameters of the quadrature method of moments (QMOM) auxiliary physics, used to solve population balance equations. The physics is enabled with ``set qmom = true`` in the :doc:`multiphysics` subsection and requires Lethe to be compiled with `QMOMKit <https://github.com/polycfd/QMOMKit>`_ (``-DLETHE_WITH_QMOMKIT=ON``).

.. note::
    For more information on the QMOMKit library, including details on the implementation and some application examples, visit the `repository <https://github.com/polycfd/QMOMKit>`_ page.

The QMOM physics transports ``number of scalars`` user-defined scalar fields and ``number of moments`` :math:`K` moments :math:`M_0, \dots, M_{K-1}` of the particle size distribution, discretized with a discontinuous Galerkin finite element of degree ``qmom degree`` (see :doc:`fem`). 

With the ``lie`` splitting scheme, each time step is split into the transport of the scalars and of the moments in physical space, independently from each other, followed by the evolution of the moments in phase space (growth, nucleation and coalescence) by QMOMKit. During the phase-space evolution, the scalars keep their values at the end of their transport; if the process couples them to the particles, they are then updated with the change of the moments. The order of these steps is detailed in :ref:`qmom-operator-splitting`. The scalars are named after the ``scalar names`` and the moments are named ``moment_0`` to ``moment_<K-1>``. The scalars may diffuse, with the ``scalar diffusivities``, while the moments are only advected.

.. seealso::
    The governing equations, the discretization and the phase-space models are described in :doc:`../../theory/multiphysics/qmom/qmom_equations`.

.. caution::
    All the parameters of the QMOM physics are in SI units.

Related parameters in other subsections:

* the time-stepping ``method`` of the :doc:`simulation_control` subsection must be ``bdf1``, ``bdf2`` or ``bdf3``, and ``bdf1`` if the growth, the nucleation or the coalescence is enabled.

.. warning::
  The steady and SDIRK methods are not supported by this auxiliary physics.

* the linear systems are solved with the parameters of the ``qmom`` subsection of the :doc:`linear_solver_control` subsection, with the GMRES method and an ``ilu`` or ``amg`` preconditioner. The tolerance of each scalar and moment is the ``relative residual`` times its initial residual, but at least the ``minimum residual``, which should therefore be much smaller than the initial residual of the smallest moment. There is no ``non-linear solver`` subsection for the QMOM physics.

Parameters
----------

The default values of the QMOM parameters are given in the text box below:

.. code-block:: text

  subsection QMOM
    set process              = generic
    set number of scalars    = 0
    set scalar names         = # empty: scalar_0, scalar_1, ...
    set scalar diffusivities = # empty: no diffusion
    set number of moments    = 4
    set splitting scheme     = lie
    set enable nucleation    = false
    set enable coalescence   = false
    set enable growth        = false

    subsection growth
      set model               = constant
      set growth coefficient  = 0
      set Function expression = 0
    end

    subsection nucleation
      set model                     = constant
      set nucleation coefficient    = 0
      set Function expression       = 0
      set nuclei mean size          = 1e-09
      set nuclei standard deviation = 1e-10
    end

    subsection coalescence
      set model                   = constant
      set coalescence coefficient = 0
      set Function expression     = 0
    end

    subsection phase evolution
      set enable adaptive node count = false
      set relative tolerance         = 1e-08
      set absolute tolerance         = 1e-20
    end

    subsection lipid nanoparticle production
      # Operating conditions
      set operating temperature     = 295
      set ion concentration mixture = 1

      # Aqueous phase
      set density aqueous                  = 997
      set molar mass aqueous               = 0.018015
      set mole fraction solubility aqueous = 2e-08
      set permittivity aqueous             = 78.2

      # Organic phase
      set density organic                  = 789
      set molar mass organic               = 0.046069
      set mole fraction solubility organic = 0.0029
      set permittivity organic             = 24.3

      # Mixture and particles
      set viscosity mixture = 0.002
      set shape factor      = 0.5235987755982988

      # Lipid
      set lipid molecular volume = 6.0173e-28
      set lipid molar mass       = 0.38665
      set lipid density          = 1067

      # Nucleation
      set nucleation factor         = 2e+22
      set surface tension           = 0.01
      set nuclei standard deviation = 2e-10

      # Growth
      set growth factor = 11.88

      # Coalescence
      set peg factor        = 6.5e-06
      set surface potential = 0.01
      set hamaker constant  = 1.7e-22
    end

    subsection physical constants
      set boltzmann constant  = 1.380649e-23
      set avogadro number     = 6.02214076e+23
      set elementary charge   = 1.602176634e-19
      set permittivity vacuum = 8.8541878188e-12
    end
  end

General
^^^^^^^

* ``process``: process simulated. The choices are:

  * ``generic``: the scalars and the moments are given by ``number of scalars``, ``scalar names`` and ``number of moments``.
  * ``lipid_nanoparticle_production``: lipid nanoparticle production, whose parameters are given in the ``lipid nanoparticle production`` subsection. The process imposes two scalars, ``ethanol_concentration`` and ``free_lipids_concentration``, and 8 moments (4-node quadrature).

  A process other than ``generic`` ignores ``number of scalars``, ``scalar names`` and ``number of moments``, and stops the simulation with an error if they contradict it.

* ``number of scalars``: number of user-defined scalars transported together with the moments.

* ``scalar names``: comma-separated list of scalars names, used for the output files (e.g. ``set scalar names = ethanol_concentration, free_lipids_concentration``). The names may only contain letters, digits and underscores. The default value is empty, in which case the scalars are named ``scalar_0``, ``scalar_1``, etc.

* ``scalar diffusivities``: comma-separated list of the diffusivities of the scalars, in :math:`\text{m}^{2} \cdot \text{s}^{-1}`, in the same order as the ``scalar names`` (e.g. ``set scalar diffusivities = 1e-9, 0`` for two scalars, of which only the first one diffuses). The default value is empty, in which case the scalars do not diffuse. Otherwise, there must be one non-negative value per scalar, including for a process that imposes the scalars. The moments of the particle size distribution never diffuse.

* ``number of moments``: number :math:`K` of transported moments, which must be even. The quadrature has :math:`N = K/2` nodes. The default value, 4, gives a 2-node quadrature, the smallest one that represents a polydisperse distribution.

* ``splitting scheme``: splitting of the physical-space transport :math:`T` and of the phase-space evolution :math:`P` over a time step :math:`\Delta t`. The choices are:

  * ``lie``: physical-space transport of the scalars and of the moments, independently from each other, over :math:`\Delta t`, then phase-space evolution over :math:`\Delta t`, which ends with the update of the scalars coupled to the particles (first-order accurate).
  * ``strang``: phase-space evolution over :math:`\Delta t / 2`, physical-space transport over :math:`\Delta t`, then phase-space evolution over :math:`\Delta t / 2` (second-order accurate).

.. warning::
  The ``strang`` splitting scheme is not available yet and stops the simulation with an error. See :doc:`../../theory/multiphysics/qmom/qmom_equations` for details on the splitting schemes.

Phase-Space Source Terms
^^^^^^^^^^^^^^^^^^^^^^^^

* ``enable nucleation``, ``enable coalescence`` and ``enable growth``: enable the corresponding phase-space mechanisms. They are disabled by default, so that the moments are only transported.

* ``phase evolution``: options of the phase-space evolution in QMOMKit.

  * ``enable adaptive node count``: whether QMOMKit may reduce the number of quadrature nodes when the moments are close to the boundary of the realizable set, instead of stopping the simulation with an error.
  * ``relative tolerance``: relative tolerance of the adaptive Runge-Kutta integration of the growth of the quadrature nodes in phase space.
  * ``absolute tolerance``: absolute tolerance of this integration, in the unit of the positions of the quadrature nodes. It must be much smaller than these positions, so that the ``relative tolerance`` controls the integration.

Growth
""""""

Growth of the particles in phase space, used if ``enable growth`` is ``true``. The growth rate :math:`G` is the velocity of the particle size :math:`\xi` in phase space, in :math:`\text{m} \cdot \text{s}^{-1}` (see :doc:`../../theory/multiphysics/qmom/qmom_equations`).

* ``model``: model of the growth rate. The choices are:

  * ``constant``: :math:`G = a`, with :math:`a` the ``growth coefficient`` in :math:`\text{m} \cdot \text{s}^{-1}`. All the sizes increase by the same amount.
  * ``linear``: :math:`G = a \xi`, with :math:`a` the ``growth coefficient`` in :math:`\text{s}^{-1}`. All the sizes are multiplied by :math:`\exp(a t)`.
  * ``function``: :math:`G` is the ``Function expression``. This generic model runs slower than the ``constant`` and ``linear`` models, since the expression is evaluated by a parser for every quadrature node, at every stage of the integration and at every degree of freedom: use it only for a growth rate that no other model provides.

* ``growth coefficient``: coefficient :math:`a` of the ``constant`` and ``linear`` models.
* ``Function expression``: growth rate of the ``function`` model, in :math:`\text{m} \cdot \text{s}^{-1}`, as an expression of the particle size ``xi``, in :math:`\text{m}`, and of the time ``t``, in :math:`\text{s}` (e.g. ``set Function expression = 2 * xi``). The expression cannot depend on the position.

Nucleation
""""""""""

Nucleation of the particles in phase space, used if ``enable nucleation`` is ``true``. The nucleation rate :math:`J` is the number of particles created per unit volume and time, in :math:`\text{m}^{-3} \cdot \text{s}^{-1}`. The sizes of the nuclei follow a log-normal distribution (see :doc:`../../theory/multiphysics/qmom/qmom_equations`).

* ``model``: model of the nucleation rate. The choices are:

  * ``constant``: :math:`J = a`, with :math:`a` the ``nucleation coefficient`` in :math:`\text{m}^{-3} \cdot \text{s}^{-1}`.
  * ``function``: :math:`J` is the ``Function expression``. This generic model runs slower than the ``constant`` model, since the expression is evaluated by a parser at every stage of the integration and at every degree of freedom.

* ``nucleation coefficient``: coefficient :math:`a` of the ``constant`` model.
* ``Function expression``: nucleation rate of the ``function`` model, in :math:`\text{m}^{-3} \cdot \text{s}^{-1}`, as an expression of the time ``t``, in :math:`\text{s}` (e.g. ``set Function expression = 1e19 * exp(-t)``). The expression cannot depend on the position.
* ``nuclei mean size``: mean size :math:`m` of the log-normal distribution of the sizes of the nuclei, in :math:`\text{m}`.
* ``nuclei standard deviation``: standard deviation :math:`\sigma` of this distribution, in :math:`\text{m}`. It must be strictly positive: nuclei of a single size would have moments that can only be inverted into a single quadrature node.

Coalescence
""""""""""""

Coalescence of the particles in phase space, used if ``enable coalescence`` is ``true``. The coalescence kernel :math:`\beta` is the rate at which two particles merge into a particle of the same volume, in :math:`\text{m}^{3} \cdot \text{s}^{-1}` (see :doc:`../../theory/multiphysics/qmom/qmom_equations`).

* ``model``: model of the coalescence kernel. The choices are:

  * ``constant``: :math:`\beta = a`, with :math:`a` the ``coalescence coefficient`` in :math:`\text{m}^{3} \cdot \text{s}^{-1}`.
  * ``function``: :math:`\beta` is the ``Function expression``. This generic model runs slower than the ``constant`` model, since the expression is evaluated by a parser for every pair of quadrature nodes, at every stage of the integration and at every degree of freedom.

* ``coalescence coefficient``: coefficient :math:`a` of the ``constant`` model.
* ``Function expression``: coalescence kernel of the ``function`` model, in :math:`\text{m}^{3} \cdot \text{s}^{-1}`, as an expression of the sizes ``xi_1`` and ``xi_2`` of the two particles, in :math:`\text{m}`, and of the time ``t``, in :math:`\text{s}` (e.g. ``set Function expression = 1e-18 * (xi_1 + xi_2) / 1e-7``). The expression cannot depend on the position.

Process Parameters
^^^^^^^^^^^^^^^^^^

* ``lipid nanoparticle production``: process parameters of the production of lipid nanoparticles (LNP), in SI units. 

  * ``operating temperature``: temperature of the mixture, in :math:`\text{K}`.
  * ``ion concentration mixture``: ion concentration in the mixture, in :math:`\text{mol} \cdot \text{m}^{-3}`.
  * ``density aqueous`` and ``density organic``: densities of the aqueous and organic phases, in :math:`\text{kg} \cdot \text{m}^{-3}`.
  * ``molar mass aqueous`` and ``molar mass organic``: molar masses of the aqueous and organic phases, in :math:`\text{kg} \cdot \text{mol}^{-1}`.
  * ``mole fraction solubility aqueous`` and ``mole fraction solubility organic``: mole fraction solubility of the lipid in the pure aqueous and organic phases, between 0 and 1.
  * ``permittivity aqueous`` and ``permittivity organic``: relative permittivities of the aqueous and organic phases.
  * ``viscosity mixture``: dynamic viscosity of the mixture of the aqueous and organic phases, in Pa s.
  * ``shape factor``: volume shape factor of the particles, which is :math:`\pi/6` for spheres.
  * ``lipid molecular volume``, ``lipid molar mass`` and ``lipid density``: molecular volume of the lipid, in :math:`\text{m}^{3}`, its molar mass, in :math:`\text{kg} \cdot \text{mol}^{-1}`, and its density in the particles, in :math:`\text{kg} \cdot \text{m}^{-3}`. The default values are those of cholesterol.
  * ``nucleation factor``: pre-exponential factor :math:`A` of the nucleation rate, in :math:`\text{m}^{-3} \cdot \text{s}^{-1}`.
  * ``surface tension``: surface tension (interfacial energy) :math:`\gamma` of the particles, in :math:`\text{N} \cdot \text{m}^{-1}`.
  * ``nuclei standard deviation``: standard deviation of the size of the nuclei, in :math:`\text{m}`.
  * ``growth factor``: kinetic coefficient :math:`k_g` of the growth rate, in :math:`\text{m}^{3} \cdot \text{kg}^{-1} \cdot \text{s}^{-1}`.
  * ``peg factor``: steric hindrance factor :math:`k_{PEG}` of the PEG lipids, the maximum attachment efficiency of the coalescence, between 0 and 1.
  * ``surface potential``: surface potential of the particles, in :math:`\text{V}`.
  * ``hamaker constant``: Hamaker constant of the particles in the mixture, in :math:`\text{J}`.

  The models that use these parameters are described in :ref:`qmom-lipid-nanoparticle-production`. The composition of the mixture is given by the ethanol concentration. The initial and inlet concentrations of ethanol and free lipids, in :math:`\text{kg} \cdot \text{m}^{-3}`, are given by the ``scalars`` of the initial conditions and of the ``dirichlet`` boundary conditions of the QMOM fields (see :ref:`qmom-functions`).

* ``physical constants``: physical constants, in SI units.

  * ``boltzmann constant``: Boltzmann constant, in J/K.
  * ``avogadro number``: Avogadro number, in 1/mol.
  * ``elementary charge``: elementary charge, in C.
  * ``permittivity vacuum``: vacuum permittivity, in F/m.

.. _qmom-functions:

Functions of the QMOM fields
----------------------------

The initial conditions (see :doc:`initial_conditions`), the values of the ``dirichlet`` boundary conditions (see :doc:`boundary_conditions_multiphysics`), the analytical solution (see :doc:`analytical_solution`) and the source terms (see :doc:`source_term`) of the QMOM fields are all given in the same way. The scalars are given as a list of expressions, while the moments are not given one by one: they are the moments of a particle size distribution. The default values are:

.. code-block:: text

  # QMOM block within initial conditions, boundary conditions, and analytical solution
  subsection QMOM
    subsection scalars
      set Function expression = 0
    end
    subsection distribution
      set type                = none
      set Function expression = 0
    end
  end

* ``scalars``: values of the scalars.

  * ``Function expression``: expressions of the scalars, as functions of space and time, separated by semicolons. There must be one expression per scalar, in the same order as the ``scalar names``. The default value, a single expression, is therefore only valid for one scalar.

  .. note::
      The function expressions are ignored if ``number of scalars`` is 0.

* ``distribution``: particle size distribution from which the moments are computed.

  * ``type``: family of the particle size distribution. The moments of each family are computed analytically by QMOMKit. The choices are:

    * ``none``: no particles. All the moments are zero and the ``Function expression`` of the distribution is ignored.
    * ``normal``: normal distribution, with parameters the mean and the standard deviation of the size, in :math:`\text{m}`. It is only meaningful when the standard deviation is much smaller than the mean, since the normal distribution extends to negative sizes.
    * ``lognormal``: log-normal distribution, with parameters the mean and the standard deviation of the natural logarithm of the size. The median size, in :math:`\text{m}`, is the exponential of the first parameter.
    * ``gamma``: gamma distribution, with parameters the shape and the scale, in :math:`\text{m}`.
    * ``weibull``: Weibull distribution, with parameters the shape and the scale, in :math:`\text{m}`.
    * ``beta``: beta distribution on the interval :math:`[0, 1]`, with shape parameters :math:`\alpha` and :math:`\beta`.
    * ``chi``: chi distribution, with parameter the number of degrees of freedom.
    * ``chi_squared``: chi-squared distribution, with parameter the number of degrees of freedom.

  * ``Function expression``: number density of the particles, **in particles per** :math:`\text{m}^{3}`, followed by the parameters of the distribution, as functions of space and time, separated by semicolons. The number of expressions is therefore 3 for the ``normal``, ``lognormal``, ``gamma``, ``weibull`` and ``beta`` distributions, and 2 for the ``chi`` and ``chi_squared`` distributions.

The moment of order :math:`k` is :math:`M_k = n \, \mathrm{E}[L^k]`, where :math:`n` is the number density and :math:`L` the size of the particles. In particular, :math:`M_0` is the number density. A number density of 0 gives zero moments whatever the distribution.

For example, with ``set number of scalars = 2``, two scalars and a log-normal distribution of :math:`10^{18}` particles per :math:`\text{m}^{3}` with a median size of :math:`100\ \text{nm}`, at the inlet of the domain:

.. code-block:: text

  subsection QMOM
    subsection scalars
      set Function expression = 1; 0
    end
    subsection distribution
      set type                = lognormal
      set Function expression = 1e18; log(1e-7); 0.3
    end
  end
