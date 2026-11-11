============================
Quadrature Method of Moments
============================

This subsection defines the parameters of the quadrature method of moments (QMOM) auxiliary physics, used to solve population balance equations. The physics is enabled with ``set qmom = true`` in the :doc:`multiphysics` subsection and requires Lethe to be compiled with `QMOMKit <https://github.com/polycfd/QMOMKit>`_ (``-DLETHE_WITH_QMOMKIT=ON``).

The QMOM physics transports, in a single finite element system:

* ``number of scalars`` regular scalars;
* the raw moments :math:`M_0, \dots, M_{K-1}` of the particle size distribution, where :math:`K` is the ``number of moments``. The quadrature has :math:`N = \lceil K/2 \rceil` nodes, since :math:`N` nodes are defined by :math:`2N` moments.

All the fields will be discretized with discontinuous Galerkin elements of degree ``qmom degree`` (see :doc:`fem`). The scalars are named after the ``scalar names`` and the moments are named ``moment_0`` to ``moment_<K-1>``. All the parameters of the QMOM physics are in SI units.

.. warning::

  The QMOM physics is under development. Only its parameters are available: enabling it in the :doc:`multiphysics` subsection stops the simulation with an error, since its solver is not implemented yet.

Parameters
----------

The default values of the QMOM parameters are given in the text box below:

.. code-block:: text

  subsection QMOM
    set number of scalars  = 0
    set scalar names       = # empty: scalar_0, scalar_1, ...
    set number of moments  = 4
    set splitting scheme   = lie
    set enable nucleation  = false
    set enable coalescence = false
    set enable growth      = false

    subsection phase evolution
      set enable adaptive node count = false
    end

    subsection pharma
      # Operating conditions
      set operating temperature         = 295
      set flow rate ratio               = 3
      set total flow rate               = 2e-07
      set mixing time scale             = 0
      set initial concentration organic = 1
      set ion concentration mixture     = 1

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
    end

    subsection physical constants
      set boltzmann constant  = 1.380649e-23
      set avogadro number     = 6.02214076e+23
      set elementary charge   = 1.602176634e-19
      set permittivity vacuum = 8.8541878188e-12
    end
  end

* ``number of scalars``: number of regular scalars transported together with the moments.

* ``scalar names``: comma-separated names of the scalars, used for the output files (e.g. ``set scalar names = ethanol_concentration, free_lipids_concentration``). There must be one name per scalar, and the names may only contain letters, digits and underscores. The default value is empty, in which case the scalars are named ``scalar_0``, ``scalar_1``, etc.

* ``number of moments``: number :math:`K` of transported raw moments. The quadrature has :math:`N = \lceil K/2 \rceil` nodes: when :math:`K` is odd, one more node is used, so that all the moments are used by the quadrature. The default value, 4, gives a 2-node quadrature, the smallest one that represents a polydisperse distribution.

* ``splitting scheme``: splitting of the physical-space transport :math:`T` and of the phase-space evolution :math:`P` over a time step :math:`\Delta t`. The choices are:

  * ``lie``: :math:`T(\Delta t)` followed by :math:`P(\Delta t)`. This scheme is first-order accurate in time.
  * ``strang``: :math:`P(\Delta t / 2)`, then :math:`T(\Delta t)`, then :math:`P(\Delta t / 2)`. This scheme is second-order accurate in time when both parts are at least second-order accurate.

* ``enable nucleation``, ``enable coalescence`` and ``enable growth``: enable the corresponding phase-space sources. They are disabled by default, so that the moments are only transported.

* ``phase evolution``: options of the phase-space evolution in QMOMKit.

  * ``enable adaptive node count``: whether QMOMKit may reduce the number of quadrature nodes when the moments are close to the boundary of the realizable set.

* ``pharma``: process parameters of the LNP formulation, in SI units. The default values correspond to a water/ethanol formulation in a microfluidic mixer at room temperature: a flow rate ratio of 3 and a total flow rate of 12 mL/min.

  * ``operating temperature``: temperature of the mixture, in K.
  * ``flow rate ratio``: ratio of the aqueous phase flow rate to the organic phase flow rate.
  * ``total flow rate``: total flow rate of the aqueous and organic phases, in m\ :sup:`3`/s.
  * ``mixing time scale``: time scale of the mixing of the aqueous and organic phases, in s.
  * ``initial concentration organic``: lipid concentration in the organic phase before mixing, in mol/m\ :sup:`3`.
  * ``ion concentration mixture``: ion concentration in the mixture, in mol/m\ :sup:`3`.
  * ``density aqueous`` and ``density organic``: densities of the aqueous and organic phases, in kg/m\ :sup:`3`.
  * ``molar mass aqueous`` and ``molar mass organic``: molar masses of the aqueous and organic phases, in kg/mol.
  * ``mole fraction solubility aqueous`` and ``mole fraction solubility organic``: mole fraction solubility of the lipid in the pure aqueous and organic phases, between 0 and 1.
  * ``permittivity aqueous`` and ``permittivity organic``: relative permittivities of the aqueous and organic phases.
  * ``viscosity mixture``: dynamic viscosity of the mixture of the aqueous and organic phases, in Pa s.
  * ``shape factor``: volume shape factor of the particles, which is :math:`\pi/6` for spheres.

* ``physical constants``: physical constants, in SI units. Their default values are the exact values of the 2019 redefinition of the SI base units and the CODATA 2022 value of the vacuum permittivity.

  * ``boltzmann constant``: Boltzmann constant, in J/K.
  * ``avogadro number``: Avogadro number, in 1/mol.
  * ``elementary charge``: elementary charge, in C.
  * ``permittivity vacuum``: vacuum permittivity, in F/m.

Functions of the QMOM fields
----------------------------

The initial conditions (see :doc:`initial_conditions`), the values of the ``dirichlet`` boundary conditions (see :doc:`boundary_conditions_multiphysics`), the analytical solution (see :doc:`analytical_solution`) and the source terms (see :doc:`source_term`) of the QMOM fields are all given in the same way. The scalars are given as a list of expressions, while the moments are not given one by one: they are the raw moments of a particle size distribution. The default values are:

.. code-block:: text

  subsection QMOM
    set Function expression = 0
    subsection distribution
      set type                = none
      set Function expression = 0
    end
  end

* ``Function expression``: expressions of the scalars, as functions of space and time, separated by semicolons. Give either one expression per scalar, or a single expression, which is applied to every scalar.

* ``distribution``: particle size distribution from which the moments are computed.

  * ``type``: family of the particle size distribution. The raw moments of each family are computed analytically by QMOMKit. The choices are:

    * ``none``: no particles. All the moments are zero and the ``Function expression`` of the distribution is ignored.
    * ``normal``: normal distribution, with parameters the mean and the standard deviation of the size, in m. It is only meaningful when the standard deviation is much smaller than the mean, since the normal distribution extends to negative sizes.
    * ``lognormal``: log-normal distribution, with parameters the mean and the standard deviation of the natural logarithm of the size. The median size, in m, is the exponential of the first parameter.
    * ``gamma``: gamma distribution, with parameters the shape and the scale, in m.
    * ``weibull``: Weibull distribution, with parameters the shape and the scale, in m.
    * ``beta``: beta distribution on the interval :math:`[0, 1]`, with shape parameters :math:`\alpha` and :math:`\beta`.
    * ``chi``: chi distribution, with parameter the number of degrees of freedom.
    * ``chi_squared``: chi-squared distribution, with parameter the number of degrees of freedom.

  * ``Function expression``: number density of the particles, **in particles per m**\ :sup:`3`, followed by the parameters of the distribution, as functions of space and time, separated by semicolons. The number of expressions is therefore 3 for the ``normal``, ``lognormal``, ``gamma``, ``weibull`` and ``beta`` distributions, and 2 for the ``chi`` and ``chi_squared`` distributions.

The moment of order :math:`k` is :math:`M_k = n \, \mathrm{E}[L^k]`, where :math:`n` is the number density and :math:`L` the size of the particles. In particular, :math:`M_0` is the number density. A number density of 0 gives zero moments whatever the distribution.

For example, with ``set number of scalars = 2``, two scalars and a log-normal distribution of :math:`10^{18}` particles per m\ :sup:`3` with a median size of 100 nm, at the inlet of the domain:

.. code-block:: text

  subsection QMOM
    set Function expression = 1; 0
    subsection distribution
      set type                = lognormal
      set Function expression = 1e18; log(1e-7); 0.3
    end
  end
