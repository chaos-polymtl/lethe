=================================
Quadrature Method of Moments
=================================

This subsection defines the parameters of the quadrature method of moments (QMOM) auxiliary physics, which models the formation of lipid nanoparticles. The physics is enabled with ``set qmom = true`` in the :doc:`multiphysics` subsection and requires Lethe to be compiled with `QMOMKit <https://github.com/polycfd/QMOMKit>`_ (``-DLETHE_WITH_QMOMKIT=ON``).

The QMOM physics transports, in a single finite element system:

* the ethanol concentration :math:`c_e`;
* the free lipids concentration :math:`c_l`;
* the raw moments :math:`M_0, \dots, M_{2N-1}` of the particle size distribution, where :math:`N` is the number of quadrature nodes. With :math:`N = 4`, the physics has ten fields.

In the output files, these fields are named ``ethanol_concentration``, ``free_lipids_concentration`` and ``moment_0`` to ``moment_<2N-1>``.

Each time step of the fluid dynamics is split into two parts:

* **Transport** (:math:`T`): every field :math:`\phi` is advected by the fluid velocity :math:`\mathbf{u}`,

  .. math::
    \partial_t \phi + \nabla \cdot (\mathbf{u} \phi) = 0,

  with an upwind discontinuous Galerkin discretization of degree ``qmom degree`` (see :doc:`fem`). It is advanced with an explicit strong-stability-preserving Runge-Kutta (SSP-RK) scheme. The transport is sub-cycled within the time step of the fluid dynamics, so that the CFL condition of the explicit scheme is satisfied, and the fluid velocity is interpolated linearly in time between the beginning and the end of the time step.

* **Phase-space evolution** (:math:`P`): nucleation, growth and coalescence change the state :math:`(c_e, c_l, M_0, \dots, M_{2N-1})` independently at every support point of the finite element space. This step is delegated to QMOMKit.

.. warning::

  The QMOM physics is under development. The phase-space sources (nucleation, growth and coalescence) are not implemented yet, so the phase-space step leaves the fields unchanged, and the concentrations are not diffused. The parameters of the sources below are declared but not used yet.

The QMOM physics requires a transient BDF time-stepping method (``bdf1``, ``bdf2`` or ``bdf3``) in the :doc:`simulation_control` subsection, a single fluid, and SI units (no ``dimensionality`` subsection). Its initial conditions, boundary conditions and analytical solution are set in :doc:`initial_conditions`, :doc:`boundary_conditions_multiphysics` and :doc:`analytical_solution`, and the properties of the aqueous and organic phases in :doc:`physical_properties`.

The default values of the QMOM parameters are given in the text box below:

.. code-block:: text

  subsection QMOM
    set number of nodes  = 1
    set splitting scheme = lie

    subsection time integration
      set scheme      = ssprk3
      set maximum cfl = 0.5
    end

    subsection phase evolution
      set enable adaptive node count = false
    end

    # Phase-space sources
    set enable nucleation  = true
    set enable coalescence = true
    set enable growth      = true

    # Operating conditions
    set operating temperature         = 295
    set flow rate ratio               = 1
    set total flow rate               = 1
    set mixing time scale             = 0
    set initial concentration organic = 1
    set ion concentration mixture     = 1

    # Physical constants
    set boltzmann constant  = 1.38065e-23
    set avogadro number     = 6.02214e+23
    set elementary charge   = 1.60218e-19
    set permittivity vacuum = 8.85419e-12
    set shape factor        = 0.523599
  end

* ``number of nodes``: number of nodes :math:`N` of the quadrature. The first :math:`2N` raw moments of the particle size distribution are transported.

* ``splitting scheme``: splitting of the transport :math:`T` and of the phase-space evolution :math:`P` over a time step :math:`\Delta t`. The choices are:

  * ``lie``: :math:`T(\Delta t)` followed by :math:`P(\Delta t)`. This scheme is first-order accurate in time.
  * ``strang``: :math:`P(\Delta t / 2)`, then :math:`T(\Delta t)`, then :math:`P(\Delta t / 2)`. This scheme is second-order accurate in time when both parts are at least second-order accurate.

* ``time integration``: explicit time integration of the transport.

  * ``scheme``: SSP-RK scheme. The choices are ``ssprk1`` (forward Euler, first order), ``ssprk2`` (Heun's method, second order) and ``ssprk3`` (Shu-Osher scheme, third order).
  * ``maximum cfl``: maximum value of the discontinuous Galerkin CFL number :math:`(2p+1)\, |\mathbf{u}|\, \Delta t_s / h` of a transport sub-step of duration :math:`\Delta t_s`, where :math:`p` is the ``qmom degree`` and :math:`h` the minimal distance between two vertices of a cell. At every time step, the transport is divided into the smallest number of sub-steps that satisfies this bound, given the largest velocity at the beginning and at the end of the time step. The default value is a conservative choice that remains stable in three dimensions.

* ``phase evolution``: options of the phase-space evolution in QMOMKit.

  * ``enable adaptive node count``: whether QMOMKit may reduce the number of quadrature nodes when the moments are close to the boundary of the realizable set.

* ``enable nucleation``, ``enable coalescence`` and ``enable growth``: enable the corresponding phase-space sources.

* Operating conditions:

  * ``operating temperature``: temperature of the mixture in K.
  * ``flow rate ratio``: ratio of the flow rates of the aqueous and organic phases.
  * ``total flow rate``: total flow rate of the mixture.
  * ``mixing time scale``: time scale of the mixing of the aqueous and organic phases.
  * ``initial concentration organic``: lipid concentration in the organic phase before mixing.
  * ``ion concentration mixture``: ion concentration in the mixture.

* Physical constants, in SI units: the Boltzmann constant (``boltzmann constant``), the Avogadro number (``avogadro number``), the elementary charge (``elementary charge``), the vacuum permittivity (``permittivity vacuum``), and the volume shape factor of the particles (``shape factor``), which is :math:`\pi/6` for spheres.

.. tip::

  The ``verbosity`` of the ``qmom`` subsection of the :doc:`linear_solver_control` prints the number of transport sub-steps and the CFL number of every time step.
