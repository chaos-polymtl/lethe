=====================
Time-Harmonic Maxwell
=====================

In this subsection, the parameters for the time-harmonic Maxwell solver are specified.

This solver is the DPG formulation used to model electromagnetic wave excitation in
waveguides and to control how the electromagnetic solution is scaled after the linear
system is solved.
The main parameter subsection is ``time harmonic maxwell``.

The nondimensional convention used by the solver is:

.. math::

  \begin{align}
     & \tilde{\mathbf{E}} = \frac{1}{E_0}\mathbf{E}_\mathrm{spatial}, &                                                     & \tilde{\mathbf{H}} = \frac{Z_0}{E_0} \mathbf{H}_\mathrm{spatial},
     &                                                                & \varepsilon_r = \frac{1}{\varepsilon_0}\varepsilon, &                                                                                 & \mu_r = \frac{1}{\mu_0}\mu, \nonumber \\
     & \tilde{\nabla} = L \nabla,                                     &                                                     & \sigma_r = \frac{1}{\omega \varepsilon_0}\sigma,                                &
     & \tilde{\omega} =  \frac{ L}{c_0}\omega,                        &                                                     & \tilde{ \mathbf{J} } = \frac{L Z_0}{E_0} \mathbf{J}_\mathrm{spatial}. \nonumber
  \end{align}

The time-harmonic Maxwell system is written in dimensionless form as:

.. math::

  \begin{subequations}\label{eq:maxwell_dimensionless}
    \begin{align}
      \nabla \times \tilde{\mathbf{E}} - i \tilde{\omega} \mu_r \tilde{\mathbf{H}}                      & = 0,                  \\
      \nabla \times \tilde{\mathbf{H}} + i \tilde{\omega} \varepsilon_{r,\text{eff}} \tilde{\mathbf{E}} & = \tilde{\mathbf{J}},
    \end{align}
  \end{subequations}

where :math:`\tilde{\mathbf{E}}` and :math:`\tilde{\mathbf{H}}` are the dimensionless electric
and magnetic fields, :math:`\tilde{\omega}` is the dimensionless angular frequency,
:math:`\mu_r` is the relative magnetic permeability, :math:`\varepsilon_{r,\text{eff}}` is
the effective relative permittivity, and :math:`\tilde{\mathbf{J}}` is the dimensionless
current density.

The solver always works in this dimensionless formulation and rescales the solution only
through the selected electromagnetic scaling type after the linear system is solved.

The time-harmonic Maxwell solver is currently supported only for 3D problems.

The parameters parsed by the solver are grouped as follows:

.. code-block:: text

  subsection time harmonic maxwell

    subsection time coupling strategy
      set type               = none
      set coupling iteration = 1
      set coupling time      = 1
      set coupling threshold = 0.1
    end

    set electromagnetic frequency   = 1
    set electromagnetic scaling type = none
    set electric field amplitude     = 1
    set magnetic field amplitude     = 1
    set number of waveguide inlets   = 0

    subsection waveguide inlet 0
      set port boundary id = 0
      set waveguide power  = 1

      subsection waveguide mode
        set mode type    = TE
        set mode order m = 1
        set mode order n = 0
      end

      # In 3D, define corner 0, corner 1, corner 2, and corner 3.
      set corner 0 = ...
      set corner 1 = ...
      set corner 2 = ...
      set corner 3 = ...
    end
  end

.. note::

  The ``waveguide inlet`` subsections are only parsed for the first ``number of waveguide inlets`` entries. At most ten inlets are currently supported by the parameter declaration.

  The solver is only available for 3D simulations, so the inlet corner coordinates always have three components.

* ``time coupling strategy``: controls how often the electromagnetic fields are recomputed.

  * ``type``: selects the coupling strategy. Choices are ``none``, ``iteration``, ``time`` or ``threshold``.

  * ``coupling iteration``: number of time iterations between two electromagnetic solves when ``type = iteration``.

  * ``coupling time``: physical time interval between two electromagnetic solves when ``type = time``.

  * ``coupling threshold``: change in the electromagnetic properties of the medium that triggers a new solve when ``type = threshold``.

* ``electromagnetic frequency``: excitation frequency of the time-harmonic wave in Hz. The value is scaled internally using the dimensionality scaling factors.

* ``electromagnetic scaling type``: selects how the solution is rescaled after solving. Choices are ``none``, ``electric field``, ``magnetic field`` or ``power``.

* ``electric field amplitude``: dimensional electric-field amplitude used for normalization of the solution.

* ``magnetic field amplitude``: dimensional magnetic-field amplitude used for normalization of the solution.

* ``number of waveguide inlets``: number of waveguide inlet subsections that will be parsed.

* ``waveguide inlet i``: describes inlet ``i``. Each inlet must define:

  * ``port boundary id``: boundary id where the inlet is applied.

  * ``waveguide power``: excitation power in Watts. This value must be strictly positive.

  * ``waveguide mode``: specifies the rectangular waveguide mode.

    * ``mode type``: either ``TE`` or ``TM``.

    * ``mode order m``: mode order in the first transverse direction.

    * ``mode order n``: mode order in the second transverse direction.

  * ``corner i``: coordinates of the inlet corners. Four corners are expected and the parser checks that they define a coplanar quadrilateral.

.. note::

  The ``mode order m`` and ``mode order n`` values cannot both be zero.

.. attention::

  The ``power`` electromagnetic scaling type requires at least one waveguide inlet to be defined.