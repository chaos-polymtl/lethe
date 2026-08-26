..
    SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
    SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

Microwave Heating
=========================

This example couples the time-harmonic Maxwell solver with heat transfer to simulate the **microwave heating** of a dielectric cylinder sitting inside a rectangular waveguide. A first case reproduces, without any fluid flow, the resonance-driven heating of a low-loss ceramic cylinder studied by Peng *et al.* [#Peng2024]_. A second case builds on the first one and shows how a fluid flow can be added to the same setup to convectively cool the heated object.

Features
--------

- Solver: ``lethe-fluid`` or ``lethe-fluid-matrix-free``
- Transient problem
- Coupling of the time-harmonic Maxwell solver with heat transfer through the ``microwave heating`` source term
- Excitation of a rectangular waveguide through a ``waveguide port`` boundary condition and absorption of outgoing waves through a matched ``impedance boundary``
- Use of the built-in ``uniform_channel_with_meshed_cylinder`` grid, in which a solid cylinder is meshed and embedded in a rectangular channel
- Comparison with literature results
- Addition of a fluid flow, of a local box refinement, and of checkpointing to the base case

Files Used in This Example
--------------------------

Both parameter files below are located in the example's folder (``examples/multiphysics/microwave-heating``).

- Parameter file, static heating of an alumina cylinder: ``filled_waveguide_cylinder_Al.prm``
- Parameter file, heating of a silicon carbide cylinder cooled by an air flow: ``filled_waveguide_cylinder_SiC.prm``

Description of the Case
-----------------------

Geometry
~~~~~~~~

Both cases share the same base geometry: a section of rectangular waveguide, with a cylindrical dielectric sample standing across its narrow dimension, meshed with the built-in ``uniform_channel_with_meshed_cylinder`` grid (see the :ref:`Uniform Channel with Meshed Cylinder <channel-cylinder>` documentation for the full description of this grid and its arguments).

.. image:: images/schematic.png
    :alt: schematic of the filled waveguide
    :align: center
    :name: schematic
    :width: 500

The waveguide's cross-section is :math:`109.2\ \mathrm{mm}\times54.6\ \mathrm{mm}` in both cases (a standard WR-430 rectangular waveguide), with the cylinder's axis parallel to the narrow (:math:`54.6\ \mathrm{mm}`) dimension, i.e., parallel to the dominant electric field of the fundamental :math:`\mathrm{TE}_{10}` mode. The cylinder itself is meshed (``mesh_obstacle = true``) so that the temperature field can be resolved inside it, and it is surrounded by a structured transition region before blending into the padded background mesh of the channel.

- In ``filled_waveguide_cylinder_Al.prm``, the channel is :math:`200\ \mathrm{mm}` long, and the cylinder has a radius of :math:`24\ \mathrm{mm}`, centered midway along the channel.
- In ``filled_waveguide_cylinder_SiC.prm``, the channel is :math:`400\ \mathrm{mm}` long to leave room for the flow to develop, and the cylinder has a radius of :math:`20\ \mathrm{mm}`.

Physical Problem
~~~~~~~~~~~~~~~~

As in the :doc:`waveguide example <../waveguide/waveguide>`, the time-harmonic Maxwell solver excites a single rectangular waveguide mode (here, the fundamental :math:`\mathrm{TE}_{10}` mode) through a ``waveguide port`` boundary condition, and absorbs the wave transmitted past the cylinder through a matched ``impedance boundary`` condition, which mimics a semi-infinite waveguide by preventing spurious reflections back toward the cylinder.

Since ``set microwave heating = true`` is used in both cases, the power dissipated by the dielectric losses of the cylinder is automatically computed from the electromagnetic solution and added as a source term to the heat transfer equation,

.. math::
    Q_\mathrm{em} = \frac{1}{2}\sigma|\mathbf{E}|^2 + \frac{1}{2}\omega\varepsilon_0\varepsilon_\mathrm{im}|\mathbf{E}|^2 + \frac{1}{2}\omega\mu_0\mu_\mathrm{im}|\mathbf{H}|^2,

see the :doc:`multiphysics <../../../parameters/cfd/multiphysics>` documentation for details. Since both materials considered here are non-magnetic and non-conductive (:math:`\mu_\mathrm{im}=\sigma=0`), only the dielectric loss term, proportional to :math:`\varepsilon_\mathrm{im}|\mathbf{E}|^2`, contributes to the heating.

.. tip::
    Because none of the physical properties used in this example depend on temperature, the electromagnetic fields do not need to be recomputed as the cylinder heats up. Both parameter files therefore rely on (or default to) ``subsection time coupling strategy`` with ``set type = none``: the electromagnetic problem is solved once, before the first time step, and the resulting heat source is then reused throughout the transient heat transfer solve. See the :doc:`../../../parameters/cfd/time_harmonic_maxwell` documentation for the other available coupling strategies, needed when the physical properties depend on the temperature.

Resonance-Driven Heating
~~~~~~~~~~~~~~~~~~~~~~~~

The amount of power absorbed by a dielectric cylinder placed in a waveguide is strongly dependent on its radius and permittivity: for specific combinations of these two parameters, the internal electromagnetic field can build up through constructive interference (a phenomenon closely related to Mie resonances in scattering theory), leading to a dramatically enhanced, and highly non-uniform, heating rate compared to what a simple volume-averaged absorption estimate would suggest. This effect is especially significant for **low-loss** materials, such as alumina, since their long internal photon lifetime allows the field to build up before being dissipated.

The ``filled_waveguide_cylinder_Al.prm`` case reproduces this behavior for a low-loss alumina cylinder and is meant to be compared against the results of Peng *et al.* [#Peng2024]_, who studied this exact resonance-driven heating mechanism, both theoretically (using Mie theory) and numerically, for low-loss cylindrical samples in a waveguide.

Parameter Files
---------------

Case 1: Static Heating of an Alumina Cylinder
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This case (``filled_waveguide_cylinder_Al.prm``) solves the coupled electromagnetics/heat transfer problem without any fluid flow.

Simulation Control
^^^^^^^^^^^^^^^^^^

.. code-block:: text

    subsection simulation control
        set method                         = bdf1
        set time step                      = 0.0001
        set adapt time step to respect CFL = true
        set adaptative time step scaling   = 1.02
        set time end                       = 60
        set output time frequency          = 0.05
        set output path                    = ./cylinder_resonance/
        set subdivision                    = 3
    end

The heat transfer equation is integrated in time with a first-order backward difference scheme (``bdf1``) up to :math:`t=60\ \mathrm{s}`. Since ``fluid dynamics = false`` (see below), no velocity field is solved for and the CFL condition is trivially satisfied at every time step; the time step therefore grows freely, by a factor of :math:`1.02` at each iteration, from its small initial value of :math:`10^{-4}\ \mathrm{s}`, allowing the transient to be resolved efficiently even over a fairly long physical time.

Mesh Adaptation
^^^^^^^^^^^^^^^

.. code-block:: text

    subsection mesh adaptation
        set type = none
    end

No adaptive mesh refinement is used in this example; the mesh resolution is instead controlled entirely by the ``mesh`` subsection below.

Multiphysics
^^^^^^^^^^^^

.. code-block:: text

    subsection multiphysics
        set fluid dynamics    = false
        set microwave heating = true
        set electromagnetics  = true
        set heat transfer     = true
    end

The fluid dynamics solver is disabled: only the electromagnetics and heat transfer physics are solved, coupled through the ``microwave heating`` source term described above.

Mesh
^^^^

.. code-block:: text

    subsection mesh
        set type               = lethe
        set grid type          = uniform_channel_with_meshed_cylinder
        set grid arguments     = 0., 0. : 0.2, 0.1092 : 0.1, 0.0546 : 0.024 : 0.04 : 1 : 1 : 2 : 2 : 0.0546 : 2 : false : true : true
        set initial refinement = 2
    end

This builds the :math:`200\ \mathrm{mm} \times 109.2\ \mathrm{mm}` channel cross-section extruded to a height of :math:`54.6\ \mathrm{mm}`, with a meshed cylinder of radius :math:`24\ \mathrm{mm}` (``mesh_obstacle = true``) centered at mid-length and mid-height of the channel, and colorized (``colorize = true``) boundaries.

Time Harmonic Maxwell
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

    subsection time harmonic maxwell
        set electromagnetic frequency    = 2.45e9
        set number of waveguide inlets   = 1
        set electromagnetic scaling type = power

        subsection waveguide inlet 0
            set port boundary id = 0

            set waveguide power = 50

            set corner 0 = 0,0,0
            set corner 1 = 0.,0.1092,0
            set corner 2 = 0., 0., 0.1092
            set corner 3 = 0., 0.1092, 0.1092

            subsection waveguide mode
                set mode type    = TE
                set mode order m = 1
                set mode order n = 0
            end
        end
    end

The waveguide is excited at the standard industrial, scientific and medical (ISM) microwave frequency of :math:`2.45\ \mathrm{GHz}` with a :math:`\mathrm{TE}_{10}` mode (:math:`m=1`, :math:`n=0`), fed through boundary id ``0`` (the left face of the channel) at a power of :math:`50\ \mathrm{W}`. Since ``electromagnetic scaling type = power`` is used, the solution is rescaled after solving so that the power flowing through this inlet matches this value.

Boundary Conditions
^^^^^^^^^^^^^^^^^^^

.. attention::
    As in the :doc:`waveguide <../waveguide/waveguide>` and :doc:`Fichera oven <../fichera-oven/fichera-oven>` examples, the ``subsection boundary conditions`` for the fluid dynamics boundaries cannot be removed from the parameter file, even though the fluid solver is disabled:

    .. code-block:: text

        subsection boundary conditions
            set number = 1
            subsection bc 0
                set id   = 0, 1, 2, 3, 4, 5
                set type = noslip
            end
        end

.. code-block:: text

    subsection boundary conditions time harmonic maxwell
        set number = 3
        subsection bc 0
            set id   = 2, 3, 4, 5
            set type = pec
        end
        subsection bc 1
            set id   = 0
            set type = waveguide port
        end
        subsection bc 2
            set id   = 1
            set type = impedance boundary
            subsection excitation x real part
                set Function expression = 0
            end
            subsection excitation x imag part
                set Function expression = 0
            end
            subsection excitation y real part
                set Function expression = 0
            end
            subsection excitation y imag part
                set Function expression = 0
            end
            subsection excitation z real part
                set Function expression = 0
            end
            subsection excitation z imag part
                set Function expression = 0
            end
            subsection surface admittance real part
                set Function expression = 0.828306014816808
            end
            subsection surface admittance imag part
                set Function expression = 0.
            end
        end
    end

- ``bc 0`` (ids ``2`` to ``5``, the four channel walls) is a ``pec`` boundary.
- ``bc 1`` (id ``0``, the left face) is the ``waveguide port`` used to excite the :math:`\mathrm{TE}_{10}` mode.
- ``bc 2`` (id ``1``, the right face) is an ``impedance boundary`` whose surface admittance is matched to the :math:`\mathrm{TE}_{10}` wave admittance of the empty guide, with no additional excitation, so that the wave transmitted past the cylinder is absorbed rather than reflected back.

.. code-block:: text

    subsection boundary conditions heat transfer
        set number = 1
        subsection bc 0
            set id   = 0, 1, 2, 3, 4, 5
            set type = noflux
        end
    end

Every wall of the channel is insulated (``noflux``): in this first case, the cylinder can only lose heat by conduction through the surrounding stagnant air, which stays fixed for the whole simulation since ``fluid dynamics = false``.

FEM
^^^

.. code-block:: text

    subsection FEM
        set temperature degree            = 3
        set electromagnetics trial degree = 2
        set electromagnetics test degree  = 3
    end

Physical Properties
^^^^^^^^^^^^^^^^^^^

.. code-block:: text

    subsection physical properties
        set number of fluids = 1
        set number of solids = 1
        subsection fluid 0
            set electric conductivity model = constant
            set electric conductivity       = 0.

            set electric permittivity model     = constant
            set electric permittivity real part = 1.
            set electric permittivity imag part = 0.

            set magnetic permeability model     = constant
            set magnetic permeability real part = 1.
            set magnetic permeability imag part = 0.

            set kinematic viscosity  = 1.48e-5
            set specific heat        = 1006
            set density              = 1.225
            set thermal conductivity = 2.6e-2
        end

        subsection solid 0
            set electric conductivity model = constant
            set electric conductivity       = 0.

            set electric permittivity model     = constant
            set electric permittivity real part = 9.2
            set electric permittivity imag part = 0.005

            set magnetic permeability model     = constant
            set magnetic permeability real part = 1.
            set magnetic permeability imag part = 0.

            set thermal conductivity = 26
            set specific heat        = 1046
            set density              = 3750
        end
    end

The waveguide is filled with air (``fluid 0``), described here only for the sake of its electromagnetic and thermal properties since it does not flow in this first case. The cylinder (``solid 0``) is alumina (:math:`\mathrm{Al_2O_3}`): a low-loss ceramic with a relative permittivity of :math:`\varepsilon_r \approx 9.2 - 0.005i`, i.e., a very small loss tangent of about :math:`5\times10^{-4}`, consistent with the low-loss regime studied by Peng *et al.* [#Peng2024]_.

Non-Linear and Linear Solver Control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

    subsection non-linear solver
        subsection fluid dynamics
            set verbosity = quiet
        end
        subsection heat transfer
            set verbosity      = verbose
            set tolerance      = 1e-4
            set max iterations = 10
        end
    end

    subsection linear solver
        subsection electromagnetics
            set verbosity         = verbose
            set relative residual = 1e-4
            set minimum residual  = 1e-8
            set preconditioner    = none
        end
        subsection heat transfer
            set max iters          = 1000
            set relative residual  = 1e-3
            set minimum residual   = 1e-6
            set max krylov vectors = 1000
        end
    end

Case 2: Adding Fluid Flow
~~~~~~~~~~~~~~~~~~~~~~~~~

This second case (``filled_waveguide_cylinder_SiC.prm``) starts from the same base setup and adds a laminar air flow along the channel, so that the cylinder is now convectively cooled while being heated by the electromagnetic field. Only the parameters that differ from, or are added to, Case 1 are detailed below.

Multiphysics
^^^^^^^^^^^^

.. code-block:: text

    subsection multiphysics
        set fluid dynamics    = true
        set microwave heating = true
        set electromagnetics  = true
        set heat transfer     = true
    end

The only difference with Case 1 is ``set fluid dynamics = true``: the incompressible Navier-Stokes equations are now solved, in addition to electromagnetics and heat transfer.

Mesh and Box Refinement
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

    subsection mesh
        set type                         = lethe
        set grid type                    = uniform_channel_with_meshed_cylinder
        set grid arguments               = 0., 0. : 40, 10.92 : 10, 5.56 : 2 : 4 : 1 : 1 : 2 : 8 : 5.46 : 3 : false : true : true
        set initial refinement           = 2
        set initial boundary refinement  = 1
        set boundaries refined           = 2,3,4,5
    end

    subsection box refinement
        set number of refinement boxes = 1
        subsection box 0
            subsection mesh
                set type                = dealii
                set grid type           = cylinder_shell
                set grid arguments      = 5.46: 2: 2.3: 20:20:false
                set initial translation = 10, 5.56, 0
            end
            set additional refinement = 2
        end
    end

The channel is now :math:`400\ \mathrm{mm}` long (to leave room for the flow to develop) and the cylinder radius is reduced to :math:`20\ \mathrm{mm}`. All lengths are given in centimeters here rather than meters (see the ``dimensionality`` subsection below), and the four channel walls (ids ``2`` to ``5``) receive one extra level of boundary refinement to better resolve the developing boundary layers.

A ``box refinement`` region, a thin cylindrical shell wrapped around the physical cylinder, is also added, providing two additional levels of refinement in the vicinity of the object to better resolve both the thermal and momentum boundary layers there. See the :doc:`../../../parameters/cfd/box_refinement` documentation for details.

Initial Conditions
^^^^^^^^^^^^^^^^^^

.. code-block:: text

    subsection initial conditions
        set type = nodal
        subsection uvwp
            set Function constants  = c_width = 10.92, c_height = 5.46, Re = 400, scaling = 1.782122724170977, percent = 0.1
            set Function expression = Re*scaling*abs((if(y < percent*c_width/2, y*(y-percent*c_width), if(y < c_width-(percent*c_width/2), (percent*c_width * percent*c_width)/4 , (y-c_width)*(y-c_width+ percent*c_width)))) * (if(z < percent*c_height/2, z*(z-percent*c_height), if(z < c_height-(percent*c_height/2), (percent*c_height * percent*c_height)/4, (z-c_height)*(z-c_height+percent*c_height)))));0;0;0
        end
        subsection temperature
            set Function expression = 0
        end
    end

The velocity field is initialized with an approximate fully-developed laminar profile for a rectangular duct, scaled to reach a Reynolds number of :math:`400` based on the channel's transverse dimensions, and rounded off near the walls (over ``percent = 10%`` of each transverse dimension) to avoid an unphysical discontinuity in the wall-normal velocity gradient at the channel edges. The same expression is reused verbatim as the inlet boundary condition below. The temperature field is initialized at :math:`0` (in the units defined by the ``dimensionality`` subsection).

Boundary Conditions
^^^^^^^^^^^^^^^^^^^

.. code-block:: text

    subsection boundary conditions
        set number = 3
        subsection bc 0
            set id   = 2, 3, 4, 5
            set type = noslip
        end
        subsection bc 1
            set id   = 0
            set type = function
            subsection u
                set Function constants  = c_width = 10.92, c_height = 5.46, Re = 400, scaling = 1.782122724170977, percent = 0.1
                set Function expression = Re*scaling*abs((if(y < percent*c_width/2, y*(y-percent*c_width), if(y < c_width-(percent*c_width/2), (percent*c_width * percent*c_width)/4 , (y-c_width)*(y-c_width+ percent*c_width)))) * (if(z < percent*c_height/2, z*(z-percent*c_height), if(z < c_height-(percent*c_height/2), (percent*c_height * percent*c_height)/4, (z-c_height)*(z-c_height+percent*c_height)))))
            end
            subsection v
                set Function expression = 0
            end
            subsection w
                set Function expression = 0
            end
        end
        subsection bc 2
            set id   = 1
            set type = outlet
        end
    end

    subsection boundary conditions heat transfer
        set number = 2
        subsection bc 0
            set id   = 1, 2, 3, 4, 5
            set type = noflux
        end
        subsection bc 1
            set id   = 0
            set type = temperature
            subsection value
                set Function expression = 0
            end
        end
    end

The left face (id ``0``) is now a velocity inlet, driven by the same laminar profile used to initialize the flow, while the right face (id ``1``) is a pressure ``outlet``. The incoming air is also imposed at a fixed temperature of :math:`0` on that same inlet face, while every other wall remains adiabatic (``noflux``); note that the electromagnetic boundary conditions (not repeated here) are unchanged from Case 1.

Physical Properties and Dimensionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

    subsection physical properties
        set number of fluids = 1
        set number of solids = 1
        subsection fluid 0
            # ... identical to Case 1 (air) ...
        end
        subsection solid 0
            set electric conductivity model = constant
            set electric conductivity       = 0.

            set electric permittivity model     = constant
            set electric permittivity real part = 9.72
            set electric permittivity imag part = 2.01

            set magnetic permeability model     = constant
            set magnetic permeability real part = 1.
            set magnetic permeability imag part = 0.

            set thermal conductivity = 120.92
            set specific heat        = 27.13
            set density              = 3210
        end
    end

    subsection dimensionality
        set length           = 0.01  #meter
        set mass             = 1     #kilogram
        set time             = 1     #second
        set temperature      = 1     #kelvin
        set electric current = 0.01  #ampere
    end

The cylinder is now silicon carbide (SiC), a much lossier ceramic than the alumina used in Case 1 (:math:`\varepsilon_r \approx 9.72 - 2.01i`). Since the mesh is built in centimeters in this case, a ``dimensionality`` subsection is added so that the physical properties above, still entered in SI units, are automatically rescaled to be consistent with the mesh: see the :doc:`../../../parameters/cfd/dimensionality` documentation for a detailed explanation of this rescaling.

Linear Solver Control
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

    subsection linear solver
        subsection fluid dynamics
            set method             = gmres
            set max iters          = 1000
            set relative residual  = 1e-4
            set minimum residual   = 1e-8
            set preconditioner     = gcmg
            set verbosity          = verbose
            set max krylov vectors = 500

            # MG parameters
            set mg verbosity                   = quiet
            set mg enable hessians in jacobian = false
            set mg coarsening type             = ph
            set mg p coarsening type           = decrease by one

            # Smoother
            set mg smoother iterations          = 3
            set mg smoother eig estimation      = true
            set mg smoother preconditioner type = inverse diagonal

            # Eigenvalue estimation parameters
            set eig estimation smoothing range = 20
            set eig estimation cg n iterations = 20
            set eig estimation verbosity       = quiet

            # Coarse-grid solver
            set mg coarse grid solver = direct
        end
        subsection heat transfer
            set max iters          = 1000
            set relative residual  = 1e-4
            set minimum residual   = 1e-8
            set max krylov vectors = 1000
        end
        subsection electromagnetics
            set verbosity         = verbose
            set relative residual = 1e-4
            set minimum residual  = 1e-8
            set preconditioner    = none
        end
    end

The fluid dynamics equations are solved with the matrix-free solver, preconditioned by a Global Coarsening Multigrid (``gcmg``) method. Refer to the :doc:`../../../parameters/cfd/linear_solver_control` documentation for more details on the multigrid parameters.

Restart
^^^^^^^

.. code-block:: text

    subsection restart
        set checkpoint = true
        set frequency  = 200
        set filename   = restart
        set restart    = false
    end

Since this second case is significantly more expensive than Case 1 (it now also solves the Navier-Stokes equations, on a locally refined mesh), checkpointing is enabled so that the simulation can be resumed if it needs to be interrupted. See the :doc:`../../../parameters/cfd/restart` documentation for details.

Running the Simulations
-----------------------

.. code-block:: text
    :class: copy-button

    mpirun -np 8 lethe-fluid filled_waveguide_cylinder_Al.prm

.. code-block:: text
    :class: copy-button

    mpirun -np 8 lethe-fluid filled_waveguide_cylinder_SiC.prm

.. warning::
    Case 2 solves a fully coupled electromagnetics/heat transfer/fluid dynamics problem on a locally refined mesh over a :math:`60\ \mathrm{s}` transient, and is therefore substantially more expensive than Case 1. Running it on a workstation with a modest number of cores may take several hours; the ``restart`` subsection described above allows the simulation to be resumed if needed.

Results and Discussion
----------------------

The following figure shows the steady-state temperature field obtained for the alumina cylinder of Case 1, once the resonant heating has raised its internal temperature:

.. image:: images/temperature_al.png
    :alt: temperature field in the alumina cylinder, without fluid flow
    :align: center
    :name: temperature-al
    :width: 500

Because the cylinder's radius and permittivity place it near a resonant condition of the waveguide-cylinder system, the internal electric field, and therefore the heating rate, are strongly enhanced and highly non-uniform, in agreement with the resonance-driven heating mechanism described by Peng *et al.* [#Peng2024]_.

The following figure shows the corresponding temperature field for the silicon carbide cylinder of Case 2, once the flow field and the temperature field have reached a statistically steady regime:

.. image:: images/temperature_sic.png
    :alt: temperature field in the SiC cylinder, cooled by an air flow
    :align: center
    :name: temperature-sic
    :width: 500

Unlike Case 1, the air flow continuously removes heat from the cylinder by forced convection, which limits its temperature rise and skews the temperature field toward the downstream side of the cylinder.

.. Once available, the animation of the transient heating and flow fields can be embedded here, following the convention used throughout Lethe's documentation, for instance:
.. .. raw:: html
..
..     <p align="center"><iframe width="720" height="405" src="https://www.youtube.com/embed/VIDEO_ID" title="Microwave heating of a cylinder in a filled waveguide" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

Possibilities for Extension
---------------------------

- **Sweep the cylinder radius:** Rerun Case 1 for a range of ``inner_radius`` values (adjusting the ``grid arguments`` accordingly) to reproduce the resonance curve reported by Peng *et al.* [#Peng2024]_, and locate the radius that maximizes the absorbed power for the alumina cylinder.
- **Temperature-dependent properties:** Replace the ``constant`` electric permittivity and conductivity models by temperature-dependent ones, and switch the ``time coupling strategy`` from ``none`` to ``iteration`` or ``threshold`` so that the electromagnetic field is periodically recomputed as the material heats up and its properties drift.
- **Increase the flow rate:** Increase the Reynolds number ``Re`` used in the initial and boundary conditions of Case 2 to study how stronger convective cooling affects the peak temperature reached by the cylinder.
- **Turbulent flow:** For higher Reynolds numbers, a turbulence model may need to be added to Case 2; see the incompressible flow examples for guidance on setting up turbulent simulations in Lethe.

References
----------

.. [#Peng2024] \Y. Peng, D. Zhou, H. Chen, M. Liu, Z. Tang, and T. Hong, "Resonance-Driven Microwave Heating of Low-Loss Cylindrical Substances in Waveguide Systems," *IEEE Transactions on Microwave Theory and Techniques*, vol. 72, no. 6, pp. 3722-3733, June 2024, doi: `10.1109/TMTT.2023.3327488 <https://doi.org/10.1109/TMTT.2023.3327488>`_\.