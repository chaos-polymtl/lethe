==============
Dimensionality
==============

This subsection enables automatic rescaling of physical properties when the simulation is not carried out in SI units. It is particularly useful when a mesh is built in non-SI units (e.g. centimetres or millimetres) but the physical properties (e.g. viscosity, thermal conductivity, permittivity) are known and tabulated in SI units. Instead of converting every property by hand, the user specifies the fundamental unit scales and Lethe rescales all properties automatically.

.. code-block:: text

  subsection dimensionality
    set length           = 1   # metre
    set mass             = 1   # kilogram
    set time             = 1   # second
    set temperature      = 1   # Kelvin
    set electric current = 1   # Ampere
  end

The default value of every parameter is ``1``, which corresponds to the SI
unit system. No rescaling is applied in that case.

Parameters
----------

``length``
  The length unit used in the simulation, expressed in metres. If the mesh
  was built in centimetres, for example, set ``length = 0.01``.

``mass``
  The mass unit used in the simulation, expressed in kilograms. For most
  simulations this can be left at its default value of ``1``.

``time``
  The time unit used in the simulation, expressed in seconds.

``temperature``
  The temperature unit used in the simulation, expressed in Kelvin.

``electric current``
  The electric current unit used in the simulation, expressed in Amperes.
  This parameter is only relevant for simulations that include
  electromagnetic physics (e.g. the time-harmonic Maxwell solver).

How rescaling works
-------------------

Each physical property has a known dimensional formula in terms of the five
base dimensions: length :math:`\mathsf{L}`, mass :math:`\mathsf{M}`, time
:math:`\mathsf{T}`, temperature :math:`\mathsf{\Theta}`, and electric current
:math:`\mathsf{I}`. When non-SI base units are chosen, a quantity expressed in
SI must be multiplied by a *scaling factor* to obtain its value in the
simulation unit system.

If the user sets a length scale :math:`L`, a mass scale :math:`M`, a time
scale :math:`T`, a temperature scale :math:`\Theta`, and an electric current
scale :math:`I` (all expressed in their respective SI units), then a property
whose SI unit is :math:`\mathsf{L}^{a}\,\mathsf{M}^{b}\,\mathsf{T}^{c}\,\mathsf{\Theta}^{d}\,\mathsf{I}^{e}`
is rescaled by

.. math::

   \text{scaling} = L^{a} \cdot M^{b} \cdot T^{c} \cdot \Theta^{d} \cdot I^{e}

and the value supplied by the user (in SI) is multiplied by this factor before
being used in the solver.

The table below lists the scaling factor applied to each property supported
by Lethe, together with its SI unit and its dimensional formula.

.. list-table:: Scaling factors for supported physical properties
   :widths: 28 22 50
   :header-rows: 1

   * - Property
     - SI unit
     - Scaling factor
   * - Density
     - :math:`\text{kg/m}^3`
     - :math:`L^3 / M`
   * - Specific gas constant
     - :math:`\text{J/(kg·K)}`
     - :math:`T^2 \Theta / L^2`
   * - Dynamic viscosity
     - :math:`\text{Pa·s}`
     - :math:`T / L^2`
   * - Specific heat capacity
     - :math:`\text{J/(kg·K)}`
     - :math:`T^2 \Theta / L^2`
   * - Thermal conductivity
     - :math:`\text{W/(m·K)}`
     - :math:`T^3 \Theta / (M \cdot L)`
   * - Specific enthalpy
     - :math:`\text{J/kg}`
     - :math:`T^2 / (M \cdot L^2)`
   * - Diffusivity
     - :math:`\text{m}^2/\text{s}`
     - :math:`T / L^2`
   * - Thermal expansion coefficient
     - :math:`\text{K}^{-1}`
     - :math:`\Theta`
   * - Surface tension
     - :math:`\text{N/m}`
     - :math:`T^2 / M`
   * - Surface tension gradient
     - :math:`\text{N/(m·K)}`
     - :math:`\Theta \cdot T^2 / M`
   * - Cahn–Hilliard mobility
     - :math:`\text{m}^3\text{·s/kg}`
     - :math:`M / (L^3 \cdot T)`
   * - Cahn–Hilliard interface thickness
     - :math:`\text{m}^{-1}`
     - :math:`1 / L`
   * - Electromagnetic frequency
     - :math:`\text{s}^{-1}`
     - :math:`L / c_0` (see note below)
   * - Electric field amplitude (E)
     - :math:`\text{V/m}`
     - :math:`I \cdot T^3 / (M \cdot L)`
   * - Magnetic field amplitude (H)
     - :math:`\text{A/m}`
     - :math:`L / I`
   * - Vacuum permittivity :math:`\varepsilon_0`
     - :math:`\text{F/m}`
     - :math:`M \cdot L^3 / (I^2 \cdot T^4)`
   * - Vacuum permeability :math:`\mu_0`
     - :math:`\text{H/m}`
     - :math:`T^2 \cdot I^2 / (M \cdot L)`

.. note::

   The electromagnetic frequency scaling uses a special formula,
   :math:`L / c_0`, where :math:`c_0 = 299\,792\,458\ \text{m/s}` is the
   speed of light in vacuum. The time scale does not appear in this factor because the time-harmonic Maxwell formulation is time-independent by construction.

.. warning:: 

   If the user chooses to do the rescaling by hand, make sure to adjust the frequency so the dimensionless wavenumber :math:`\kappa = \omega \sqrt{\varepsilon \mu}L` remains unchanged when changing the length scale :math:`L`. This is necessary to preserve the physics of the problem as the electromagnetic solver always solves the system in dimensionless form and does not know about the length scale.

Examples
--------

Mesh in centimetres, properties in SI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A mesh is built in centimetres so that a physical height of :math:`0.10\ \text{m}`
appears as :math:`10` in the mesh file. All physical properties (viscosity,
thermal conductivity, etc.) are available in SI units.

.. code-block:: text

  subsection dimensionality
    set length = 0.01
  end

With this setting, a kinematic viscosity of :math:`1 \times 10^{-6}\ \text{m}^2/\text{s}`
is automatically rescaled by :math:`T/L^2 = 1/(0.01)^2 = 10^4`, yielding the
equivalent value :math:`10^{-2}\ \text{cm}^2/\text{s}` inside the solver.

Electromagnetic simulation in centimetres
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A time-harmonic Maxwell simulation is carried out on a centimetre-scale
geometry. The electric field, permittivity, and permeability are all rescaled
consistently.

.. code-block:: text

  subsection dimensionality
    set length = 0.01
  end

The key scalings are:

* **Electric field** :math:`[\text{V/m}]`: scaling :math:`= I T^3/(M L) = 1/(0.01) = 100`.
  A field of :math:`1\ \text{V/m}` becomes :math:`100` in simulation units (which are not :math:`[\text{V/m}]` anymore, but rather :math:`[\text{V_{cm}/cm}]`, a custom basis for the voltage).
* **Magnetic field** :math:`[\text{A/m}]`: scaling :math:`= L/I = 0.01/1 = 0.01`.
  :math:`H = 1\ \text{A/m}` becomes :math:`0.01` in simulation units.

