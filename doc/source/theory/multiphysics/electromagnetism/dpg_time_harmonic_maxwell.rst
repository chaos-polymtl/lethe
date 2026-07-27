==================================================
DPG formulation for time-harmonic Maxwell problems
==================================================

The discontinuous Petrov-Galerkin (DPG) method is a residual-minimizing finite
element strategy that is particularly effective for wave propagation and
high-frequency electromagnetic problems. Rather than using the same space for
trial and test functions, DPG uses different spaces and seeks the discrete
solution by minimizing the residual of the PDE in a norm induced by the test
space.

For an abstract problem

.. math::

   b(v,u) = l(v) \quad \forall v \in V,

DPG seeks an approximate trial function :math:`u_h \in U_h` that minimizes the
residual in the dual norm of the test space,

.. math::

   u_h = \arg\min_{w_h \in U_h} \frac{1}{2} \| l - B w_h \|_{V'}^2,

where :math:`B : U \to V'` is the operator associated with the bilinear form
:math:`b`. To make this practical, one introduces the Riesz operator
:math:`R_V : V \to V'` and the corresponding error representation function
:math:`\Psi = R_V^{-1}(l - B u_h)`. The quantity :math:`\Psi` is a test-space
function that measures the residual and provides the basis for the DPG
mixed formulation.

In practice, the test space is chosen as an enriched broken space so that the
Riesz operator becomes block-diagonal and can be inverted locally on each mesh
element. Trace unknowns on the mesh skeleton are then introduced to recover the
conformity of the solution across element interfaces. This is the key idea behind
ultra-weak DPG formulations and is the setting used for time-harmonic Maxwell
problems.

Ultra-weak formulation for the Maxwell problem
----------------------------------------------

The starting point is the time-harmonic Maxwell system written in the form

.. math::

   \begin{aligned}
   \nabla \times \mathbf{E} - i\omega \mu \mathbf{H} &= 0, \\
   \nabla \times \mathbf{H} + i\omega \varepsilon_{eff} \mathbf{E} &= \mathbf{J},
   \end{aligned}

where :math:`\mathbf{E}` is the electric field, :math:`\mathbf{H}` the magnetic
field, :math:`\mu` the permeability, :math:`\varepsilon_{eff}` the effective
permittivity, and :math:`\mathbf{J}` a prescribed current density. After
non-dimensionalization, one writes

.. math::

   \begin{aligned}
   \nabla \times \mathbf{E} - i\tilde{\omega} \mu_r \mathbf{H} &= 0, \\
   \nabla \times \mathbf{H} + i\tilde{\omega} (\varepsilon_r + i\sigma_r) \mathbf{E}
   &= \mathbf{\tilde{J}}.
   \end{aligned}

In the DPG framework, these equations are weakened by transferring derivatives to
the test functions and by introducing trace unknowns on the mesh skeleton. The
resulting formulation is naturally suited to broken finite element spaces and to
boundary conditions such as port, absorbing, or impedance conditions. This is
especially useful for waveguide and scattering problems, where the boundary
operators play an important role in the stability of the discrete system.

Waveguide setting and boundary conditions
-----------------------------------------

A typical application is the propagation of a guided mode in a rectangular
waveguide. The geometry is assumed to be a box of dimensions :math:`a` and
:math:`b` in the :math:`x` and :math:`y` directions, respectively, and the
fields are excited by a transverse electric (TE) mode. The TE\ :math:`_{mn}`
mode is characterized by

.. math::

   k_x = \frac{m\pi}{a}, \qquad k_y = \frac{n\pi}{b},

and the cutoff wavenumber is

.. math::

   k_c = \sqrt{k_x^2 + k_y^2}.

The longitudinal wavenumber is then

.. math::

   k_z = \sqrt{\omega^2 \mu \varepsilon_{eff} - k_x^2 - k_y^2}.

The boundary conditions used in the waveguide problem are of three different
types. On :math:`\Gamma_1` one imposes a perfect electric conductor (PEC)
condition,

.. math::

   \mathbf{n} \times \mathbf{E} = 0,

on :math:`\Gamma_2` an incident field is prescribed through an impedance-type
port condition,

.. math::

   \mathbf{n} \times \mathbf{H} + \frac{k_z}{\omega \mu_r} \mathbf{E}_{\parallel} = g,

and on :math:`\Gamma_3` an absorbing condition is used with :math:`g = 0`.
These boundary operators are naturally handled by the ultra-weak DPG formulation
because the boundary terms are treated explicitly through the traces.

For completeness, the incident TE fields can be written in compact form as

.. math::

   \mathbf{E}_{inc} = i\frac{\omega \mu_r}{k_c^2}
   \begin{bmatrix}
   -k_y \cos(k_x x) \sin(k_y y) \\
   k_x \sin(k_x x) \cos(k_y y) \\
   0
   \end{bmatrix} e^{i k_z z},

and

.. math::

   \mathbf{H}_{inc} =
   \begin{bmatrix}
   - i \dfrac{k_z k_x}{k_c^2} \sin(k_x x) \cos(k_y y) \\
   - i \dfrac{k_z k_y}{k_c^2} \cos(k_x x) \sin(k_y y) \\
   \cos(k_x x) \cos(k_y y)
   \end{bmatrix} e^{i k_z z}.

These expressions are useful for defining the incoming wave on the port and for
comparing numerical solutions with the expected guided-mode behavior.
