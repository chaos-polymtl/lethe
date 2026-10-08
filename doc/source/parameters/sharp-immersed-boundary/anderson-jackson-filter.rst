..
  SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#######################
Anderson-Jackson Filter
#######################

The ``lethe-fluid-sharp-filter`` application computes the phase-averaged fields of Anderson and Jackson [#anderson1967]_ from a resolved CFD-DEM snapshot of ``lethe-fluid-sharp``. These fields are the fluid and solid volume fractions and the fluid velocity and pressure averaged over a filter of finite width. They are the fields that an unresolved (volume-averaged) CFD-DEM model describes, and can thus be used to develop and verify such models.

The theory of the filter and its parallel algorithm are described in :doc:`../../theory/multiphase/cfd_dem/filtering_resolved_cfd-dem`.

The application reads the same parameter file as ``lethe-fluid-sharp``, extended with the ``anderson jackson filter`` subsection. It restores the mesh, the velocity-pressure solution and the particles from the checkpoint of the simulation when ``restart = true`` in the ``restart`` subsection (see :doc:`../cfd/restart`), filters them, and writes the filtered fields. It does not advance the simulation in time. When ``restart = false``, the initial condition of the simulation is filtered instead.

.. code-block:: text

  mpirun -np 8 lethe-fluid-sharp-filter snapshot.prm

***************
Filtered fields
***************

For a filter centered at :math:`\mathbf{x}`, the application computes

.. math::

  M(\mathbf{x}) = \int_\Omega g(|\mathbf{x}-\mathbf{y}|) \, \mathrm{d}V_\mathbf{y}, \qquad
  E(\mathbf{x}) = \int_\Omega I_f(\mathbf{y}) \, g(|\mathbf{x}-\mathbf{y}|) \, \mathrm{d}V_\mathbf{y},

.. math::

  \bar{\mathbf{u}}_f(\mathbf{x}) = \frac{1}{E(\mathbf{x})} \int_\Omega I_f \, \mathbf{u} \, g \, \mathrm{d}V, \qquad
  \bar{p}_f(\mathbf{x}) = \frac{1}{E(\mathbf{x})} \int_\Omega I_f \, p \, g \, \mathrm{d}V,

where :math:`g` is a kernel of unit mass with compact support and :math:`I_f` is the fluid indicator, which is one in the fluid and zero in the particles. The fluid volume fraction is :math:`\epsilon_f = E` and the solid volume fraction is :math:`\epsilon_s = M - E`. Where the support of the kernel does not reach a non-periodic boundary of the domain, the kernel mass :math:`M` is exactly one, and both volume fractions are divided by its computed value, which removes the quadrature error of the kernel and ensures that :math:`\epsilon_s = 1 - \epsilon_f`. Near a wall, the kernel is truncated by the boundary and :math:`M < 1`: the ``normalize at domain boundaries`` and ``extend velocity beyond walls`` parameters select how the part of the kernel beyond the boundary is treated.

Three kernels are available:

* The ``gaussian`` kernel of standard deviation :math:`\sigma` is truncated at the radius :math:`R = c\sigma`, where :math:`c` is the ``gaussian cutoff``, and renormalized so that its integral over the ball of radius :math:`R` is exactly one:

  .. math::

    g(r) = \frac{\exp\left(-r^2/2\sigma^2\right)}{(2\pi\sigma^2)^{d/2} \, m(c)} \quad \text{if } r < R, \qquad g(r) = 0 \quad \text{otherwise},

  where :math:`m(c)` is the mass of the untruncated gaussian contained in the ball of radius :math:`R`: :math:`m(c) = 1 - e^{-c^2/2}` in 2D and :math:`m(c) = \mathrm{erf}(c/\sqrt{2}) - \sqrt{2/\pi} \, c \, e^{-c^2/2}` in 3D. Without the renormalization, the kernel mass would only be 0.971 in 3D for a cutoff of 3, which would bias every volume fraction by the same amount.

* The ``top-hat`` kernel is the normalized indicator function of the ball of radius :math:`R`: :math:`g(r) = 1/|B_R|` if :math:`r < R`.

* The ``wendland`` kernel is the Wendland C2 polynomial of support radius :math:`R`:

  .. math::

    g(r) = \alpha_d \left(1 - \frac{r}{R}\right)^4 \left(1 + 4 \frac{r}{R}\right) \quad \text{if } r < R, \qquad g(r) = 0 \quad \text{otherwise},

  with :math:`\alpha_2 = 7/(\pi R^2)` and :math:`\alpha_3 = 21/(2\pi R^3)`. It vanishes smoothly at its support radius, along with its first two derivatives, which makes its quadrature much more accurate than that of the two other kernels and makes the gradients of the filtered fields continuous. It has the second moment of a gaussian of standard deviation :math:`0.264 R` in 2D and :math:`0.258 R` in 3D.

The filtered fields are computed at the support points (the filter centers) of a continuous ``FE_Q`` space of arbitrary degree defined on the mesh of the simulation, and are thus finite element fields on this mesh. The values at the hanging nodes are interpolated from the constraints.

**********
Parameters
**********

.. code-block:: text

  subsection anderson jackson filter
    set kernel type                    = gaussian
    set filter width                   = 1
    set gaussian cutoff                = 3
    set output polynomial degree       = 0
    set quadrature points              = 0
    set cut cell subdivisions          = 4
    set filter pressure                = true
    set minimum fluid fraction         = 1e-12
    set normalize at domain boundaries = false
    set extend velocity beyond walls   = true
    set output folder                  = ./filter_output/
    set output name                    = filtered
    set verbosity                      = quiet
  end

* The ``kernel type`` parameter selects the kernel of the filter. Choices are ``gaussian``, ``top-hat`` and ``wendland``.

* The ``filter width`` parameter is the standard deviation :math:`\sigma` of the gaussian kernel, or the support radius :math:`R` of the top-hat and Wendland kernels. It must be strictly positive.

* The ``gaussian cutoff`` parameter is the truncation radius of the gaussian kernel, expressed in number of standard deviations. The kernel is discontinuous at the truncation radius, where its value relative to its maximum is :math:`e^{-c^2/2}`, which increases the quadrature error of the kernel mass. A cutoff between 3 and 4 is recommended: with the default number of quadrature points, increasing it from 3 to 4 reduces the error of the kernel mass by about two orders of magnitude, for a support whose volume is 2.4 times larger in 3D. This parameter is only used by the gaussian kernel.

* The ``output polynomial degree`` parameter is the degree of the ``FE_Q`` space of the filtered fields. If it is 0, the velocity degree of the simulation is used.

* The ``quadrature points`` parameter is the number of Gauss points per direction used to integrate the cells of the mesh. If it is 0, the velocity degree plus one is used. Increasing it reduces the quadrature error of the kernel, at a cost that grows with its power ``dim``. The reduction is fast for the Wendland kernel, but slow for the top-hat kernel and limited by the cutoff for the gaussian kernel, since these two kernels are discontinuous at their support radius.

* The ``cut cell subdivisions`` parameter is the number of subdivisions per direction of the Gauss quadrature used on the cells cut by a particle. The fluid indicator is discontinuous in these cells and is evaluated at every point of the subdivided quadrature. The error on the volume of the particles decreases with the number of subdivisions, but the cost of a cut cell grows with its power ``dim``.

* The ``filter pressure`` parameter enables the filtering of the pressure. The velocity is always filtered.

* The ``minimum fluid fraction`` parameter is the ratio :math:`E/M` below which the phase-averaged velocity and pressure are considered undefined, for example at a filter center whose kernel lies entirely inside a particle. The undefined averages are set to zero and are flagged by the ``valid`` field.

* The ``normalize at domain boundaries`` parameter divides the fluid and solid volume fractions by the kernel mass :math:`M` where the kernel is truncated by a non-periodic boundary, which renormalizes the kernel. With this option, :math:`\epsilon_f + \epsilon_s = 1` everywhere, but the filter is no longer a convolution near the boundaries. The phase-averaged velocity and pressure are not affected, since they are ratios of integrals of the same kernel. Where the kernel does not reach a non-periodic boundary, or where its computed mass exceeds one, the volume fractions are always divided by the kernel mass, whatever the value of this parameter, since this only removes the quadrature error of the kernel.

* The ``extend velocity beyond walls`` parameter fills the part of the kernel located outside of the domain, beyond the walls where the velocity is imposed (``noslip``, ``function`` and ``function weak`` boundary conditions), with fluid moving at the velocity of the wall. This fluid is counted both in the fluid volume fraction, which becomes :math:`\epsilon_f = E + M_e` where :math:`M_e` is the mass of the kernel beyond the walls, and in the phase-averaged velocity, so that the product :math:`\epsilon_f \bar{\mathbf{u}}_f` remains the flux of fluid seen by the kernel. Near a domain bounded by walls only, :math:`\epsilon_f + \epsilon_s = 1`. The velocity of the wall is evaluated at the time of the snapshot. Without this extension, the kernel is truncated by the walls, and the averaged velocity near a wall is the velocity at the centroid of the truncated kernel, which lies inside the domain. For a linear velocity profile, the extension halves the deviation of the averaged velocity from the exact one at the wall, and raises the gradient of the averaged velocity on the wall from about a third to one half of the exact gradient (see :doc:`../../theory/multiphase/cfd_dem/filtering_resolved_cfd-dem`). The solid volume fraction and the averaged pressure are not affected. The other boundaries, such as outlets and slip walls, always truncate the kernel. When ``normalize at domain boundaries`` is also enabled, the volume fractions are divided by :math:`M + M_e`, which only differs from one near these other boundaries.

* The ``output folder`` and ``output name`` parameters define the location and the prefix of the vtu files and of the pvd file of the filtered fields. They must differ from the ``output folder`` and ``output name`` of the ``simulation control`` subsection, so that the output of the simulation is not overwritten. The number of vtu files is set by the ``group files`` parameter of the ``simulation control`` subsection.

* The ``verbosity`` parameter controls the diagnostics printed by the filter. With ``verbose``, the parameters of the kernel and global statistics of the filtered fields are printed. With ``extra verbose``, the distribution of the filter centers among the processes is printed as well. The choices are ``quiet``, ``verbose`` and ``extra verbose``.

The time spent in each stage of the filter is printed when the ``type`` of the ``timer`` subsection is not ``none``.

*******
Outputs
*******

The following fields are written:

* ``fluid_volume_fraction`` and ``solid_volume_fraction``: the volume fractions :math:`\epsilon_f` and :math:`\epsilon_s`. They sum to one wherever the kernel does not reach a non-periodic boundary;
* ``kernel_mass``: the kernel mass :math:`M`, which is one away from the non-periodic boundaries, up to the quadrature error of the kernel. Its deviation from one away from these boundaries measures the accuracy of the quadrature of the kernel;
* ``filtered_velocity``: the phase-averaged fluid velocity :math:`\bar{\mathbf{u}}_f`;
* ``filtered_velocity_gradient``: the gradient :math:`\nabla \bar{\mathbf{u}}_f` of the phase-averaged fluid velocity, computed from its finite element interpolant. Its components are numbered as those of the ``velocity_gradient`` of the fluid solvers (0 = xx, 1 = xy, ...), where xy is the derivative of the x component along y. This is the gradient of the averaged velocity, which differs from the phase average of the velocity gradient by the integral of the velocity over the surface of the particles within the kernel. Near the filter centers where the averages are not defined (``valid`` = 0), the averaged velocity drops to zero and its gradient is meaningless. Within the support radius of a wall, the kernel is truncated and off-centered, which underestimates this gradient, down to about a third of the exact gradient on the wall for a gaussian truncated at three standard deviations, or one half with ``extend velocity beyond walls = true`` (see :doc:`../../theory/multiphase/cfd_dem/filtering_resolved_cfd-dem`);
* ``filtered_pressure``: the phase-averaged fluid pressure :math:`\bar{p}_f`, if the pressure is filtered;
* ``valid``: one where the phase-averaged velocity and pressure are defined, zero otherwise.

With ``verbosity = verbose``, the filter also prints the integral of the solid volume fraction over the domain. For a periodic domain, this integral equals the volume of the particles up to discretization errors, since the kernel has unit mass. Comparing it with the solid volume of the source cells, which is also printed, verifies that every part of the domain contributed to every filter center within its support.

******************
Algorithm and cost
******************

The convolution is computed without assembling a filtering matrix and without evaluating the solution at remote points. Every process integrates the cells it owns once and adds the contribution of each batch of quadrature points to the filter centers within the support radius of the kernel, which it finds with an R-tree. Beforehand, every filter center is sent to the processes that own cells within the support radius of the kernel, and the partial integrals are returned to the owner of the filter center afterwards.

The cost of the filter is dominated by the evaluations of the kernel. It scales as the number of filter centers times the number of quadrature points within the support of the kernel, i.e. roughly as :math:`N_c \, n_q^d \, \frac{4}{3}\pi (R/h)^3` in 3D, where :math:`N_c` is the number of filter centers, :math:`n_q` is the number of quadrature points per direction and :math:`h` is the size of the cells. The cost therefore grows as the cube of the ratio between the support radius and the cell size. For example, doubling the filter width multiplies the cost by eight. The memory used by the halo of filter centers received from the other processes also grows with the support radius. It can be monitored with ``verbosity = extra verbose``.

***********
Limitations
***********

* Periodic boundaries are supported by filtering the periodic images of the filter centers. The support radius of the kernel may exceed the period of the domain: the kernel is then periodized, and its cost grows with the number of periodic images within its support. The periodic boundaries must be planes normal to their direction of periodicity. Like ``lethe-fluid-sharp``, the filter does not consider the periodic images of the particles themselves.
* The cells are classified as fluid, solid or cut with the signed distance functions of the particles, assuming that these functions are 1-Lipschitz, which holds for the exact signed distances of the analytical shapes. Shapes defined from IGES files are not supported.
* Only quadrilateral and hexahedral meshes are supported.

*********
Reference
*********

.. [#anderson1967] \T. B. Anderson and R. Jackson, "Fluid mechanical description of fluidized beds. Equations of motion," *Industrial & Engineering Chemistry Fundamentals*, vol. 6, no. 4, pp. 527–539, 1967, doi: `10.1021/i160024a007 <https://doi.org/10.1021/i160024a007>`_\.
