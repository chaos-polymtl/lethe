..
  SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#######################
Anderson-Jackson Filter
#######################

The ``lethe-fluid-sharp-filter`` application computes the phase-averaged fields of Anderson and Jackson [#anderson1967]_ from a resolved CFD-DEM snapshot of ``lethe-fluid-sharp``. These fields are the fluid and solid volume fractions and the fluid velocity and pressure averaged over a filter of finite width. They are the fields that an unresolved (volume-averaged) CFD-DEM model describes, and can thus be used to develop and verify such models.

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

where :math:`g` is a kernel of unit mass with compact support and :math:`I_f` is the fluid indicator, which is one in the fluid and zero in the particles. The fluid volume fraction is :math:`\epsilon_f = E` and the solid volume fraction is :math:`\epsilon_s = M - E`. Away from the boundaries of the domain, the kernel mass :math:`M` is one and :math:`\epsilon_s = 1 - \epsilon_f`. Near a wall, the kernel is truncated by the boundary and :math:`M < 1`.

Two kernels are available:

* The ``gaussian`` kernel of standard deviation :math:`\sigma` is truncated at the radius :math:`R = c\sigma`, where :math:`c` is the ``gaussian cutoff``, and renormalized so that its integral over the ball of radius :math:`R` is exactly one:

  .. math::

    g(r) = \frac{\exp\left(-r^2/2\sigma^2\right)}{(2\pi\sigma^2)^{d/2} \, m(c)} \quad \text{if } r < R, \qquad g(r) = 0 \quad \text{otherwise},

  where :math:`m(c)` is the mass of the untruncated gaussian contained in the ball of radius :math:`R`: :math:`m(c) = 1 - e^{-c^2/2}` in 2D and :math:`m(c) = \mathrm{erf}(c/\sqrt{2}) - \sqrt{2/\pi} \, c \, e^{-c^2/2}` in 3D. Without the renormalization, the kernel mass would only be 0.971 in 3D for a cutoff of 3, which would bias every volume fraction by the same amount.

* The ``top-hat`` kernel is the normalized indicator function of the ball of radius :math:`R`: :math:`g(r) = 1/|B_R|` if :math:`r < R`.

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
    set output folder                  = ./filter_output/
    set output name                    = filtered
    set verbosity                      = quiet
  end

* The ``kernel type`` parameter selects the kernel of the filter. Choices are ``gaussian`` and ``top-hat``.

* The ``filter width`` parameter is the standard deviation :math:`\sigma` of the gaussian kernel, or the radius :math:`R` of the top-hat kernel. It must be strictly positive.

* The ``gaussian cutoff`` parameter is the truncation radius of the gaussian kernel, expressed in number of standard deviations. The kernel is discontinuous at the truncation radius, where its value relative to its maximum is :math:`e^{-c^2/2}`, which increases the quadrature error of the kernel mass. A cutoff between 3 and 4 is recommended. This parameter is not used by the top-hat kernel.

* The ``output polynomial degree`` parameter is the degree of the ``FE_Q`` space of the filtered fields. If it is 0, the velocity degree of the simulation is used.

* The ``quadrature points`` parameter is the number of Gauss points per direction used to integrate the cells of the mesh. If it is 0, the velocity degree plus one is used.

* The ``cut cell subdivisions`` parameter is the number of subdivisions per direction of the Gauss quadrature used on the cells cut by a particle. The fluid indicator is discontinuous in these cells and is evaluated at every point of the subdivided quadrature. The error on the volume of the particles decreases with the number of subdivisions, but the cost of a cut cell grows with its power ``dim``.

* The ``filter pressure`` parameter enables the filtering of the pressure. The velocity is always filtered.

* The ``minimum fluid fraction`` parameter is the ratio :math:`E/M` below which the phase-averaged velocity and pressure are considered undefined, for example at a filter center whose kernel lies entirely inside a particle. The undefined averages are set to zero and are flagged by the ``valid`` field.

* The ``normalize at domain boundaries`` parameter divides the fluid and solid volume fractions by the kernel mass :math:`M`, which renormalizes the kernel where it is truncated by a wall. With this option, :math:`\epsilon_f + \epsilon_s = 1` everywhere, but the filter is no longer a convolution near the walls. The phase-averaged velocity and pressure are not affected, since they are ratios of integrals of the same kernel.

* The ``output folder`` and ``output name`` parameters define the location and the prefix of the vtu files and of the pvd file of the filtered fields. They must differ from the ``output folder`` and ``output name`` of the ``simulation control`` subsection, so that the output of the simulation is not overwritten. The number of vtu files is set by the ``group files`` parameter of the ``simulation control`` subsection.

* The ``verbosity`` parameter controls the diagnostics printed by the filter. With ``verbose``, the parameters of the kernel and global statistics of the filtered fields are printed. With ``extra verbose``, the distribution of the filter centers among the processes is printed as well. The choices are ``quiet``, ``verbose`` and ``extra verbose``.

The time spent in each stage of the filter is printed when the ``type`` of the ``timer`` subsection is not ``none``.

*******
Outputs
*******

The following fields are written:

* ``fluid_volume_fraction`` and ``solid_volume_fraction``: the volume fractions :math:`\epsilon_f` and :math:`\epsilon_s`;
* ``kernel_mass``: the kernel mass :math:`M`, which is one away from the walls, up to the quadrature error of the kernel;
* ``filtered_velocity``: the phase-averaged fluid velocity :math:`\bar{\mathbf{u}}_f`;
* ``filtered_pressure``: the phase-averaged fluid pressure :math:`\bar{p}_f`, if the pressure is filtered;
* ``valid``: one where the phase-averaged velocity and pressure are defined, zero otherwise.

With ``verbosity = verbose``, the filter also prints the integral of the solid volume fraction over the domain. For a periodic domain without renormalization, this integral equals the volume of the particles, since the kernel has unit mass. Comparing it with the solid volume of the source cells, which is also printed, verifies that every part of the domain contributed to every filter center within its support.

******************
Algorithm and cost
******************

The convolution is computed without assembling a filtering matrix and without evaluating the solution at remote points. Every process integrates the cells it owns once and adds the contribution of each batch of quadrature points to the filter centers within the support radius of the kernel, which it finds with an R-tree. Beforehand, every filter center is sent to the processes that own cells within the support radius of the kernel, and the partial integrals are returned to the owner of the filter center afterwards.

The cost of the filter is dominated by the evaluations of the kernel. It scales as the number of filter centers times the number of quadrature points within the support of the kernel, i.e. roughly as :math:`N_c \, n_q^d \, \frac{4}{3}\pi (R/h)^3` in 3D, where :math:`N_c` is the number of filter centers, :math:`n_q` is the number of quadrature points per direction and :math:`h` is the size of the cells. The cost therefore grows as the cube of the ratio between the support radius and the cell size. For example, doubling the filter width multiplies the cost by eight. The memory used by the halo of filter centers received from the other processes also grows with the support radius. It can be monitored with ``verbosity = extra verbose``.

***********
Limitations
***********

* Periodic boundaries are supported by filtering the periodic images of the filter centers. The support radius of the kernel must not exceed half of the period of the domain in each periodic direction. Like ``lethe-fluid-sharp``, the filter does not consider the periodic images of the particles themselves.
* The cells are classified as fluid, solid or cut with the signed distance functions of the particles, assuming that these functions are 1-Lipschitz, which holds for the exact signed distances of the analytical shapes. Shapes defined from IGES files are not supported.
* Only quadrilateral and hexahedral meshes are supported.

*********
Reference
*********

.. [#anderson1967] \T. B. Anderson and R. Jackson, "Fluid mechanical description of fluidized beds. Equations of motion," *Industrial & Engineering Chemistry Fundamentals*, vol. 6, no. 4, pp. 527–539, 1967, doi: `10.1021/i160024a007 <https://doi.org/10.1021/i160024a007>`_\.
