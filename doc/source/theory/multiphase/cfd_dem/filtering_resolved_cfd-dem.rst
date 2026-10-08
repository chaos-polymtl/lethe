..
  SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

==================================
Filtering Resolved CFD-DEM Results
==================================

A resolved CFD-DEM simulation (see :doc:`resolved_cfd-dem`) describes the flow between the particles at the scale of the pores. An unresolved CFD-DEM model (see :doc:`unresolved_cfd-dem`) describes the same flow at a coarser scale, through fields averaged over a region containing several particles: the void fraction and the phase-averaged velocity and pressure of the fluid. The closures of the unresolved models, such as the drag models, are commonly derived from resolved simulations. For such a comparison to be meaningful, the resolved fields must be averaged with the same operator as the one that defines the Volume-Averaged Navier-Stokes (VANS) equations.

The ``lethe-fluid-sharp-filter`` application computes these averaged fields from a snapshot of ``lethe-fluid-sharp``. This section presents the averaging operator, its discretization, and the parallel algorithm that makes it affordable for large simulations. The parameters of the application are described in :doc:`../../../parameters/sharp-immersed-boundary/anderson-jackson-filter`.

Anderson-Jackson averaging
--------------------------

Following Anderson and Jackson [#anderson1967]_ [#jackson2000]_, the local average of a quantity is defined by a convolution with a weighting function, or kernel, :math:`g`. The kernel is radially symmetric, has a unit mass and, in Lethe, a compact support of radius :math:`R`:

.. math::
    \int_{\mathbb{R}^d} g(\lVert \mathbf{x} \rVert) \, \mathrm{d}\mathbf{x} = 1, \qquad g(r) = 0 \quad \text{for } r \geq R.

Let :math:`I_f` be the fluid indicator, which is one in the fluid and zero in the particles. At a point :math:`\mathbf{x}`, called a filter center, the fluid volume fraction (or void fraction) and the phase average of a fluid quantity :math:`a` are

.. math::
    \varepsilon_f(\mathbf{x}) &= \int_\Omega I_f(\mathbf{y}) \, g(\lVert \mathbf{x}-\mathbf{y} \rVert) \, \mathrm{d}\mathbf{y} \\
    \varepsilon_f(\mathbf{x}) \, \langle a \rangle_f(\mathbf{x}) &= \int_\Omega I_f(\mathbf{y}) \, a(\mathbf{y}) \, g(\lVert \mathbf{x}-\mathbf{y} \rVert) \, \mathrm{d}\mathbf{y}

These are the definitions used in the unresolved CFD-DEM section, where the kernel is denoted :math:`k_r`. The application computes the phase-averaged fluid velocity :math:`\bar{\mathbf{u}}_f = \langle \mathbf{u} \rangle_f` and pressure :math:`\bar{p}_f = \langle p \rangle_f`, as well as the kernel mass and the solid volume fraction:

.. math::
    M(\mathbf{x}) = \int_\Omega g(\lVert \mathbf{x}-\mathbf{y} \rVert) \, \mathrm{d}\mathbf{y}, \qquad \varepsilon_s(\mathbf{x}) = M(\mathbf{x}) - \varepsilon_f(\mathbf{x}).

Where the support of the kernel does not reach a non-periodic boundary of the domain, :math:`M = 1` and :math:`\varepsilon_s = 1 - \varepsilon_f`. The discrete kernel mass however differs from one by the quadrature error of the kernel (see below). Since :math:`M` is known exactly there, the application divides both volume fractions by the discrete kernel mass, :math:`\varepsilon_f = E / M` and :math:`\varepsilon_s = (M - E) / M`. This does not change their definition, removes the quadrature error of the kernel mass from the volume fractions, and ensures that they sum to one, in the same way as the phase averages, which are ratios of two integrals of the same kernel. The same division is applied wherever the discrete mass of the kernel exceeds one, since the mass of the kernel is at most one and any excess is thus a quadrature error. The volume fractions therefore never exceed one.

Near a wall, the support of the kernel leaves the domain and :math:`M < 1`: the averages of a truncated kernel are still well defined, since :math:`\langle a \rangle_f` is a ratio of two integrals of the same kernel, but the volume fractions no longer sum to one. Three treatments of the part of the kernel beyond the boundaries are available:

* the kernel is truncated, and the volume fractions :math:`\varepsilon_f = E` and :math:`\varepsilon_s = M - E` sum to :math:`M`;
* the volume fractions are divided by :math:`M`, which renormalizes the kernel near the boundaries;
* the part of the kernel beyond the walls is filled with fluid moving at the velocity of the wall (see below), which is the default.

:math:`\langle a \rangle_f` is undefined where the kernel contains no fluid, for example at a filter center deep inside a large particle. The application reports these filter centers and sets their averages to zero.

Integrating the solid volume fraction over a periodic domain gives

.. math::
    \int_\Omega \varepsilon_s(\mathbf{x}) \, \mathrm{d}\mathbf{x} = \int_\Omega \left(1 - I_f(\mathbf{y})\right) \int_\Omega g(\lVert \mathbf{x}-\mathbf{y} \rVert) \, \mathrm{d}\mathbf{x} \, \mathrm{d}\mathbf{y} = V_s,

where :math:`V_s` is the volume of the particles, since the kernel has unit mass. The application prints both quantities, which provides a global check of the filtering.

Gradients of averaged quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Differentiating the definition of the phase average with respect to :math:`\mathbf{x}` and integrating by parts yields, when the support of the kernel lies inside the domain (or the domain is periodic), the averaging theorem [#jackson2000]_:

.. math::
    \nabla \left( \varepsilon_f \langle a \rangle_f \right) = \varepsilon_f \langle \nabla a \rangle_f + \int_{S_p} a(\mathbf{y}) \, \mathbf{n}(\mathbf{y}) \, g(\lVert \mathbf{x}-\mathbf{y} \rVert) \, \mathrm{d}S_\mathbf{y},

where :math:`S_p` is the surface of the particles and :math:`\mathbf{n}` its unit normal pointing into the fluid. The gradient of an averaged quantity is thus not the average of its gradient: the two differ by the surface integral, which is the origin of the interphase momentum exchange terms of the VANS equations. The velocity gradient written by the application is the gradient of the averaged velocity, :math:`\nabla \bar{\mathbf{u}}_f`, obtained by differentiating the finite element representation of :math:`\bar{\mathbf{u}}_f`.

Near a wall, the integration by parts also produces a term on the part of the wall within the support of the kernel, and the wall plays the same role as the surface of a particle. Geometrically, the kernel truncated by the wall is no longer centered on the filter center: its centroid lies at a distance :math:`\delta` inside the domain, which decreases from about :math:`0.8\sigma` on the wall to zero at the distance :math:`R` from the wall. For a linear velocity profile :math:`u(y) = \dot{\gamma} y` along the normal to the wall, the averaged velocity is the velocity at the centroid of the kernel, :math:`\bar{u}_f(y) = \dot{\gamma} \left(y + \delta(y)\right)`, and its gradient is

.. math::
    \frac{\mathrm{d}\bar{u}_f}{\mathrm{d}y} = \dot{\gamma} \left( 1 + \frac{\mathrm{d}\delta}{\mathrm{d}y} \right),

which is smaller than :math:`\dot{\gamma}` within the distance :math:`R` of the wall, since :math:`\delta` decreases. For a 2D gaussian kernel truncated at :math:`3\sigma`, the gradient of the averaged velocity on the wall is only 37% of the exact gradient, and it recovers 89% of it at :math:`2\sigma` and 100% at :math:`3\sigma`, while the averaged velocity itself never differs from the exact one by more than :math:`0.8\sigma\dot{\gamma}`. The phase average of the gradient, :math:`\langle \nabla \mathbf{u} \rangle_f`, is not affected: it is exactly :math:`\dot{\gamma}` everywhere for this profile. Renormalizing the kernel by its mass does not change this behavior, since the phase averages are already ratios of integrals of the same kernel.

Extension of the fluid beyond the walls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Optionally, the part of the kernel outside of the domain, beyond a wall where the velocity :math:`\mathbf{u}_w` is imposed (Dirichlet boundary condition), is considered as fluid moving at the velocity of the wall. This exterior fluid is counted both in the fluid volume fraction and in the phase-averaged velocity:

.. math::
    \varepsilon_f(\mathbf{x}) = E(\mathbf{x}) + M_e(\mathbf{x}), \qquad
    \bar{\mathbf{u}}_f(\mathbf{x}) = \frac{\int_\Omega I_f \, \mathbf{u} \, g \, \mathrm{d}V + M_e(\mathbf{x}) \, \bar{\mathbf{u}}_w(\mathbf{x})}{E(\mathbf{x}) + M_e(\mathbf{x})},

where :math:`E = \int_\Omega I_f \, g \, \mathrm{d}V`, :math:`M_e` is the mass of the kernel beyond the walls and :math:`\bar{\mathbf{u}}_w` is the velocity of the walls within the support of the kernel. The product :math:`\varepsilon_f \bar{\mathbf{u}}_f` thus remains the flux of fluid seen by the kernel, :math:`\int_\Omega I_f \, \mathbf{u} \, g \, \mathrm{d}V + M_e \bar{\mathbf{u}}_w`, which reduces to the flux of the truncated kernel near a stationary wall. The latter is divergence-free for fixed particles and impermeable walls, a property that the averaged fields would lose if the exterior fluid was only counted in the velocity. Near a domain bounded by walls only, :math:`\varepsilon_f + \varepsilon_s = M + M_e = 1`. Since the kernel has unit mass, the mass of the kernel outside of the domain is :math:`1 - M`. When the support of the kernel also crosses boundaries without imposed velocity (e.g. outlets), this mass is attributed to the walls in proportion to the integral of the kernel over the walls, :math:`S_w`, and over all the non-periodic boundaries, :math:`S`:

.. math::
    M_e = (1 - M) \frac{S_w}{S}, \qquad
    S_w(\mathbf{x}) = \int_{\Gamma_w} g(\lVert \mathbf{x}-\mathbf{y} \rVert) \, \mathrm{d}S_\mathbf{y}, \qquad
    \bar{\mathbf{u}}_w(\mathbf{x}) = \frac{1}{S_w(\mathbf{x})} \int_{\Gamma_w} \mathbf{u}_w(\mathbf{y}) \, g(\lVert \mathbf{x}-\mathbf{y} \rVert) \, \mathrm{d}S_\mathbf{y}.

This attribution is exact near a single wall, or near several walls moving at the same velocity. The surface integrals are accumulated by the same source-centric sweep as the volume integrals, from the quadrature points of the boundary faces of the locally owned cells, and they are communicated with the other moments. The integral :math:`S` is always accumulated, since it also identifies the filter centers whose kernel reaches a non-periodic boundary (:math:`S > 0`), where the volume fractions are not divided by the kernel mass unless the renormalization is requested. With both the extension and the renormalization, the volume fractions are divided by :math:`M + M_e`, which only differs from one near the boundaries that are not walls. The solid volume fraction and the averaged pressure are not affected by the extension.

For the linear profile above, the extended velocity has a kink at the wall, where its gradient jumps from :math:`\dot{\gamma}` to zero. The gradient of the averaged velocity is then the average of this gradient, i.e. :math:`\dot{\gamma}` times the fraction of the kernel inside the domain, which is 1/2 on the wall. The extension thus raises the gradient of the averaged velocity on the wall from 37% to 50% of the exact gradient for a 2D gaussian truncated at :math:`3\sigma`, it recovers 84% of it at :math:`\sigma` and 98% at :math:`2\sigma`, and it halves the deviation of the averaged velocity from the exact one on the wall, down to :math:`0.4\sigma\dot{\gamma}`. In any case, the gradient of the averaged velocity must be interpreted with care within the distance :math:`R` of the walls.

Kernels
~~~~~~~

Three kernels are available. The top-hat kernel is the normalized indicator function of the ball :math:`B_R` of radius :math:`R`, :math:`g(r) = 1/|B_R|` for :math:`r < R`. The gaussian kernel of standard deviation :math:`\sigma` is truncated at :math:`R = c\sigma` and renormalized:

.. math::
    g(r) = \frac{\exp\left(-r^2 / 2\sigma^2\right)}{(2\pi\sigma^2)^{d/2} \, m_d(c)} \quad \text{for } r < R, \qquad
    m_2(c) = 1 - e^{-c^2/2}, \qquad
    m_3(c) = \mathrm{erf}\left(\frac{c}{\sqrt{2}}\right) - \sqrt{\frac{2}{\pi}} \, c \, e^{-c^2/2},

where :math:`m_d(c)` is the mass of the untruncated gaussian contained in :math:`B_R`. The renormalization matters: in 3D, a gaussian truncated at three standard deviations retains only 97.1% of its mass, and would bias every volume fraction by 2.9% without it. The truncated kernel is discontinuous at :math:`r = R`, where its value relative to its maximum is :math:`e^{-c^2/2}`, which limits the accuracy of its quadrature. A cutoff :math:`c` between 3 and 4 is a good compromise between this discontinuity and the size of the support, which drives the cost of the filter.

The Wendland C2 kernel [#wendland1995]_ of support radius :math:`R` is the polynomial

.. math::
    g(r) = \alpha_d \left(1 - \frac{r}{R}\right)^4 \left(1 + 4 \frac{r}{R}\right) \quad \text{for } r < R, \qquad
    \alpha_2 = \frac{7}{\pi R^2}, \qquad
    \alpha_3 = \frac{21}{2 \pi R^3}.

Unlike the two other kernels, it vanishes at :math:`r = R` along with its first two derivatives, so it is twice continuously differentiable. The averaging theorem above, which is obtained by integrating by parts, assumes such a differentiable kernel, and the gradients of the fields filtered with this kernel are continuous. Its shape is close to a gaussian: it has the second moment of a gaussian of standard deviation :math:`\sqrt{5/72} \, R \approx 0.264 R` in 2D and :math:`R/\sqrt{15} \approx 0.258 R` in 3D, i.e. it is comparable to a gaussian truncated at about 3.8 standard deviations.

The smoothness of the kernel at its support radius controls the accuracy of its quadrature. The error of the kernel mass of the top-hat kernel, which jumps from its maximum to zero, only decreases slowly and irregularly with the number of quadrature points. The error of the truncated gaussian decreases until it reaches the level set by its jump :math:`e^{-c^2/2}`, so that increasing the cutoff from 3 to 4 is more effective than adding quadrature points. The error of the Wendland kernel decreases steadily with the number of quadrature points and with the size of the cells. For example, on a 3D mesh of Q1 elements with :math:`R/h = 2.8` and two Gauss points per direction, the error of the kernel mass is :math:`5 \times 10^{-2}` for the top-hat kernel and :math:`7 \times 10^{-4}` for the Wendland kernel of the same support radius, and the latter drops to :math:`10^{-6}` with four Gauss points per direction.

Discretization
--------------

Filter centers
~~~~~~~~~~~~~~

The averaged fields are represented in a continuous ``FE_Q`` space of degree :math:`k` defined on the mesh of the simulation (by default, the degree of the velocity). The filter centers are the support points of this space. Each process owns the filter centers of its locally owned and unconstrained degrees of freedom. The values at the hanging nodes are interpolated from the hanging node constraints, and the degrees of freedom located on two periodic boundaries are two filter centers that evaluate the same average. The averaged fields are thus genuine finite element fields: they can be visualized, interpolated, integrated and differentiated like the solution itself.

Quadrature of the convolution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The integrals over the domain are split over the cells :math:`K` of the mesh and evaluated with a quadrature rule of points :math:`\mathbf{y}_q` and weights :math:`w_q` (including the Jacobian of the mapping). For example, the fluid volume fraction at the filter center :math:`\mathbf{x}_t` is

.. math::
    \varepsilon_f(\mathbf{x}_t) \approx \sum_K \sum_q I_f(\mathbf{y}_q) \, g(\lVert \mathbf{x}_t-\mathbf{y}_q \rVert) \, w_q,

and the velocity and pressure are evaluated from the finite element solution at the quadrature points. Every filtered quantity is a combination of four moments accumulated at each filter center: :math:`M`, :math:`E`, :math:`\int I_f \mathbf{u} g` and :math:`\int I_f p g`. The kernel is evaluated from the squared distance, which avoids a square root, and returns zero beyond its support. The quadrature of the kernel is accurate when its support contains many cells, i.e. when :math:`R \gg h`, where :math:`h` is the size of the cells.

Cells cut by the particles
~~~~~~~~~~~~~~~~~~~~~~~~~~

The fluid indicator is discontinuous across the surface of the particles, and a regular Gauss quadrature of a cut cell would misplace the interface by a fraction of the cell. Each cell is therefore classified as fluid, solid or cut from the signed distance functions :math:`\phi_i` of the particles (negative inside). The cell is enclosed in the ball of center :math:`\mathbf{c}` and radius :math:`r` that circumscribes its bounding box. Since a signed distance function is 1-Lipschitz:

* if :math:`\phi_i(\mathbf{c}) < -r` for a particle :math:`i`, the cell lies entirely inside that particle;
* if :math:`\phi_i(\mathbf{c}) > r` for every particle, the cell lies entirely in the fluid;
* otherwise, the cell may be cut by the particles for which :math:`|\phi_i(\mathbf{c})| \leq r`.

This classification is conservative: a cell is only treated as fluid or solid if this is guaranteed. The candidate particles of a cell are found with an R-tree of the bounding boxes of the particles, so that the cost of the classification does not grow with the number of particles. Fluid cells use a Gauss quadrature with :math:`I_f = 1`. Solid cells use the same quadrature with :math:`I_f = 0`, since they still contribute to the kernel mass. Cut cells use an iterated Gauss quadrature, made of :math:`n_s^d` copies of the Gauss rule on a regular subdivision of the cell, and the fluid indicator is evaluated at each of its points with the signed distance functions of the candidate particles only. The error on the volume of the particles decreases with the number of subdivisions :math:`n_s`, while the cost of a cut cell grows as :math:`n_s^d`.

Parallel algorithm
------------------

The naive way to compute the filtered fields is target-centric: for each filter center, integrate the solution over the ball of radius :math:`R` around it. On a distributed mesh, this ball spans many subdomains when :math:`R` is much larger than the cells, so each filter center would need to evaluate the solution at points owned by other processes, or each process would need ghost layers of width :math:`R`. Alternatively, the filter could be assembled as a sparse matrix mapping the solution to the filtered fields, but each of its rows would contain one entry per degree of freedom within the support, i.e. :math:`\mathcal{O}\left((R/h)^d\right)` entries, which quickly exceeds the memory used by the simulation itself.

The filter instead uses a source-centric formulation. The integrals are sums of independent contributions of the cells, so each process can integrate the cells it owns, whose solution values it already has, and compute *partial* moments for every filter center within the support of the kernel of these cells. The partial moments of a filter center are then summed by its owner. The solution never leaves its process: only the positions of the filter centers and four partial moments per filter center are communicated. The algorithm proceeds in the following stages.

.. code-block:: text

    1. Generate the owned filter centers (owned, unconstrained degrees of freedom)
    2. Describe the owned cells of every process by a few bounding boxes,
       and gather these boxes on all processes in a global R-tree
    3. Send each owned filter center, and its periodic images, to the processes
       whose boxes lie within the support radius             (sparse exchange)
    4. Build a local R-tree of the filter centers received by the process
    5. Sweep the owned cells once, in batches of quadrature points, and
       accumulate partial moments at the filter centers within the support
    6. Return the non-zero partial moments to the owners      (sparse exchange)
    7. Sum the moments, normalize them, and fill the hanging nodes and ghosts

Description of the subdomains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each process computes the bounding box of each of its cells, as given by the mapping and slightly enlarged to account for curved cells, and groups them into a regular grid of at most :math:`4^d` bins covering its subdomain. Each bin is described by the bounding box of its cells, so a process is described by at most 16 boxes in 2D and 64 boxes in 3D, which follow the shape of its subdomain even when it is not convex or not connected. These boxes are gathered on every process and stored in an R-tree, which indexes the region of the domain covered by the quadrature points of every process.

Distribution of the filter centers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A filter center :math:`\mathbf{x}_t` receives contributions from a process if a quadrature point of that process lies within the distance :math:`R` of :math:`\mathbf{x}_t`. Since every quadrature point lies in the boxes of its process, it is sufficient to send :math:`\mathbf{x}_t` to the processes that have a box within the distance :math:`R`. These processes are found by querying the global R-tree with the box of half-width :math:`R` around :math:`\mathbf{x}_t`, and discarding the boxes whose exact distance to :math:`\mathbf{x}_t` is larger than :math:`R`. The filter centers are then exchanged with a sparse point-to-point communication, in which a process only communicates with the processes that share filter centers with it. Each request contains the position of the filter center and its index on its owner, so that the partial moments can later be returned without any search. The filter centers of a process that fall within the support of its own cells are kept locally without communication.

This distribution is complete, as every quadrature point that can contribute to a filter center reaches it, and it only sends a filter center to the processes that can contribute to it, up to the difference between the boxes and the subdomains. Each process then holds a halo of filter centers around its subdomain, of thickness :math:`R`.

Periodic boundaries
~~~~~~~~~~~~~~~~~~~

Across a periodic boundary, a filter center must gather the quadrature points located on the other side of the domain. This is achieved by distributing periodic images of the filter centers. The translation :math:`\mathbf{L}` of each pair of periodic boundaries is computed from the coarse mesh, which every process stores, so that it is known on every process without communication, including the processes that do not own any cell at the periodic boundaries. A filter center has an image translated by :math:`k \mathbf{L}` for every non-zero integer :math:`k` such that the image lies within the distance :math:`R` of the domain, and the images along several periodic directions are combined to reach the edges and corners of the domain. Each image is distributed like a regular filter center, and its partial moments are returned to the owner of the original filter center.

The support radius is not limited by the period. Since the fields are periodic, the average over the infinite periodic medium is the integral over the domain of the periodized kernel:

.. math::
    \int_{\mathbb{R}^d} I_f(\mathbf{y}) \, a(\mathbf{y}) \, g(\lVert \mathbf{x}-\mathbf{y} \rVert) \, \mathrm{d}\mathbf{y} = \int_\Omega I_f(\mathbf{y}) \, a(\mathbf{y}) \sum_{\mathbf{k}} g(\lVert \mathbf{x} + \mathbf{k} \mathbf{L} - \mathbf{y} \rVert) \, \mathrm{d}\mathbf{y}.

When :math:`R \leq \lVert \mathbf{L} \rVert / 2`, a quadrature point contributes to a filter center through at most one of its images. When the support radius is larger, as is common for the small periodic domains used to derive closures, a quadrature point contributes through several images: these contributions are the terms of the sum above and are not counted twice. The mass of the periodized kernel over the domain remains one.

Batched sweep of the source cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each process builds an R-tree of the filter centers it received, including its own, and sweeps its owned cells once. Querying the tree for every quadrature point, or even for every cell, would make the search as expensive as the evaluation of the kernel when the cells are small compared to the support of the kernel. The quadrature points of consecutive cells are therefore gathered in batches. The cells are visited in the order of the triangulation, which follows the space-filling curve of the partitioner [#burstedde2011]_, so that consecutive cells are close to each other. A batch is processed when it reaches 1024 points, or when adding a cell would make it wider than :math:`0.1R`. For each batch:

1. the R-tree is queried once for the filter centers within the box of the batch enlarged by :math:`R`;
2. the candidates whose exact distance to the box of the batch is larger than :math:`R` are discarded;
3. for each remaining filter center, the kernel is evaluated at every point of the batch and the four moments are accumulated.

The cost of a search is thus shared by all the points of a batch. The points of a batch are stored as a structure of arrays, so that the innermost loop over the points is contiguous and can be vectorized. The width of a batch is limited because every filter center within :math:`R` of the batch evaluates the kernel at every point of the batch, including the points beyond its support. The ratio between the number of candidate filter centers and the number of filter centers within the support of a single point is the ratio of the volume of the box of width :math:`e` enlarged by :math:`R` to the volume of the ball of radius :math:`R`, i.e. approximately :math:`1 + 9e/4R` in 3D and :math:`1 + 4e/\pi R` in 2D. With :math:`e = 0.1R`, at most about 24% (3D) and 13% (2D) of the kernel evaluations are wasted, in exchange for one search per batch instead of one search per point. No memory is allocated and no communication occurs during the sweep.

Reduction and normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The partial moments of the filter centers received from other processes are returned to their owners with a second sparse exchange, which only contains the filter centers that received a non-zero contribution. The owner sums the contributions in the order of the ranks of the processes, so that the result does not depend on the order in which the messages arrive. The filtered fields are then computed from the moments, the values of the hanging nodes are interpolated, and the ghost values are updated.

Cost and scalability
~~~~~~~~~~~~~~~~~~~~

The cost of the filter is dominated by the evaluations of the kernel. For :math:`N_c` filter centers and a quadrature of :math:`n_q^d` points per cell, it scales as

.. math::
    W \approx N_c \, n_q^d \, \frac{|B_R|}{h^d},

i.e. as :math:`(R/h)^d` per filter center, with an additional factor :math:`n_s^d` for the cut cells. Doubling the width of the filter in 3D multiplies the cost by eight. The work of a process is proportional to the number of quadrature points it owns times the number of filter centers around them, so the partition of the mesh, which balances the number of cells, also balances the work, apart from the cut cells.

The communication consists of one gathering of the boxes of the subdomains and two sparse exchanges, whose volume is proportional to the number of halo filter centers. For a subdomain of width :math:`L_s`, the halo contains about :math:`\left((L_s + 2R)^d - L_s^d\right)/h^d` filter centers. When the support radius becomes larger than the subdomains, i.e. when many processes are used for a wide filter, the halo dominates the owned filter centers and bounds the strong scaling of the filter. The halo can be monitored with the ``extra verbose`` verbosity of the application.

Verification
------------

The unit tests of the filter verify properties that hold exactly for the discrete filter, on one, two and three processes:

* On a periodic mesh of uniform Q1 elements, every filter center sees the same arrangement of quadrature points, so the kernel mass must be identical at every filter center. A contribution missed across the subdomains or across the periodic boundaries breaks this uniformity. This also holds for a kernel wider than the period, whose filter centers have several periodic images per direction.
* On a uniform mesh with walls, the kernel mass of a filter center on a face, an edge or a corner is exactly 1/2, 1/4 or 1/8 of the kernel mass at the middle of the domain, by symmetry.
* With a top-hat kernel centered on a particle smaller than the kernel, the solid volume fraction times the kernel mass and the volume of the kernel equals the volume of the particle integrated by the quadratures of the source cells.
* The phase averages of a uniform field recover the field exactly wherever they are defined.
* On an adaptively refined mesh, the value of every filtered field at a hanging node is the interpolation of the values at the filter centers that constrain it, and a uniform field is recovered at every degree of freedom.
* For a Couette flow between a stationary and a moving wall, the averaged velocity is exact beyond the distance :math:`R` of the walls and antisymmetric about the middle of the channel, with and without the extension of the fluid beyond the walls. With the extension, the fluid and solid volume fractions sum to one up to the walls.

Limitations
-----------

* The classification of the cells assumes that the signed distance functions of the particles are 1-Lipschitz, which holds for the exact signed distances of analytical shapes.
* Like the sharp-interface immersed boundary solver, the filter does not consider the periodic images of the particles.
* The periodic boundaries must be planes normal to their direction of periodicity.
* Within the distance :math:`R` of a non-periodic boundary, the exact mass of the truncated kernel is not known, and the volume fractions carry the quadrature error of the kernel mass, unless they are renormalized. This error is small when the support of the kernel contains many cells.
* Within the distance :math:`R` of a wall, the kernel is truncated and off-centered, which reduces the gradient of the averaged velocity, even when the fluid is extended beyond the walls (see above).
* Only quadrilateral and hexahedral meshes are supported.

References
----------

.. [#anderson1967] \T. B. Anderson and R. Jackson, "Fluid mechanical description of fluidized beds. Equations of motion," *Industrial & Engineering Chemistry Fundamentals*, vol. 6, no. 4, pp. 527–539, 1967, doi: `10.1021/i160024a007 <https://doi.org/10.1021/i160024a007>`_\.
.. [#jackson2000] \R. Jackson, *The Dynamics of Fluidized Particles*. Cambridge University Press, 2000.
.. [#wendland1995] \H. Wendland, "Piecewise polynomial, positive definite and compactly supported radial functions of minimal degree," *Advances in Computational Mathematics*, vol. 4, pp. 389–396, 1995, doi: `10.1007/BF02123482 <https://doi.org/10.1007/BF02123482>`_\.
.. [#burstedde2011] \C. Burstedde, L. C. Wilcox, and O. Ghattas, "p4est: Scalable algorithms for parallel adaptive mesh refinement on forests of octrees," *SIAM Journal on Scientific Computing*, vol. 33, no. 3, pp. 1103–1133, 2011, doi: `10.1137/100791634 <https://doi.org/10.1137/100791634>`_\.
