================================
The Volume of Fluid (VOF) Method
================================

Numerous examples of flow encountered in engineering involve multiple fluids: sloshing of fuel in aircraft tanks, mixing of bread dough, and motion of droplets and bubbles to name a few. In these cases, the involved fluids can be immiscible, and we are interested in the evolution of the interfaces between those fluids.

Let :math:`\Omega = \Omega_0 \cup \Omega_1` be the domain formed by two fluids, namely fluid :math:`0` and :math:`1`, with :math:`\Gamma` denoting their interface and :math:`\partial \Omega`, the remaining boundaries, as illustrated in the figure below. In the VOF method [#hirt1981]_, we define the scalar function :math:`\phi` as a phase indicator such that:

.. math::
  \phi =
  \begin{cases}
    0 \quad \forall \mathbf{x} \in \Omega_0\\
    1 \quad \forall \mathbf{x} \in \Omega_1
  \end{cases}

This phase indicator (or phase fraction) changes rapidly but smoothly from :math:`0` to :math:`1` at the interface such that :math:`\Gamma` is located at the iso-contour :math:`\phi=0.5`, as illustrated below.

.. image:: images/vof.png
    :alt: Schematic
    :align: center
    :width: 600

The evolution of :math:`\Gamma` (or the iso-contour :math:`\phi=0.5`) in the time interval :math:`[0,T]` due to the action of velocity field :math:`\mathbf{u}` on :math:`\Omega` is described by the advection equation of the field :math:`\phi`:

.. math::
  \frac{\partial \phi}{\partial t} + \nabla \cdot \left( \mathbf{u} \phi \right) = 0 \quad \forall (\mathbf{x},t)\in \Omega\times[0,T]

or using `Einstein notation <https://en.wikipedia.org/wiki/Einstein_notation>`_:

.. math::
  \partial_t \phi + \partial_i (\phi u_i) = 0 \quad \forall (x_i,t)\in \Omega\times[0,T]

Developing the second term gives:

.. math::
  \partial_t \phi + \phi\partial_i u_i + u_i\partial_i\phi = 0 \quad \forall (x_i,t)\in \Omega\times[0,T]

Typically, the term :math:`\phi\partial_i u_i` (or :math:`\phi \nabla \cdot \mathbf{u}`) is zero due to mass conservation in the Navier-Stokes equations. However, previous work done in Lethe showed that while :math:`\nabla \cdot \mathbf{u}=0` is globally respected, it is not locally respected, especially around the interface, so lets keep it for now.

To complete the strong formulation of the problem, let's impose a no flux boundary condition on :math:`\partial \Omega`:

.. math::
  (\partial_i \phi) n_i= 0 \quad \forall (x_i,t)\in \partial \Omega\times[0,T]

where :math:`n_i` represent the outward pointing unit normal vector of :math:`\partial \Omega`, i.e., :math:`\mathbf{n}`.

Finite Element Formulation
---------------------------

To obtain the finite element formulation, we first need the weak formulation. Therefore, let :math:`v` be an arbitrary scalar function on :math:`\Omega`. To obtain the weak form, we multiply the strong problem by :math:`v` and integrate over :math:`\Omega`:

.. math::
  \int_\Omega v \left( \partial_t \phi + \phi\partial_i u_i + u_i\partial_i\phi\right) d \Omega = 0

To ensure that those integrals are well defined in :math:`\Omega`, we chose the appropriate solution spaces:

.. math::
  \phi \in \Phi(\Omega) = \mathcal{H}^1(\Omega)

.. math::
  v \in V(\Omega) = \mathcal{L}^2(\Omega)

Thus, the weak problem is:

Find :math:`\phi \in \Phi(\Omega) \times [0,T]` such that

.. math::
  \int_\Omega v \left( \partial_t \phi + \phi\partial_i u_i + u_i\partial_i\phi\right) d \Omega = 0 \quad \forall v\in V

Using Petrov-Galerkin method, the finite element formulation reads:

Find :math:`\phi^h \in \Phi^h \times [0,T]` such that

.. math::
  \int_\Omega v^h \left( \partial_t \phi^h + \phi^h\partial_i u_i + u_i\partial_i\phi^h\right) d \Omega = 0 \quad \forall v^h\in V^h\times [0,T]

where :math:`\Phi^h` and :math:`V^h` are finite element spaces, and :math:`\phi^h(\mathbf{x},t) = \sum_{j=1}^N \phi_j(t)\psi_j(\mathbf{x})`. In standard notation, this formulation corresponds to:

Find :math:`\phi^h \in \Phi^h \times [0,T]` such that

.. math::
  \int_\Omega v^h \left(\frac{\partial \phi^h}{\partial t} + \phi^h \nabla \cdot \mathbf{u} + \mathbf{u}\cdot \nabla \phi^h \right) d \Omega = 0 \quad \forall v^h\in V^h\times [0,T]

Stabilization
--------------

The numerical resolution of the advection equation requires stabilization because of its purely advective character, which makes the equation hyperbolic. Furthermore, a second stabilization term is added to improve the capturing of the interface due to sharp gradient across :math:`\Gamma`. Since SUPG only adds diffusion along the streamlines, crosswind oscillations may occur if no appropriate shock capturing scheme is used. To that end, a Discontinuity-Capturing Directional Dissipation (DCDD) shock capturing scheme is used [#tezduyar2003]_:

.. math::

  &\int_\Omega v^h \left( \partial_t \phi^h + \phi^h\partial_i u_i + u_i\partial_i\phi^h\right) d \Omega \\
  &\quad + \sum_k \int_{\Omega_k}\tau_\mathrm{SUPG} (u_i\partial_i v^h)\left(\partial_t \phi^h + \phi^h\partial_i u_i + u_i\partial_i\phi^h \right) d \Omega_k \\
  &\qquad + \sum_k \int_{\Omega_k}v_\mathrm{DCDD} (\partial_i v^h)( f_{\mathrm{DCDD}_ij}\partial_i \phi^h)  d \Omega_k  = 0

where the first element-wise summation represents the SUPG stabilization term and the second is the shock capturing scheme. The same SUPG stabilization as in the Navier-Stokes finite element formulation is used (see :doc:`../../multiphysics/fluid_dynamics/stabilization`). The terms of the DCDD scheme are:  

.. math::

  &v_\mathrm{DCDD} = \frac{1}{2} h^2 \|\mathbf{u}\| \| \nabla \phi^h_\mathrm{old} \| \\
  &\mathbf{f}_\mathrm{DCDD} =
    \frac{\mathbf{u}}{\|\mathbf{u}\| } \otimes \frac{\mathbf{u}}{\|\mathbf{u}\|}
    - \left(\frac{\nabla \phi^h_\mathrm{old}}{\| \nabla \phi^h_\mathrm{old} \|} \cdot \frac{\mathbf{u}}{\|\mathbf{u}\| } \right)^2
    \frac{\nabla \phi^h_\mathrm{old}}{\| \nabla \phi^h_\mathrm{old} \|} \otimes \frac{\nabla \phi^h_\mathrm{old}}{\| \nabla \phi^h_\mathrm{old} \|}

The term :math:`v_\mathrm{DCDD}` ensures that diffusivity is added only where there is a large phase gradient and a non-zero velocity, i.e., where the interface :math:`\Gamma` is in motion. The term :math:`\mathbf{f}_\mathrm{DCDD}` adds diffusivity only in the crosswind direction, since streamline diffusion is already added by the SUPG stabilization.

To avoid a non-linear finite element formulation, the phase gradient of the previous time step :math:`(\phi^h_\mathrm{old})` is used.

Interface Diffusion and Regularization
--------------------------------------

The VOF method tends to diffuse the interface, i.e., over time, the interface becomes blurry instead of a sharp definition, and the change from :math:`\phi = 0` to :math:`1` occurs on a larger distance.

Thus, we use regularization methods to keep the change in :math:`\phi` sharp at the interface. Three methods are currently available: projection-based interface sharpening, algebraic interface reinitialization and interface filtration.

""""""""""""""""""""""""""""""""""""""
Projection-Based Interface Sharpening
""""""""""""""""""""""""""""""""""""""

The current projection-based interface sharpening method consists of two steps:


1. Phase fraction limiter

.. math::

    \phi = \min \left( \max \left(\phi_\mathrm{old},0 \right),1 \right)

The phase fraction limiter above will update the phase fraction if it failed to respect these bounds.


2. Interface sharpening using a projection

.. math::

    \phi =
    \begin{cases}
    c^{1-\alpha} \phi^{\alpha} &  (0 \leq \phi < c  ) \\
    1-(1-c)^{1-\alpha}(1-\phi)^{\alpha} & (c \leq \phi \leq 1  )
    \end{cases}

where :math:`c` denotes the sharpening threshold, which defines
a phase fraction threshold (generally :math:`0.5`), and :math:`\alpha` corresponds to the interface sharpness, which is a model parameter generally in the range of :math:`(1,2]`. This projection-based interface sharpening method was proposed by Aliabadi and Tezduyar (2000) [#aliabadi2000]_.

""""""""""""""""""""""""""""""""""""""""
Geometric Interface Reinitialization
""""""""""""""""""""""""""""""""""""""""

The geometric interface reinitialization implemented in Lethe uses the signed distance :math:`d` from the interface to regularize the phase fraction field.  The method is based on the work of Ausas *et al.* (2011) [#ausas2011]_, originaly proposed in a level-set framework. Once computed, the signed distance is transformed into a phase fraction field using a transformation function :math:`g` such as :math:`\phi = g(d)`.

To compute the signed distance, the interface is linearly reconstructed from the iso-contour :math:`\phi=0.5` using the `Marching Cube algorithm implemented in deal.II <https://dealii.org/current/doxygen/deal.II/classGridTools_1_1MarchingCubeAlgorithm.html>`_. Then, the signed distance is computed layer-by-layer, from the interface until the user-defined maximum distance :math:`d_\mathrm{max}` is reached on each side of the interface. 

For the first layer, the analytical minimum distance between the DoFs of the cells cut by the interface and the reconstructed interface is computed. It comes down to the computation of the distance between points (i.e., the DoFs) and the line segments in 2D or triangles in 3D that approximate the interface.

For the subsequent layers, the signed distance for the remaining DoFs is computed iteratively in the narrow band defined by :math:`d \in [-d_\mathrm{max}, d_\mathrm{max}]`, by solving a minimization problem. We refer the reader to the work of Ausas *et al.* (2011) [#ausas2011]_ for more details. The signed distance for the DoFs outside of the narrow band is set to :math:`\pm d_\mathrm{max}`, where the sign depends on the side of the interface on which they are located.

Finally, the signed distance field is transformed to a phase fraction field. Here, we want that:

.. math::
  \phi \approx
  \begin{cases}
    1 \quad \text{if } d < -d_\mathrm{max} \\
    0.5 \quad \text{if } d = 0\\
    0 \quad \text{if } d > d_\mathrm{max} \\
  \end{cases}

In Lethe, two functions are available to achieve that: a hyperbolic tangent function or a 4th degree, piecewise polynomial. 

* hyperbolic tangent: the regularized phase fraction is given by

  .. math::
    \phi = 0.5-0.5\tanh(d/\varepsilon)
    
  where :math:`\varepsilon` is a measure of the interface thickness. Note that this transformation does not ensure that :math:`\phi=0` or :math:`1` when :math:`d = \pm d_\mathrm{max}`. The value of :math:`\phi` at :math:`\pm d_\mathrm{max}` depends on :math:`d_\mathrm{max}` and :math:`\varepsilon`.

* piece-wise polynomial: this transformation takes the form

  .. math::
    \phi =
    \begin{cases}
      0.5 - 0.5(4d' + 6d'^2 + 4d'^3 + d'^4) \text{ if } d' < 0.0 \\
      0.5 - 0.5(4d' - 6d'^2 + 4d'^3 - d'^4) \text{ if } d' > 0.0
    \end{cases}
  
  where :math:`d' = d/d_\mathrm{max}` is the dimensionless distance. Contrary to the hyperbolic tangent function, it ensures that :math:`\phi=0` or :math:`1` when :math:`d = \pm d_\mathrm{max}`.

""""""""""""""""""""""""""""""""""""""""
Algebraic Interface Reinitialization
""""""""""""""""""""""""""""""""""""""""

The algebraic interface reinitialization method consists of compressing and diffusing the interface in its normal direction. This is done by solving the following transient Partial Differential Equation (PDE) until steady-state is reached using a pseudo-time-stepping scheme as proposed by Olsson and coworkers (2007) [#olsson2007]_:

.. math::

    \underbrace{\frac{\partial\phi_\text{reinit}}{\partial \tau}}_\text{transient}
     + \underbrace{\nabla \cdot \left[ \phi_\text{reinit} (1-\phi_\text{reinit})\mathbf{n}\right]}_\text{compression}
    = \underbrace{\varepsilon \nabla \cdot \left[  (\nabla \phi_\text{reinit} \cdot \mathbf{n}) \mathbf{n}\right]}_\text{diffusion}


where:

- :math:`\phi_\text{reinit}` is the reinitialized phase fraction;

- :math:`\tau` is the pseudo-time independent variable. It is different from the time independent variable :math:`t` of the actual simulation.

- :math:`\mathbf{n} = \frac{\nabla \psi}{\lVert \nabla \psi \rVert}` is the normal vector of the interface with :math:`\nabla \psi` the :ref:`projected VOF phase gradient<Normal and curvature computations>`, and;

- :math:`\varepsilon = C \, h_\text{min}^d` is the diffusion coefficient with  :math:`h_\text{min}` the smallest cell-size, :math:`C` a constant factor multiplying :math:`h_\text{min}`, and :math:`d` a constant power to which :math:`h_\text{min}` is elevated. As default, :math:`C` and :math:`d` are set to :math:`1`.

.. note::

    :math:`\nabla \psi` is computed with the VOF phase fraction gradient field and remains constant through the interface reinitialization process of a same global time iteration.

.. note::

    Here, we define the cell-size as being the diameter of:

    - a disk of equivalent area in 2D, and;
    - a sphere of equivalent volume in 3D.

The equation is solved using the finite element method:

Considering :math:`\upsilon \in V(\Omega)` an arbitrary scalar function on :math:`\Omega` as the test function and :math:`\phi_\text{reinit} \in \Phi(\Omega)`, the weak problem goes as follows:

Find :math:`\phi_\text{reinit} \in \Phi(\Omega) \times [0, \mathrm{\tau}_\text{end}]`

.. math::

    \int_\Omega \upsilon \left(\partial_\tau \phi_\text{reinit} + \nabla \cdot \left[ \phi_\text{reinit} (1-\phi_\text{reinit})\mathbf{n}\right] - \varepsilon \nabla \cdot \left[  (\nabla \phi_\text{reinit} \cdot \mathbf{n}) \mathbf{n}\right]\right)\, \mathrm{d}\Omega = 0 \quad \forall \, \upsilon \in V

As the equation is non-linear, we use the `Newton-Raphson method <https://en.wikipedia.org/wiki/Newton%27s_method>`_ and solve the linear system :math:`\mathcal{J} \, \delta\phi_\text{reinit} = -\mathcal{R}`. We obtain the Residual (:math:`\mathcal{R}`) and Jacobian (:math:`\mathcal{J}`) below:

.. math::
    \begin{split}
    \mathcal{R} = & \int_\Omega \upsilon \, \partial_\tau \phi_\text{reinit}^{n} \, \mathrm{d}\Omega \\
                   & + \int_\Omega\upsilon
                  \left[ \mathbf{n} \cdot \nabla \phi_\text{reinit}^{n} - \mathbf{n} \cdot \left(  2 \phi_\text{reinit}^{n} \nabla \phi_\text{reinit}^{n}\right) +(\phi_\text{reinit}^{n} -{\left(\phi_\text{reinit}^{n}\right)}^2) (\nabla \cdot \mathbf{n}) \right] \mathrm{d}\Omega \\
                   & + \int_\Omega \nabla\upsilon \cdot \varepsilon (\nabla \phi_\text{reinit}^{n} \cdot \mathbf{n}) \mathbf{n} \, \mathrm{d} \Omega

    \end{split}

where, :math:`\phi_\text{reinit}^{n}` is the reinitialized phase fraction value of the previous Newton iteration.

Considering,

.. math::

    \delta \phi_\text{reinit} = \sum_{j=1}^N \delta \phi_{\text{reinit},j} \, \xi_j

where :math:`\xi_j` is the :math:`j\text{th}` interpolation function of the reinitialized phase fraction field, the Jacobian reads:

.. math::
    \begin{split}
    \mathcal{J} = & \int_\Omega \upsilon \, \partial_\tau \xi_j \, \mathrm{d}\Omega \\
                   & + \int_\Omega\upsilon
                  \left[ \mathbf{n} \cdot \left( \nabla \xi_j - 2 \phi_\text{reinit}^{n} \nabla \xi_j -2 \xi_j \nabla\phi_\text{reinit}^{n} \right) + \left( \xi_j -2\phi_\text{reinit}^{n}\xi_j \right) \right] \mathrm{d}\Omega \\
                   & + \int_\Omega \nabla\upsilon \cdot \varepsilon (\nabla \xi_j \cdot \mathbf{n}) \mathbf{n} \,\mathrm{d}\Omega

    \end{split}

and the reinitialized phase fraction is given by:

.. math::
    \phi_\text{reinit}^{n+1} = \phi_\text{reinit}^{n} + \delta \phi_\text{reinit}


""""""""""""""""""""""""""""""""
Interface Filtration
""""""""""""""""""""""""""""""""

In the interface filtration method, the following filter function is applied to the phase fraction :math:`\phi` in order to get a better definition of the interface between the fluids:

.. math::
    \phi' = 0.5 \tanh[\beta(\phi-0.5)] + 0.5

where :math:`\phi'` is the filtered phase fraction value, and :math:`\beta` is a model parameter that enables sharper definition when increased. Recommended value is :math:`\beta=20`.

Surface Tension
---------------

When two immiscible fluids are in contact, surface tension tends to deform their interface (also called the free surface) into a shape that ensures a minimal energy state. An example would be the force that drives a droplet into its spherical shape [#brackbill1992]_.

Resolution of the interface motion via the advection equation allows to compute the surface tension term and add its effect to the Navier-Stokes momentum equations.

As its name suggests, the surface tension :math:`\bf{f_{\sigma}}` is a surface force. It is applied at the interface between two immiscible fluids and is given by:

.. math::

    {\bf{f_{\sigma}}} = \sigma \kappa {\bf{n}}

where :math:`\sigma` is the surface tension coefficient, :math:`\kappa` is the curvature and :math:`\bf{n}` is the unit normal vector of the free surface. Here, :math:`{\bf{f_{\sigma}}}` is a force per unit of area. To account for its effect in the Navier-Stokes equations, the surface force is transformed in a volumetric surface force :math:`\bf{F_{\sigma}}` using the continuous surface force (CSF) model [#brackbill1992]_, that is:

.. math::

    {\bf{F_{\sigma}}} = \bf{f_{\sigma}} \delta = \sigma \kappa {\bf{n}}\delta

where :math:`\delta` is a Dirac delta measure with support on the interface. A good approximation for the term :math:`{\bf{n}}\delta` is :math:`{\bf{n}}\delta = \nabla \phi`, where :math:`\phi` is the phase fraction. Thus, the volumetric surface force is given by:

.. math::

    {\bf{F_{\sigma}}} =  \sigma \kappa \nabla \phi

where the curvature :math:`\kappa` is computed according to:

.. math::

    \kappa = - \nabla \cdot \bf{n}

and the unit normal vector of the free surface, pointing from fluid 0 to 1, is obtained with:

.. math::

    \bf{n} = \frac{\nabla \phi}{|\nabla \phi|}

When including the surface tension force in the resolution of the Navier-Stokes equations, the numerical computation of the curvature can give rise to parasitic flows near the interface between the two fluids. To avoid such spurious currents, the phase fraction gradient and curvature are filtered using projection steps [#zahedi2012]_, as presented in section :ref:`Normal and curvature computations`.

.. _Normal and curvature computations:

"""""""""""""""""""""""""""""""""
Normal and Curvature Computations
"""""""""""""""""""""""""""""""""

The following equations are used to compute the filtered phase fraction gradient and filtered curvature. They correspond to the projection steps previously mentioned.

.. math::

    \int_\Omega \left( {\bf{v}} \cdot {\bf{\psi}} + \eta_n \nabla {\bf{v}} \cdot \nabla {\bf{\psi}} \right) d\Omega = \int_\Omega \left( {\bf{v}} \cdot \nabla {\phi} \right) d\Omega

where :math:`{\bf{v}}` is a vector test function, :math:`\bf{\psi}` is the filtered phase fraction gradient, :math:`\eta_n` is the phase fraction gradient filter value, and :math:`\phi` is the phase fraction.

.. math::

    \int_\Omega \left( v \kappa + \eta_\kappa \nabla v \cdot \nabla \kappa \right) d\Omega = \int_\Omega \left( \nabla v \cdot \frac{\bf{\psi}}{|\bf{\psi}|} \right) d\Omega

where :math:`\kappa` is the filtered curvature, and :math:`\eta_\kappa` is the curvature filter value, and :math:`v` is a test function.

The phase fraction gradient filter :math:`\eta_n` and the curvature filter value :math:`\eta_\kappa` are respectively computed according to:

.. math::

  \eta_n = \alpha h^2

  \eta_\kappa = \beta h^2

where :math:`\alpha` and :math:`\beta` are user-defined factors, and :math:`h` is the cell size. Recommended values are :math:`\alpha = 4.0` and :math:`\beta = 1.0`.


References
-----------

.. [#hirt1981] \C. W. Hirt and B. D. Nichols, “Volume of fluid (VOF) method for the dynamics of free boundaries,” *J. Comput. Phys.*, vol. 39, no. 1, pp. 201–225, Jan. 1981, doi: `10.1016/0021-9991(81)90145-5 <https://doi.org/10.1016/0021-9991(81)90145-5>`_\.

.. [#tezduyar2003] \T. E. Tezduyar, “Computation of moving boundaries and interfaces and stabilization parameters,” *Int. J. Numer. Methods Fluids*, vol. 43, no. 5, pp. 555–575, 2003, doi: `10.1002/fld.505 <https://doi.org/10.1002/fld.505>`_\.

.. [#aliabadi2000] \S. Aliabadi and T. E. Tezduyar, “Stabilized-finite-element/interface-capturing technique for parallel computation of unsteady flows with interfaces,” *Comput. Methods Appl. Mech. Eng.*, vol. 190, no. 3, pp. 243–261, Oct. 2000, doi: `10.1016/S0045-7825(00)00200-0 <https://doi.org/10.1016/S0045-7825(00)00200-0>`_\.

.. [#olsson2007] \E. Olsson, G. Kreiss, and S. Zahedi, "A conservative level set method for two phase flow II," *J. Comput. Phys.*, vol. 225, no. 1, pp. 785–807, Jul. 2007, doi: `10.1016/j.jcp.2006.12.027 <https://doi.org/10.1016/j.jcp.2006.12.027>`_\.

.. [#brackbill1992] \J. U. Brackbill, D. B. Kothe, and C. Zemach, “A continuum method for modeling surface tension,” *J. Comput. Phys.*, vol. 100, no. 2, pp. 335–354, Jun. 1992, doi: `10.1016/0021-9991(92)90240-Y <https://doi.org/10.1016/0021-9991(92)90240-Y>`_\.

.. [#zahedi2012] \S. Zahedi, M. Kronbichler, and G. Kreiss, “Spurious currents in finite element based level set methods for two-phase flow,” *Int. J. Numer. Methods Fluids*, vol. 69, no. 9, pp. 1433–1456, 2012, doi: `10.1002/fld.2643 <https://doi.org/10.1002/fld.2643>`_\.

.. [#ausas2011] \R.F. Ausas, E.A. Dari, and G.C. Buscaglia, “A geometric mass-preserving redistancing scheme for the level set function,” *Int. J. Numer. Methods Fluids*, vol. 65, no. 8, pp. 989-1010, 2011, doi: `10.1002/fld.2227 <https://doi.org/10.1002/fld.2227>`_\.
