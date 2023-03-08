==========================
Static bubble
==========================

This example simulates a `two-dimensional static bubble`_.

.. _two-dimensional static bubble: https://onlinelibrary.wiley.com/doi/full/10.1002/fld.2643


----------------------------------
Features
----------------------------------
- Solver: ``gls_navier_stokes_2d``
- Two phase flow handled by the Volume of fluids (VOF) approach with surface tension force
- Calculation of filtered phase fraction gradient and curvature fields
- Unsteady problem handled by an adaptive BDF1 time-stepping scheme

---------------------------
Files used in this example
---------------------------
- Parameter file: ``examples/multiphysics/static-bubble/static_bubble.prm``
- Python file to generate plot: ``examples/multiphysics/static-bubble/static_bubble.py``

-----------------------------
Description of the case
-----------------------------

A circular bubble of radius :math:`R=0.5` is in equilibrium in the center of a two-dimensional squared domain of side length :math:`L=5.0` filled with air. The gravitational force is neglected, such as in a microgravity environment, and the ratio of density between the droplet and the air is 1, meaning that buoancy is also neglected. Therefore, without any external force, the bubble and the air are at rest and only the surface tension effects are involved, maintaining the droplet in its circular shape. The following schematic describes the geometry and dimensions of the simulation in the :math:`(x,y)` plane:

.. image:: images/static-bubble.png
    :alt: Schematic
    :align: center
    :width: 400


""""""""""""""""""""""""""""""""
Surface tension force
""""""""""""""""""""""""""""""""
As its name suggests, the surface tension is a surface force. It is applied at the interface between two immisible fluids and is given by :

.. math::

    {\bf{f_{\sigma}}} = \sigma \kappa {\bf{n}}

where :math:`\sigma`, :math:`\kappa` is the curvature and :math:`\bf{n}` is the normal vector of the free surface. Here, :math:`{\bf{f_{\sigma}}}` is a force per unit of area. To account for its effect in the Navier-Stokes equations, the surface force is transformed in a volumetric surface force using the continuous surface force (CSF) model [`1 <https://doi.org/10.1016/0021-9991(92)90240-Y>`_], that is :

.. math::

    {\bf{F_{\sigma}}} = \bf{f_{\sigma}} \delta = \sigma \kappa {\bf{n}}\delta

where :math:`\delta` is a Dirac delta measure with support on the interface. A good approximation for the term :math:`{\bf{n}}\delta` is :math:`{\bf{n}}\delta = \nabla \phi`, where :math:`\phi` is the phase fraction. Thus, the volumetric surface force is given by :

.. math::

    {\bf{F_{\sigma}}} =  \sigma \kappa \nabla \phi

where the curvature :math:`\kappa` is computed according to:

.. math::

    \kappa = - \nabla \cdot \bf{n}

and the unit normal vector of the free surface is obtained with:

.. math::

    \bf{n} = \frac{\nabla \phi}{|\phi|}

When including the surface tension force in the resolution of the Navier-Stokes equations, the numerical computation of the curvature can give rise to parasitic flows near the interface between the two fluids. To avoid such spurious currents, the phase fraction gradient and curvature are filtered using projection steps, as presented in section :ref:`Normal and curvature computations`.

The static bubble case is a relevant case to study the spurious currents, since the analyical solution is zero for the velocity. Therefore, nonzero velocities in the computed velocity field are considered as spurious currents [`2 <https://doi.org/10.1002/fld.2643>`_]. The analytical pressure drop between the interior (:math:`p_{int}`) and exterior (:math:`p_{ext}`) of the bubble is given by the Young-Laplace relation:

.. math::

    \Delta p = p_{int} - p_{ext} = \sigma \kappa

with the analyical curvature of the bubble : :math:`\kappa = 1/R`. This example is based on the static droplet case reported in [`2 <https://doi.org/10.1002/fld.2643>`_].

.. _Normal and curvature computations:

"""""""""""""""""""""""""""""""""
Normal and curvature computations
"""""""""""""""""""""""""""""""""
The following equations calculate the filtered phase fraction gradient and filtered curvature, respectively.


.. math::

    \int_\Omega \left( {\bf{v}} \cdot {\bf{\psi}} + \eta_n \nabla {\bf{v}} \cdot \nabla {\bf{\psi}} \right) d\Omega = \int_\Omega \left( {\bf{v}} \cdot \nabla {\phi} \right) d\Omega

where :math:`{\bf{v}}` is a vector test function, :math:`\bf{\psi}` is the filtered phase fraction gradient, :math:`\eta_n` is the phase fraction gradient filter value, and :math:`\phi` is the phase fraction.

.. math::

    \int_\Omega \left( v k + \eta_k \nabla v \cdot \nabla k \right) d\Omega = \int_\Omega \left( \nabla v \cdot \frac{\bf{\psi}}{|\bf{\psi}|} \right) d\Omega

where :math:`k` is the filtered curvature, and :math:`\eta_k` is the curvature filter value, and :math:`v` is a test function.

--------------
Parameter file
--------------

Time integration is handled by a 1st order backward differentiation scheme `(bdf1)`, for a :math:`3~\text{s}` simulation time with an constant time step of :math:`0.005~\text{s}`.

.. note::
    This example uses an adaptive time-stepping method, where the
    time-step is modified during the simulation to keep the maximum value of the CFL condition below a given threshold. Using ``output frequency = 20`` ensures that the results are written every 20 iterations. Consequently, the time increment between each vtu file is not constant.

.. code-block:: text

    # --------------------------------------------------
    # Simulation Control
    #---------------------------------------------------
    subsection simulation control
      set method           = bdf1
      set time end         = 3
      set time step        = 0.005
      set output name      = static-bubble
      set output frequency = 20
      set output path      = ./output/
    end

The ``multiphysics`` subsection enables to turn on `(true)`
and off `(false)` the physics of interest. Here ``VOF`` is chosen. The ``surface tension force`` are enabled in the VOF subsection.


.. code-block:: text

    #---------------------------------------------------
    # Multiphysics
    #---------------------------------------------------
    subsection multiphysics
      set VOF = true
    end

""""""""""""""""""""""""""""""""
Volume of Fluid (VOF)
""""""""""""""""""""""""""""""""

.. code-block:: text

    #---------------------------------------------------
    # VOF
    #---------------------------------------------------

    subsection VOF
      subsection surface tension force
        set enable                         = true
        set surface tension coefficient    = 1
        set phase fraction gradient filter factor = 4
        set curvature filter factor               = 1
        set output auxiliary fields        = true
      end
    end



""""""""""""""""""""""""""""""""
Initial conditions
""""""""""""""""""""""""""""""""
In the ``initial condition``, the initial velocity and initial position
of the droplet are defined. The droplet is initially
defined as a circle with a radius :math:`R= 0.5` at :math:`(x,y)=(0.0, 0.0)`. We enable the use of a projection step to ensure that the initial phase distribution
sufficiently smooth and avoid a staircase respresentation of the interface.

.. code-block:: text

    #---------------------------------------------------
    # Initial condition
    #---------------------------------------------------

    subsection initial conditions
      set type = nodal
      subsection uvwp
        set Function expression = 0; 0; 0
      end
      subsection VOF
        set Function expression = if ((x-0) * (x-0) + (y-0) * (y-0) < 0.50 * 0.5 , 1, 0)
        subsection projection step
          set enable = true
          set diffusion factor = 1
        end
      end
    end



-----------
References
-----------
`[1] <https://doi.org/10.1016/0021-9991(92)90240-Y>`_ Brackbill, J.U., Kothe, D.B. and Zemach, C., 1992. A continuum method for modeling surface tension. Journal of computational physics, 100(2), pp.335-354.

`[2] <https://doi.org/10.1002/fld.2643>`_ Zahedi, S., Kronbichler, M. and Kreiss, G., 2012. Spurious currents in finite element based level set methods for two‐phase flow. International Journal for Numerical Methods in Fluids, 69(9), pp.1433-1456.
