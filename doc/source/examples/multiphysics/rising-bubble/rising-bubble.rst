==========================
Rising bubble
==========================

This example simulates a `two-dimensional rising bubble`_. 

.. _two-dimensional rising bubble: https://onlinelibrary.wiley.com/doi/full/10.1002/fld.2643


----------------------------------
Features
----------------------------------
- Solver: ``gls_navier_stokes_2d`` 
- Two phase flow handled by the Volume of fluids (VOF) approach with phase filtering, phase sharpening, and surface tension force
- Calculation of filtered phase fraction gradient and curvature fields
- Unsteady problem handled by an adaptive BDF1 time-stepping scheme 
- Post-processing of a fluid barycentric coordinate and velocity


---------------------------
Files used in this example
---------------------------
- Parameter file: ``examples/multiphysics/rising-bubble/rising_bubble.prm``
- Python file to generate plot: ``examples/multiphysics/rising-bubble/rising_bubble.py``


-----------------------------
Description of the case
-----------------------------

A circular bubble with density of 100 and kinematic viscosity of 0.01 (all the units in this example are dimensionless) is defined at an initial location (0.5, 0.5) in a rectangular column filled with a denser fluid (with a density of 1000 and kinematic viscosity of 0.01). At :math:`t = 0` the bubble is released to rise inside the denser fluid column. The corresponding parameter file is 
``rising-bubble.prm``.

The following schematic describes the geometry and dimensions of the simulation in the :math:`(x,y)` plane:

.. image:: images/bubble-initial-configuration.png
    :alt: Schematic
    :align: center
    :width: 400

.. note:: 
    All the four boundary conditions are ``noslip``, and an external 
    gravity field of :math:`-0.98` is applied in the y direction.


--------------
Parameter file
--------------

Time integration is handled by a 1st order backward differentiation scheme `(bdf1)`, for a :math:`3~\text{s}` simulation time with an initial time step of :math:`0.001~\text{s}`.

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
      set time step        = 0.001
      set adapt            = true
      set max cfl          = 0.8
      set output name      = rising-bubble
      set output frequency = 20
      set output path      = ./output/
    end

The ``multiphysics`` subsection enables to turn on `(true)` 
and off `(false)` the physics of interest. Here ``VOF`` is chosen. The ``phase filtration``, ``interface sharpening``, and ``surface tension force`` are enabled in the VOF subsection.


.. code-block:: text

    #---------------------------------------------------
    # Multiphysics
    #---------------------------------------------------
    subsection multiphysics
      set VOF = true
    end 

The ``source term`` subsection defines the gravitational acceleration:

.. code-block:: text
    
    #---------------------------------------------------
    # Source term
    #---------------------------------------------------
    subsection source term
      set enable = true
      subsection xyz
        set Function expression = 0; -0.98; 0
      end
    end
    
""""""""""""""""""""""""""""""""
Volume of Fluid (VOF)
""""""""""""""""""""""""""""""""

In Lethe, the surface tension force (:math:`{\bf{F_{\sigma}}}`) is calculated using the continuous surface force (CSF) [1, 2]:

.. math::

    {\bf{F_{\sigma}}} = \sigma k \nabla {\phi}

where :math:`\sigma`, :math:`k` and :math:`\nabla {\phi}` denote respectively the surface tension coefficient, the filtered curvature and the phase fraction gradient. :math:`\rho`, :math:`\rho_1`, and :math:`\rho_2` are the density of the flow, the density of phase 0, and the density of phase 1, respectively.

The curvature :math:`k` is computed according to:

.. math::

    k = - \nabla \cdot \bf{n}

where :math:`\bf{n}` is the unit normal vector of the free surface. The latter is obtained with:

.. math::

    \bf{n} = \frac{\nabla \phi}{|\phi|}

When including the surface tension force in the resolution of the Navier-Stokes equations, the numerical computation of the curvature can give rise to parasitic flows near the interface between the two fluids. To avoid such spurious currents, the phase fraction gradient and curvature are filtered using L2-projections.
The following equations calculate the filtered phase fraction gradient and filtered curvature, respectively.

.. math:: 

    \int_\Omega \left( {\bf{v}} \cdot {\bf{\psi}} + \eta_n \nabla {\bf{v}} \cdot \nabla {\bf{\psi}} \right) d\Omega = \int_\Omega \left( {\bf{v}} \cdot \nabla {\phi} \right) d\Omega

where :math:`{\bf{v}}` is a vector test function, :math:`\bf{\psi}` is the filtered phase fraction gradient, :math:`\eta_n = \alpha h^2` is the phase fraction gradient filter value with :math:`h` denoting the cell size, and :math:`\phi` is the phase fraction.

.. math::

    \int_\Omega \left( v k + \eta_k \nabla v \cdot \nabla k \right) d\Omega = \int_\Omega \left( \nabla v \cdot \frac{\bf{\psi}}{|\bf{\psi}|} \right) d\Omega

where :math:`k` is the filtered curvature, and :math:`\eta_k = \beta h^2` is the curvature filter value, and :math:`v` is a test function.

.. tip::

  The phase fraction gradient filter value (:math:`\eta_n = \alpha h^2`) and curvature filter value (:math:`\eta_k = \beta h^2`) must be small values larger than 0. The values of :math:`\alpha` and :math:`\beta` are controlled respectively by the parameters ``phase fraction gradient filter factor`` and ``curvature filter factor``  in the parameter file.
  We recommend the following procedure to choose a proper value for these parameters:

  1. Enable ``output auxiliary fields`` to write filtered phase fraction gradient and filtered curvature fields.
  2. Choose a value close to 1, for example, the default values  :math:`\alpha = 4` and :math:`\beta = 1`.
  3. Run the simulation and check whether the filtered phase fraction gradient and filtered curvature fields are smooth and without oscillation.
  4. If the filtered phase fraction gradient and filtered curvature fields show oscillations, increase the value :math:`\alpha` and :math:`\beta` to larger values, and repeat this process until reaching smooth filtered phase fraction gradient and filtered curvature fields without oscillations. Generally, the default values should be sufficient.

The interface sharpening method and its parameters are explained in the :doc:`../dam-break/dam-break` example. We also enable phase filtration. This filters the phase field used for the calculation of physical properties by stiffening the value of the phase fraction. We refer the reader to the :doc:`../../../../parameters/cfd/volume_of_fluid` documentation for more explanation on the phase filtration.

.. code-block:: text

  #---------------------------------------------------
  # VOF
  #---------------------------------------------------
  subsection VOF
    subsection interface sharpening
      set enable              = true
      set threshold           = 0.5
      set interface sharpness = 1.5
      set frequency           = 50
    end

    subsection phase filtration
      set type      = tanh
      set verbosity = quiet
      set beta      = 10
    end

    subsection surface tension force
      set enable                                = true
      set surface tension coefficient           = 24.5
      set phase fraction gradient filter factor = 4
      set curvature filter factor               = 1
      set output auxiliary fields               = true
    end
  end

""""""""""""""""""""""""""""""""
Initial condition
""""""""""""""""""""""""""""""""
In the ``initial condition``, the initial velocity and initial position 
of the liquid phase are defined. The light phase is initially 
defined as a circle with a radius :math:`r= 0.25` at :math:`(x,y)=(0.5, 0.5)`. We enable the use of a projection step to ensure that the initial phase distribution 
sufficiently smooth.

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
        set Function expression = if ((x-0.5) * (x-0.5) + (y-0.5) * (y-0.5) < 0.25 * 0.25 , 1, 0)
      
        subsection projection step
          set enable           = true
          set diffusion factor = 1
        end
      end
    end


""""""""""""""""""""""""""""""""
Physical Properties
""""""""""""""""""""""""""""""""
We define two fluids here simply by setting the number of fluids to be :math:`2`.
In ``subsection fluid 0``, we set the density and the kinematic viscosity for the phase associated with a VOF indicator of 0. 
A similar procedure is done for the phase associated with a VOF indicator of 1 in ``subsection fluid 1``:


.. code-block:: text

    #---------------------------------------------------
    # Physical Properties
    #---------------------------------------------------
    subsection physical properties
      set number of fluids = 2
      subsection fluid 0
        set density             = 1000
        set kinematic viscosity = 0.01
      end
      subsection fluid 1
        set density             = 100
        set kinematic viscosity = 0.01
      end
    end



""""""""""""""""""""""""""""""""
Mesh
""""""""""""""""""""""""""""""""

We start off with a rectangular mesh that spans the domain defined by the corner points situated at the origin and at point
:math:`[1,2]`. The first :math:`1,2` couple defines that number of initial grid subdivisions along the length and height of the rectangle. 
This makes our initial mesh composed of perfect squares. We proceed then to redefine the mesh globally six times by setting
``set initial refinement = 6``. 

.. code-block:: text
        
    #---------------------------------------------------
    # Mesh
    #---------------------------------------------------
    subsection mesh
      set type               = dealii
      set grid type          = subdivided_hyper_rectangle
      set grid arguments     = 1, 2 : 0, 0 : 1, 2 : true
      set initial refinement = 6
    end
    
In the ``mesh adaptation subsection``, adaptive mesh refinement is 
defined for ``phase``. ``min refinement level`` and ``max refinement level`` are 6 and 9, respectively. Since the bubble rises and changes its location, we choose a rather large ``fraction refinement`` (0.99) and moderate ``fraction coarsening`` (0.01).
To capture the bubble adequately, we set ``initial refinement steps = 5`` so that the initial mesh is adapted to ensure that the initial condition is imposed for the VOF phase with maximal accuracy.

.. code-block:: text

    #---------------------------------------------------
    # Mesh Adaptation
    #---------------------------------------------------
    subsection mesh adaptation
      set type                     = kelly
      set variable                 = phase
      set fraction type            = fraction
      set max refinement level     = 9
      set min refinement level     = 6
      set frequency                = 1
      set fraction refinement      = 0.99
      set fraction coarsening      = 0.01
      set initial refinement steps = 5
    end

""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Post-processing: Fluid barycenter position and velocity
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To compare our simulation results to the literature, we extract the position and the velocity of the barycenter of the bubble. This generates a ``vof_barycenter_information.dat`` file which contains the position and the velocity of the barycenter of the bubble.

.. code-block:: text

    #---------------------------------------------------
    # Post-processing
    #---------------------------------------------------

    subsection post-processing
      set verbosity                = quiet
      set calculate VOF barycenter = true
    end


---------------------------
Running the simulation
---------------------------

Call the gls_navier_stokes_2d by invoking:  

``mpirun -np 8 gls_navier_stokes_2d rising-bubble.prm``

to run the simulation using eight CPU cores. Feel free to use more.


.. warning:: 
    Make sure to compile lethe in `Release` mode and 
    run in parallel using mpirun. This simulation takes
    :math:`\approx` 10 mins on 8 processes.


-------
Results
-------

The following image shows the shape and dimensions of the bubble after 3 seconds of simulation, and compares it with results of [`2 <https://doi.org/10.1002/fld.2643>`_, `3 <https://doi.org/10.1002/fld.1934>`_ ].

.. image:: images/bubble.png
    :alt: bubble
    :align: center
    :width: 400

A python post-processing code `(rising-bubble.py)` is added to the example folder to post-process the data files generated by the barycenter post-processing.
Run ``python3 ./rising-bubble.py output`` to execute this post-processing code, where ``output`` is the directory that 
contains the simulation results. The results for the barycenter position and velocity of the bubble are compared with the simulations of Zahedi, Kronbichler, and Kreiss [`2 <https://doi.org/10.1002/fld.2643>`_] and Hysing et al. [`3 <https://doi.org/10.1002/fld.1934>`_]. The following images show the results of these comparisons. The agreement between the two simulations is remarkable considering the coarse mesh used within this example.

.. image:: images/ymean-t.png
    :alt: ymean_t
    :align: center
    :width: 500

.. image:: images/bubble-rise-velocity.png
    :alt: bubble_rise_velocity
    :align: center
    :width: 500

Animation of the rising bubble example:

.. raw:: html

    <iframe width="800" height="450" src="https://www.youtube.com/embed/o73WJ36-2zo"  frameborder="0" allowfullscreen></iframe>

-----------
References
-----------
`[1] <https://doi.org/10.1016/0021-9991(92)90240-Y>`_ Brackbill, J.U., Kothe, D.B. and Zemach, C., 1992. A continuum method for modeling surface tension. Journal of computational physics, 100(2), pp.335-354.

`[2] <https://doi.org/10.1002/fld.2643>`_ Zahedi, S., Kronbichler, M. and Kreiss, G., 2012. Spurious currents in finite element based level set methods for two‐phase flow. International Journal for Numerical Methods in Fluids, 69(9), pp.1433-1456.

`[3] <https://doi.org/10.1002/fld.1934>`_ Hysing, S., Turek, S., Kuzmin, D., Parolini, N., Burman, E., Ganesan, S., & Tobiska, L. (2009). Quantitative benchmark computations of two‐dimensional bubble dynamics. International Journal for Numerical Methods in Fluids, 60(11), 1259-1288.