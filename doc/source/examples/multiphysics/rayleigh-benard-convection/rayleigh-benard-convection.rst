==========================
Rayleigh-Bénard Convection
==========================

This example simulates `two-dimensional Rayleigh–Benard convection`_ at Rayleigh numbers of :math:`10^4` and :math:`2.5 \times 10^4` . 

.. _two-dimensional Rayleigh–Benard convection: https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/stochastic-bifurcation-analysis-of-rayleighbenard-convection/019773F174C453F84E7EB179CB1C89F1


----------------------------------
Features
----------------------------------
- Solver: ``gls_navier_stokes_2d`` 
- Buoyant force (natural convection)
- Unsteady problem handled by an adaptive BDF1 time-stepping scheme 


------------------------
Location of the example
------------------------
``examples/multiphysics/rayleigh_benard_convection/rayleigh_benard_convection_Ra10k.prm``
``examples/multiphysics/rayleigh_benard_convection/rayleigh_benard_convection_Ra25k.prm``


-----------------------------
Description of the case
-----------------------------

In this example, we evaluate the performance of the ``gls_navier_stokes_2d`` solver in the simulation of the stability of natural convection within a two-dimensional rectangular domain. The following schematic describes the geometry and dimensions of the simulation in the :math:`(x,y)` plane:

.. image:: images/geometry.png
    :alt: Schematic
    :align: center
    :width: 400


The incompressible Navier-Stokes equations with a Boussinesq approximation for the buoyant force are:
    .. math::
        \nabla \cdot {\bf{u}} = 0

    .. math::
        \rho \frac{\partial {\bf{u}}}{\partial t} + \rho ({\bf{u}} \cdot \nabla) {\bf{u}} = -\nabla p + \nabla \cdot {\bf{\tau}} - \rho \beta {\bf{g}} (T - T_0)

where :math:`\beta` and :math:`T_0` denote thermal expansion coefficient and a reference temperature, respectively.

A two-dimensional block of fluid is heated from its bottom wall at :math:`t = 0` s. The temperature of the bottom wall is equal to :math:`T_h=50`, the temperature of the top wall is equal to :math:`T_c=0`, and the left and right walls are insulated. By heating the fluid from the bottom wall, the buoyant force (natural convection) creates vortices inside the fluid. The shape and number of these vortices mainly depend on the Rayleigh number [1, 2]:

    .. math::
        \text{Ra} = \frac{\rho^2 \beta g (T_h - T_c) H^3 c_p}{k \mu}


where :math:`\rho` is the fluid density, :math:`g` is the magnitude of gravitational acceleration, :math:`H` denotes the characteristic length, :math:`k` is the thermal conduction coefficient, and :math:`\mu` is the dynamic viscosity, and :math:`c_p` is the specific thermal capacity.

In this example, we simulate the Rayleigh-Bénard convection problem at two Rayleigh numbers of 10000 and 25000. According to the literature [1, 2], we should see different numbers (2 and 3 vortices at :math:`Ra=10^4` and :math:`2.5 \times 10^4`, respectively) of vortices in the fluid at these two Rayleigh numbers. The gravity magnitude is set to -10 for both simulations. Additionally, :math:`\rho = 100`, :math:`\beta = 0.0002`, :math:`H = 0.25`, :math:`c_p = 100` and :math:`\mu = 1`. Thus, the Rayleigh number is controlled only by the thermal conduction coefficient for this example. In other words, we change the Rayleigh number by changing the thermal conduction coefficient of the fluid.

.. note:: 
    All four boundary conditions are ``noslip``, and an external 
    gravity field of :math:`-10` is applied in the :math:`y` direction.


--------------
Parameter file
--------------

Time integration is handled by a 1st order backward differentiation scheme 
`(bdf1)`, for a :math:`10000` s simulation time with an initial 
time step of :math:`0.01` second.

.. note::   
    This example uses an adaptive time-stepping method, where the 
    time-step is modified during the simulation to keep the maximum value of the CFL condition below a given threshold (0.5 here). Using ``output control = time``, and ``output time = 25`` the simulation results are written every 25 s.

.. note::   
    Note that the heating process is slow, and the velocity magnitudes are small inside the fluid. Hence, we expect large time-steps and a long simulation.

.. code-block:: text

    # --------------------------------------------------
    # Simulation Control
    #---------------------------------------------------
    subsection simulation control
        set method                          = bdf1
        set time end                        = 10000
        set time step                       = 0.01
        set adapt                           = true
        set max cfl                         = 0.5
        set stop tolerance                  = 1e-5
        set adaptative time step scaling    = 1.3
        set output name                     = rayleigh-benard_convection
        set output control                  = time
        set output time                     = 25
        set output path                     = ./output/     
    end


The ``multiphysics`` subsection enables to turn on ``true`` and off ``false`` the physics of interest. Here ``heat transfer``, ``buoyancy force``, and ``fluid dynamics`` are chosen.

.. code-block:: text

    #---------------------------------------------------
    # Multiphysics
    #---------------------------------------------------
    subsection multiphysics
        set buoyancy force      = true
        set heat transfer       = true
        set fluid dynamics      = true
    end 
    
The ``source term`` subsection defines gravitational acceleration.

.. code-block:: text
    
    #---------------------------------------------------
    # Source term
    #---------------------------------------------------
    subsection source term
        set enable                      = true
        subsection xyz
            set Function expression     = 0 ; -10 ; 0
        end
    end


The ``physical properties`` subsection defines the physical properties of the fluid. Since we simulate the Rayleigh-Bénard convection at two Rayleigh numbers (:math:`Ra=10^4` and :math:`2.5 \times 10^4`), we use different thermal conductivities to reach mentioned Rayleigh numbers. We change the thermal conductivity of the fluid in the two simulations. Note that any other physical property (that is present in the Rayleigh number equation defined above) can be used instead of thermal conductivity. Both thermal conductivity values (:math:`k=0.15625` for :math:`Ra=10^4`, and :math:`k=0.0625` for :math:`Ra=2.5 \times 10^4`) are added to the parameter handler file. However, only one of them should be uncommented for each simulation.


.. code-block:: text

    #---------------------------------------------------
    # Physical Properties
    #---------------------------------------------------
    subsection physical properties
        set number of fluids            = 1
        subsection fluid 0
            set density                 = 100
            set kinematic viscosity     = 0.01
            set thermal expansion       = 0.0002
            set thermal conductivity    = 0.15625	# for Ra = 10000
            #set thermal conductivity   = 0.0625	# for Ra = 25000
            set specific heat           = 100
        end
    end

---------------------------
Running the simulation
---------------------------

Call the gls_navier_stokes_2d by invoking:  

``mpirun -np 8 gls_navier_stokes_2d rayleigh_benard_convection_Ra10k.prm``

and

``mpirun -np 8 gls_navier_stokes_2d rayleigh_benard_convection_Ra25k.prm``

to run the simulations using eight CPU cores. Feel free to use more. Note that the first and second commands belong to the simulations at :math:`Ra=10^4` and :math:`Ra=2.5 \times 10^4`, repectively.


.. warning:: 
    Make sure to compile lethe in `Release` mode and 
    run in parallel using mpirun. This simulation takes
    :math:`\approx` 20 minutes on 8 processes.


-------
Results
-------

The following animation shows the results of this simulation:

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/tEg5M-wiCp8" frameborder="0" allowfullscreen></iframe>


Note that at Ra=10000, two vortices exist in the fluid, while an extra (relatively small) vortex appears near the right wall. The velocity magnitude in the vortices is larger at smaller Rayleigh number.

-----------
References
-----------
[1] Venturi, D., Wan, X. and Karniadakis, G.E., 2010. Stochastic bifurcation analysis of Rayleigh–Bénard convection. Journal of fluid mechanics, 650, pp.391-413.

[2] `https://www.mis.mpg.de/applan/research/rayleigh.html`_

.. _https://www.mis.mpg.de/applan/research/rayleigh.html: https://www.mis.mpg.de/applan/research/rayleigh.html