==========================
Rayleigh-Bénard Convection
==========================

This example simulates two-dimensional Rayleigh–Benard convection [#venturi2010]_ [#mpi2022]_ at Rayleigh numbers of :math:`10^4` and :math:`2.5 \times 10^4` .


----------------------------------
Features
----------------------------------

- Solver: ``lethe-fluid`` 
- Buoyant force (natural convection)
- Unsteady problem handled by an adaptive BDF1 time-stepping scheme 


---------------------------
Files Used in This Example
---------------------------

Both files mentioned below are located in the example's folder (``examples/multiphysics/rayleigh-benard-convection``).

- Parameter file for :math:`Ra=10\, 000`: ``rayleigh-benard-convection-Ra10k.prm``
- Parameter file for :math:`Ra=25\, 000`: ``rayleigh-benard-convection-Ra25k.prm``


-----------------------------
Description of the Case
-----------------------------

In this example, we evaluate the performance of the ``lethe-fluid`` solver in the simulation of the stability of natural convection within a two-dimensional rectangular domain. The following schematic describes the geometry and dimensions of the simulation in the :math:`(x,y)` plane:

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
Parameter File
--------------

Simulation Control
~~~~~~~~~~~~~~~~~~

Time integration is handled by a 1st order backward differentiation scheme 
`(bdf1)`, for a :math:`10000` s simulation time with an initial 
time step of :math:`0.01` second.

.. note::   
    This example uses an adaptive time-stepping method, where the 
    time-step is modified during the simulation to keep the maximum value of the CFL condition below a given threshold (0.5 here). Using ``output control = iteration``, and ``output frequency = 100`` the simulation results are written every 100 iteration regardless of the time steps.

.. note::   
    Note that the heating process is slow, and the velocity magnitudes are small inside the fluid. Hence, we expect large time-steps and a long simulation.

.. code-block:: text

    subsection simulation control
      set method                       = bdf1
      set time end                     = 10000
      set time step                    = 0.01
      set adapt                        = true
      set max cfl                      = 0.5
      set stop tolerance               = 1e-5
      set adaptative time step scaling = 1.3
      set number mesh adapt            = 0
      set output name                  = rayleigh-benard_convection
      set output control               = iteration
      set output frequency             = 100
      set output path                  = ./output/
    end

Multiphysics
~~~~~~~~~~~~

The ``multiphysics`` subsection enables to turn on ``true`` and off ``false`` the physics of interest. Here ``heat transfer``, ``buoyancy force``, and ``fluid dynamics`` are chosen.

.. code-block:: text

    subsection multiphysics
      set buoyancy force = true
      set heat transfer  = true
      set fluid dynamics = true
    end

Source Term
~~~~~~~~~~~

The ``source term`` subsection defines gravitational acceleration.

.. code-block:: text
    
    subsection source term
      subsection fluid dynamics
        set Function expression = 0 ; -10 ; 0
      end
    end

Physical Properties
~~~~~~~~~~~~~~~~~~~

The ``physical properties`` subsection defines the physical properties of the fluid. Since we simulate the Rayleigh-Bénard convection at two Rayleigh numbers (:math:`Ra=10^4` and :math:`2.5 \times 10^4`), we use different thermal conductivities to reach mentioned Rayleigh numbers. We change the thermal conductivity of the fluid in the two simulations. Note that any other physical property (that is present in the Rayleigh number equation defined above) can be used instead of thermal conductivity. Both thermal conductivity values (:math:`k=0.15625` for :math:`Ra=10^4`, and :math:`k=0.0625` for :math:`Ra=2.5 \times 10^4`) are added to the parameter handler file. However, only one of them should be uncommented for each simulation.


.. code-block:: text

    subsection physical properties
      set number of fluids = 1
      subsection fluid 0
        set density              = 100
        set kinematic viscosity  = 0.01
        set thermal expansion    = 0.0002
        set thermal conductivity = 0.15625 # for Ra = 10000
        #set thermal conductivity = 0.0625 # for Ra = 25000
        set specific heat        = 100
      end
    end


---------------------------
Running the Simulation
---------------------------

Call the ``lethe-fluid`` by invoking:

.. code-block:: text
  :class: copy-button

  mpirun -np 8 lethe-fluid rayleigh-benard-convection-Ra10k.prm

and

.. code-block:: text
  :class: copy-button

  mpirun -np 8 lethe-fluid rayleigh-benard-convection-Ra25k.prm

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


Note that at :math:`Ra=10000`, two vortices exist in the fluid, while an extra (relatively small) vortex appears near the right wall for :math:`Ra=25000`. The velocity magnitude in the vortices is larger at smaller Rayleigh number.


-----------
References
-----------

.. [#venturi2010] \D. Venturi, X. Wan, and G. E. Karniadakis, “Stochastic bifurcation analysis of Rayleigh–Bénard convection,” *J. Fluid Mech.*, vol. 650, pp. 391–413, May 2010, doi: `10.1017/S0022112009993685 <https://doi.org/10.1017/S0022112009993685>`_\.

.. [#mpi2022] \“Rayleigh-Bénard Convection” *Max Planck Institute*, Accessed: 17 Jul. 2024, Available: https://archive.ph/XrJXx\.
