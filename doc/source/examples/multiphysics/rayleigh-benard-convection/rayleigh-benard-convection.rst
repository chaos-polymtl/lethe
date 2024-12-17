==========================
Rayleigh-Bénard Convection
==========================

This example simulates two-dimensional Rayleigh–Benard convection [#ouertatani][#venturi2010]_ [#mpi2022]_ at Rayleigh numbers of :math:`10^4` and :math:`10^6` .


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

- Parameter file for :math:`Ra=10^4`: ``rayleigh-benard-convection-Ra10k.prm``
- Parameter file for :math:`Ra=10^6`: ``rayleigh-benard-convection-Ra25k.prm``


-----------------------------
Description of the Case
-----------------------------

In this example, we evaluate the performance of the ``lethe-fluid`` solver in the simulation of the stability of natural convection within a two-dimensional square domain. Our results are tested against the benchmark presented by Ouertatani et al. [#ouertatani]. The following schematic describes the geometry and dimensions of the simulation in the :math:`(x,y)` plane:

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

A two-dimensional block of fluid is heated from its bottom wall at :math:`t = 0` s. The temperature of the bottom wall is equal to :math:`T_h=10`, the temperature of the top wall is equal to :math:`T_c=0`, and the left and right walls are insulated. By heating the fluid from the bottom wall, the buoyant force (natural convection) creates vortices inside the fluid. The shape and number of these vortices depend on the Rayleigh number and the Prandtl [1, 2]:

    .. math::
        \text{Ra} = \frac{\rho^2 \beta g (T_h - T_c) H^3 c_p}{k \mu}
        \text{Pr} = \frac{\mu c_p}{k}


where :math:`\rho` is the fluid density, :math:`g` is the magnitude of gravitational acceleration, :math:`H` denotes the characteristic length, :math:`k` is the thermal conduction coefficient, and :math:`\mu` is the dynamic viscosity, and :math:`c_p` is the specific heat capacity.

In this example, we simulate the Rayleigh-Bénard convection problem at two Rayleigh numbers of :math:`Ra=10^4` and :math:`Ra=10^6` with a Prandtl number of :math:`Pr=0.71` which correspond to air. According to the literature [1], we should see one big convective cell at steady-state for both :math:`Ra=10^4` and :math:`Ra=10^6`, but for the latter, there should also be two small voriticies in opposite corners rotating in the reverse direction of the big vortex. The gravity magnitude is set to -10 for both simulations for simplicity. Additionally, because the two adimensionale number above are the only thing that characterize the flow we may choose the remaining parameter as we want. Here we chose for simplicity :math:`\rho = 1`, :math:`\beta = 0.0002`, :math:`H = 1`, :math:`c_p = 100` and :math:`\mu = 0.071`. Thus, the Rayleigh number is controlled only by the thermal expension coefficient (:math:`\beta = 0.71` or :math:`\beta = 71`) for this example. In other words, we change the Rayleigh number by changing the thermal expansion coefficient of the fluid.

.. note:: 
    All four boundary conditions are ``noslip``, and an external 
    gravity field of :math:`-10` is applied in the :math:`y` direction.


--------------
Parameter File
--------------

Simulation Control
~~~~~~~~~~~~~~~~~~

Time integration is handled by a 1st order backward differentiation scheme 
`(bdf1)`, for a :math:`24`s and :math:`10`s simulation time with an initial 
time step of :math:`0.01` second. Those are obtained by simulating for longer time periodes using coarser grid.

.. note::   
    This example uses an adaptive time-stepping method, where the 
    time-step is modified during the simulation to keep the maximum value of the CFL condition below a given threshold (0.5 here). Using ``output control = time``, and ``output time frequency = 0.5`` the simulation results are written every 0.5 seconds regardless of the time steps. This is chosen 

.. note::   
    Note that at the resolution to match the article (256x256) result in really small time step which result in a rather long simulation. If desired, you can choose to reduced to 7 or 6 ``initial refinement level`` of the prm file to reduce the simulation time without loosing too much accuracy.

.. code-block:: text

    subsection simulation control
      set method                       = bdf1
      set time end                     = 24
      set time step                    = 0.01
      set adapt                        = true
      set max cfl                      = 0.5
      set stop tolerance               = 1e-5
      set adaptative time step scaling = 1.3
      set number mesh adapt            = 0
      set output name                  = rayleigh-benard_convection
      set output control               = time
      set output time frequency        = 0.5
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

The ``physical properties`` subsection defines the physical properties of the fluid. Since we simulate the Rayleigh-Bénard convection at two Rayleigh numbers (:math:`Ra=10^4` and :math:`Ra=10^6`), we use different thermal expension coefficient to reach mentioned Rayleigh numbers. Thermal expension coefficient values (:math:`k=0.71` for :math:`Ra=10^4`, and :math:`k=71` for :math:`Ra=10^6`) are adjusted in the respective parameter handler file. Here we present only the parameter file for :math:`Ra=10^4`.


.. code-block:: text

    subsection physical properties
      set number of fluids = 1
      subsection fluid 0
        set density              = 1
        set kinematic viscosity  = 0.071
        set thermal expansion    = 0.71
        set thermal conductivity = 10
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

  mpirun -np 8 lethe-fluid rayleigh-benard-convection-Ra1M.prm

to run the simulations using eight CPU cores. Feel free to use more. Note that the first and second commands belong to the simulations at :math:`Ra=10^4` and :math:`Ra=10^6`, repectively.


.. warning:: 
    Make sure to compile lethe in `Release` mode and 
    run in parallel using mpirun. The first simulation takes
    :math:`\approx` 20 minutes on 8 processes and the second :math:`\approx` 2 days.


-------
Results
-------

The following animation shows the results of the simulation at :math:`Ra=10^6`:

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/tEg5M-wiCp8" frameborder="0" allowfullscreen></iframe>



-----------
References
-----------

.. [#ouertatani] \N. Ouertatani, N. Ben Cheikh, B. Ben Beya, T. Lili, "Numerical simulation of two-dimensional Rayleigh-Bénard convection in an enclosure," Comptes Rendus – Mec. 2008;336(5):464–70. `10.1016/j.crme.2008.02.004 <https://comptes-rendus.academie-sciences.fr/mecanique/articles/10.1016/j.crme.2008.02.004/>`_\.

.. [#venturi2010] \D. Venturi, X. Wan, and G. E. Karniadakis, “Stochastic bifurcation analysis of Rayleigh–Bénard convection,” *J. Fluid Mech.*, vol. 650, pp. 391–413, May 2010, doi: `10.1017/S0022112009993685 <https://doi.org/10.1017/S0022112009993685>`_\.

.. [#mpi2022] \“Rayleigh-Bénard Convection” *Max Planck Institute*, Accessed: 17 Jul. 2024, Available: https://archive.ph/XrJXx\.
