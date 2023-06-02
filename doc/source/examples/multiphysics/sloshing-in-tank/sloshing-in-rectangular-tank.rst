================================
Sloshingin a rectangular tank
================================

This example simulates the damping of a small amplitude wave for Reynolds number of 2, 20, 200 and 2000. The problem is inspired by the test case of s


----------------------------------
Features
----------------------------------
- Solver: ``gls_navier_stokes_2d`` 
- Volume of fluid (VOF)
- Unsteady problem handled by an adaptive BDF2 time-stepping scheme 
- Usage of a python script for post-processing data


---------------------------
Files used in this example
---------------------------
``examples/multiphysics/sloshing-in-rectangular-tank/sloshing_0002/sloshing-in-rectangular-tank_Re0002.prm``
``examples/multiphysics/sloshing-in-rectangular-tank/sloshing_0020/sloshing-in-rectangular-tank_Re0020.prm``
``examples/multiphysics/sloshing-in-rectangular-tank/sloshing_0200/sloshing-in-rectangular-tank_Re0200.prm``
``examples/multiphysics/sloshing-in-rectangular-tank/sloshing_2000/sloshing-in-rectangular-tank_Re2000.prm``
``examples/multiphysics/sloshing-in-rectangular-tank/sloshing_post_processing.py``

-----------------------------
Description of the case
-----------------------------

Predicting the dynamics of free surface waves is essential for many industrial applications (e.g. transport of liquified natural gas). Yet, simulating their dynamics is difficult, especially for high values of the Reynolds number. Indeed, in this case, the amplitude of the waves dampen very slowly, which lead to problem with complex dynamics which are highly sensitive to the time integration scheme and the coupling between the Volume of Fluid solver and the Navier-Stokes solver. 

In this problem, we simulate the damping of a small amplitude wave in a cavity of dimension :math:`(-1,-1)X(1,1)` . The initial height of the wave :math:`\xi (x)` is given by:

.. math::

  \xi = 1+0.01*sin(\pi(x+0.5))

Four values of the Reynolds number are investigated: 2, 20, 200 and 2000. 

--------------
Parameter file
--------------

The results for this problem are highly sensitive to the accuracy of the time-stepping scheme. For this reason, we use a 2nd order backward differentiation scheme 
with a variable time step. The

.. code-block:: text

    # --------------------------------------------------
    # Simulation Control
    #---------------------------------------------------
    subsection simulation control
      set method                       = bdf2
      set time end                     = 40000
      set time step                    = 0.1
      set max cfl                      = 0.5
      set adaptative time step scaling = 1.3
      set output name                  = melting
      set output control               = time
      set output time                  = 100
      set output path                  = ./output/      
    end


The ``multiphysics`` subsection is used to enable the VOF solver.

.. code-block:: text

    #---------------------------------------------------
    # Multiphysics
    #---------------------------------------------------
    subsection multiphysics
      set heat transfer  = true
      set buoyancy force = true
      set fluid dynamics = true
    end 
    

In the ``initial condition``, we define the initial height of the wave. We define 

.. code-block:: text

    #---------------------------------------------------
    # Initial condition
    #---------------------------------------------------
    subsection initial conditions
      set type = nodal
      subsection uvwp
        set Function expression = 0; 0; 0
      end
      subsection temperature
        set Function expression = 104.9
      end
    end



---------------------------
Running the simulation
---------------------------

Call the gls_navier_stokes_2d by invoking:  

``mpirun -np 8 gls_navier_stokes_2d melting-cavity.prm``

to run the simulation using eight CPU cores. Feel free to use more.


.. warning:: 
    Make sure to compile lethe in `Release` mode and 
    run in parallel using mpirun. This simulation takes
    :math:`\approx` 0.5 to 6 hours on 8 processes depending on the Reynolds number. 


-------
Results
-------


-----------
References
-----------
`[1] <https://doi.org/10.1016/j.compfluid.2018.03.037>`_ Blais, B. and Ilinca, F., 2018. Development and validation of a stabilized immersed boundary CFD model for freezing and melting with natural convection. Computers & Fluids, 172, pp.564-581.

`[2] <https://doi.org/10.1115/1.3246884>`_ Gau, C. and Viskanta, R., 1986. Melting and solidification of a pure metal on a vertical wall.
