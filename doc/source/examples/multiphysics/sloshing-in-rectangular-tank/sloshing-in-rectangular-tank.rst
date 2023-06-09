================================
Sloshing in a rectangular tank
================================

This example simulates the damping of a small amplitude wave for Reynolds number of (2, 20, 200 and 2000). The problem is inspired by the test case of `Carrica et al. [1] <https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.1279>`_. 


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
- ``examples/multiphysics/sloshing-in-rectangular-tank/sloshing_0002/sloshing-in-rectangular-tank_Re0002.prm``
- ``examples/multiphysics/sloshing-in-rectangular-tank/sloshing_0020/sloshing-in-rectangular-tank_Re0020.prm``
- ``examples/multiphysics/sloshing-in-rectangular-tank/sloshing_0200/sloshing-in-rectangular-tank_Re0200.prm``
- ``examples/multiphysics/sloshing-in-rectangular-tank/sloshing_2000/sloshing-in-rectangular-tank_Re2000.prm``

-----------------------------
Description of the case
-----------------------------

Predicting the dynamics of free surface waves is essential for many industrial applications (e.g. transport of liquified natural gas). Yet, simulating their dynamics is difficult, especially for high Reynolds number values . Indeed, in this case, the amplitude of the waves dampen very slowly. This leads to an oscillatory wave problem which is highly sensitive to the time integration scheme and the coupling between the VOF solver and the Navier-Stokes solver. 

In this problem, we simulate the damping of a small amplitude wave in a rectangular cavity defined from  :math:`(-1,-1)` to :math:`(1,0.1)`. The initial height of the wave :math:`\xi (x)` is given by:

.. math::

  \xi = 1+0.01*sin(\pi(x+0.5))

Four values of the Reynolds number are investigated: 2, 20, 200 and 2000. 

--------------
Parameter file
--------------

The results for this problem are highly sensitive to the accuracy of the time-stepping scheme. For this reason, we use a 2nd order backward differentiation scheme (``bdf2``) with a variable time step. The ``adaptive time step scaling`` is set to 1.025 to ensure that the time-step does not rise too quickly during wave oscillations.

.. code-block:: text

    # --------------------------------------------------
    # Simulation Control
    #---------------------------------------------------
    subsection simulation control
      set method                       = bdf2
      set time end                     = 50
      set time step                    = 0.01
      set adapt                        = true
      set max cfl                      = 0.25
      set output name                  = sloshing-in-rectangular-tank_Re20
      set output path                  = ./output_Re20/
      set output frequency             = 1
      set adaptative time step scaling = 1.025
    end

The ``multiphysics`` subsection is used to enable the VOF solver.

.. code-block:: text

    #---------------------------------------------------
    # Multiphysics
    #---------------------------------------------------
    subsection multiphysics
      set VOF  = true
    end 
    

In the ``initial condition``, we define the initial height of the wave, such that the interface (:math:`\phi = 0.5` isocurve) lies at the right height.

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
        set Function expression =  if (y<=(0.01*sin(3.1416*(x+0.5))), min(0.5-(y-0.01*sin(3.1416*(x+0.5)))/0.0025,1), max(0.5-(y-0.01*sin(3.1416*(x+0.5)))/0.0025,0))
      end
    end

In the ``mesh`` subsection, we define a hyper rectangle with appropriate dimensions. The mesh is initially refined 6 times to ensure adequate definition of the interface.

.. code-block:: text

  #---------------------------------------------------
  # Mesh
  #---------------------------------------------------

  subsection mesh
    set type               = dealii
    set grid type          = subdivided_hyper_rectangle
    set grid arguments     = 5, 2 : -1, -1 : 1, 0.1 : true
    set initial refinement = 6
  end


The ``physical properties`` are mainly used to establish the Reynolds number of the sloshing liquid. For the air, however, the work of Carrica *et al.* `[1]  <https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.1279>`_ does not give any physical properties. We thus fix the air to be significantly less dense than the liquid, but we keep its viscosity at a certain reasonable viscosity to ensure numerical stability.

.. code-block:: text

  #---------------------------------------------------
  # Physical Properties
  #---------------------------------------------------
  
  subsection physical properties
    set number of fluids = 2
    subsection fluid 0
      set density             = 0.001
      set kinematic viscosity = 0.001
    end
    subsection fluid 1
      set density             = 1
      set kinematic viscosity = 0.5
    end
  end


The ``source term`` subsection is used to enable the gravitational acceleration along the :math:`y` direction.

.. code-block:: text

  # --------------------------------------------------
  # Source term
  #---------------------------------------------------
  
  subsection source term
    set enable = true
    subsection xyz
      set Function expression = 0 ; -1 ; 0
    end
  end
    


-----------------------
Running the simulation
-----------------------

We can call the gls_navier_stokes_2d for each Reynolds number. For :math:`Re=20`, this can be done by invoking the following command:

.. code-block:: text

  mpirun -np 8 gls_navier_stokes_2d sloshing-in-rectangular-tank_Re0020.prm

to run the simulation using eight CPU cores. Feel free to use more.


.. warning:: 
    Make sure to compile lethe in `Release` mode and 
    run in parallel using mpirun. This simulation takes
    :math:`\approx` 8 minutes (Re=2) to 6 hours (Re=2000) on 8 processes.


-------
Results
-------

We compare the relative height of the free surface at :math:`x=0` with an analytical solution proposed by `Wu et al. [2] <https://link.springer.com/article/10.1023/A:1017558826258>`_. For the Reynolds number of 2, 20 and 200, data were directly extracted from `[1] <https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.1279>`_, whereas for the Reynolds of 2000, the simplified analytical expression of Wu *et al.* `[2] <https://link.springer.com/article/10.1023/A:1017558826258>`_ is used. The results for Reynolds number of 2, 20, 200 and 2000 can be post-processed by invoking the following command from the folder of the Reynolds number of interest (Re=20 in the example below):

.. code-block:: text

  python3 ../sloshing_post_processing.py . sloshing-in-rectangular-tank_Re0020.prm

.. important::
    You need to ensure that the ``lethe_pyvista_tools`` module included within Lethe is in your Python path.


The following table presents a comparison between the analytical results and the simulation results for all Reynolds numbers mentioned above. A very good agreement is obtained for each of them, demonstrating the accuracy of the VOF solver.

+------+--------------------------------------+
| Re   | Results                              |
+======+======================================+
| 2    | .. image:: images/Re2.png            |
+------+--------------------------------------+
| 20   | .. image:: images/Re20.png           |
+------+--------------------------------------+
| 200  | .. image:: images/Re200.png          |
+------+--------------------------------------+
| 2000 | .. image:: images/Re2000.png         |
+------+--------------------------------------+

-----------
References
-----------
`[1] <https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.1279>`_ Carrica, P. M., Wilson, R. V., & Stern, F. (2007). An unsteady single‐phase level set method for viscous free surface flows. International Journal for Numerical Methods in Fluids, 53(2), 229-256.


`[2] <https://link.springer.com/article/10.1023/A:1017558826258>`_ Wu, G. X., Eatock Taylor, R., & Greaves, D. M. (2001). The effect of viscosity on the transient free-surface waves in a two-dimensional tank. Journal of Engineering Mathematics, 40, 77-90.