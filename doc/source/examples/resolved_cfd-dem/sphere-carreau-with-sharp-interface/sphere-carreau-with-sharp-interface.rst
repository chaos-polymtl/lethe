================================
Non-Newtonian flow past a sphere
================================

This example showcases a laminar non-Newtonian flow around a sphere, with an *a priori* Reynolds number Re = 50, using the Carreau rheological model.

Features
----------------------------------
- Solvers: ``gls_sharp_navier_stokes_3d`` (with Q1-Q1) 
- Steady-state problem
- Non-Newtonian behavior
- Ramping initial condition
- Displays the use of non-uniform mesh adaptation 

Location of the example
------------------------

- Parameter file: ``/examples/resolved_cfd-dem/sphere_carreau_with_sharp_inferface/sphere_carreau_with_sharp_inferface.prm``


Description of the case
-----------------------

In this example, we study the flow around a static sphere using the sharp-interface method to represent the sphere. The geometry of the flow is the following, with a particle of diameter D = 1.0 located at (0,0,0)
and the flow domain located between (-18,-15,-15) and (42,15,15).

.. image:: images/sharp_carreau_case.png
    :alt: Simulation schematic
    :align: center

Parameter file
-----------------------

Mesh
~~~~~

The mesh is defined using the following subsection.

.. code-block:: text

	subsection mesh
	  set type                 = dealii
	  set grid type            = subdivided_hyper_rectangle
	  set grid arguments       = 2,1,1 : -18,-15,-15 : 42,15,15 : true
	  set initial refinement   = 4
	end
	
Using an ``initial refinement`` of 4, the initial size of the cubic cells is 1.875. Since the particle size is small in regards to the mesh size, a refinement zone is generated around the particle to better capture it .

.. code-block:: text

	subsection  box refinement
	  set initial refinement   = 3
	  subsection mesh
	  set type                 = dealii
	  set grid type            = subdivided_hyper_rectangle
	  set grid arguments       = 1,1,1: -2,-2,-2 : 6,2,2 : true
	  set initial refinement   = 0
	  end
	end

Boundary conditions
~~~~~~~~~~~~~~~~~~~~
We define the boundary conditions in order to have an inlet velocity of 1 m/s on the left, ``slip`` boundary conditions parallel to the flow direction, and an outlet on the right of the domain.

.. code-block:: text

	subsection boundary conditions
	  set number                  = 5
	  subsection bc 0
		set id 		= 0
		set type    = function
		subsection u
		  set Function expression = 1
		end
		subsection v
		  set Function expression = 0
		end
		subsection w
		  set Function expression = 0
		end
	  end
	  subsection bc 1
		set id 		= 2
		set type    = slip
	  end    
	  subsection bc 2
		set id 		= 3
		set type    = slip
	  end
	  subsection bc 3
		set id 		= 4
		set type    = slip
	  end
	  subsection bc 4
		set id 		= 5
		set type    = slip
	end
	end

Physical properties
~~~~~~~~~~~~~~~~~~~~

This example showcases a shear-thinning flow, for which the viscosity decreases when the local shear rate increases. The Carreau model is being used. For more information on rheological models, see :doc:`../../../parameters/cfd/physical_properties`

.. code-block:: text

	subsection physical properties
	  set number of fluids = 1
	  subsection fluid 0
		set rheological model	= carreau
		subsection non newtonian
		  subsection carreau
			set n 		   		= 0.5
			set viscosity_0    	= 0.063403
			set viscosity_inf  	= 0
			set lambda	   		= 10
			set a	           	= 2.0
		  end
		end
	  end
	end

With ``viscosity_inf = 0`` (3-parameter Carreau model), the *a priori* Reynolds number can be estimated using :

.. math::

	 Re = \frac{u_{\infty}D(1+(\lambda(\frac{u_\infty}{D}))^2)^{\frac{1-n}{2}}}{\eta_0}

We use an *a priori* Reynolds number, since it is not possible, *a priori*, to know the effective viscosity of the flow. For the given parameters, the *a priori* Reynolds number is 50. 

Initial conditions
~~~~~~~~~~~~~~~~~~~~

This examples uses a ramping initial condition that first ramps on the ``n`` parameter, and after on the ``viscosity_0`` parameter. This allows for a smooth transition of regime and of non-Newtonian level.

.. code-block:: text

	subsection initial conditions
	  set type = ramp
	  subsection ramp
		subsection n
		  set initial n = 1.0
		  set iterations = 2
		  set alpha = 0.5
		end
		subsection viscosity
		  set initial viscosity = 1.0
		  set iterations = 2
		  set alpha = 0.5
		end
	  end
	end
	
The first initial condition simulation solves for ``n=1.0``, ``viscosity_0 = 1.0``, ``viscosity_inf = 0``, ``lambda=10`` and ``a=2``. The subsequent initial simulations are:

* (Second ``n`` iteration) ``n=0.75``, ``viscosity_0 = 1.0``, ``viscosity_inf = 0``, ``lambda=10`` and ``a=2`` ;
* (First ``viscosity`` iteration) ``n=0.5``, ``viscosity_0 = 1.0``, ``viscosity_inf = 0``, ``lambda=10`` and ``a=2`` ;
* (Second ``viscosity`` iteration) ``n=0.5``, ``viscosity_0 = 0.531702``, ``viscosity_inf = 0``, ``lambda=10`` and ``a=2`` 

and the first simulation uses the parameters in the Physical Properties section. 

Particle
~~~~~~~~~~~~~~~~~~~~

In this case, we want to define a spherical boundary of radius 0.5 center at (0,0,0) that has no velocity. For more information on particle immersed boundary conditions using shar interface, see :doc:`../../../parameters/resolved_cfd-dem/resolved_cfd-dem`.

.. code-block:: text

	subsection particles
	  set number of particles = 1
	  set stencil order = 2
	  set length ratio  = 1
	  set refine mesh inside radius factor = 0.85
	  set refine mesh outside radius factor = 1.3
	  set initial refinement = 2
	  set integrate motion = false
	  set assemble Navier-Stokes inside particles = false    
	  subsection particle info 0
	  subsection position
		set Function expression =0;0;0
	  end
	  subsection velocity
		set Function expression =0;0;0
	  end
		subsection omega
		set Function expression =0;0;0
	  end
		set pressure x =0.00001
		set pressure y =0.00001
		set pressure z =0.00001
		set radius = 0.5
	  end
	end

The hypershell around the boundary between ``refine mesh inside radius factor`` (r = 0.425) and ``refine mesh outside radius factor`` (r = 0.65) will initialy be refined twice. 

Simulation control
~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation is solved at steady-state with 2 mesh adaptation.

.. code-block:: text

	subsection simulation control
  	  set method                  = steady
	  set number mesh adapt       = 2
	  set output name             = sharp-carreau-output
	  set output frequency        = 1
	  set subdivision             = 1
	end

Mesh Adaptation Control
~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to the hypershell refinement zone around the immersed boundary, the ``mesh adaptation`` ``type`` must be set to ``kelly``. During both of the mesh refinement steps, 40% of the cells with be slip in 8 equal cubes (``fraction refinement   = 0.4``) using a velocity-gradient kelly operator.

.. code-block:: text

	subsection mesh adaptation
	  set type                  = kelly
	  set fraction coarsening   = 0.0
	  set fraction refinement   = 0.4
	  set fraction type	      	= number
	  set frequency             = 1
	  set max number elements   = 8000000
	  set min refinement level  = 0
	  set max refinement level  = 11
	  set variable		      	= velocity
	end

Results
---------------
The simulation of this case results in the following solution for the velocity and pressure field. 




We get the following force applied on the particle for each of the mesh refinements, which is similar to the one obtained with a conformal mesh in :doc:`../../incompressible-flow/2d-flow-around-cylinder/2d-flow-around-cylinder`. With the conformal mesh drag force applied to the particle is 7.123. The difference between the 2 can mostly be attributed to the discretization error.

.. code-block:: text

    particle_ID    T_z      f_x       f_y    
          0 -0.033177 5.698080  0.016542 
          0 -0.006670 6.438133  0.004265 
          0 -0.000349 6.773126 -0.000063 
          0  0.000040 6.905268 -0.000170 
          0 -0.000014 6.962307  0.000057 
          
.. note:: 
	The drag coefficient obtained in this case is higher than the drag coefficient for a cylinder at a Reynolds number of 1 as the size of the domain is not large enough relative to the diameter of the cylinder. The flow around the cylinder is then constrained by the lateral boundaries, and this incrases the drag coefficient.
	
	
	
