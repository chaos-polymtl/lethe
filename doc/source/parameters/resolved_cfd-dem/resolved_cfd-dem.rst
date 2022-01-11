***********************************************
Resolved CFD-DEM
***********************************************

This subsection contains the parameters related to the resolved CFD-DEM around particles using a sharp interface immersed boundary method. This part of the parameter file is used will using the ``gls_sharp_navier_stokes_2d`` or ``gls_sharp_navier_stokes_3d`` . These solvers can also be used to simulate the flow around static particles. In that case, using this solver eliminates the need to define a conformal mesh for the fluid between the particles.

.. code-block:: text

	subsection particles
		set number of particles                     = 1
		set stencil order                           = 2
		set length ratio                            = 4
		set assemble Navier-Stokes inside particles = false
		set calculate force                         = true
		set ib force output file                    = ib_force
		set ib particles pvd file                   = ib_particles_data
		set initial refinement                      = 0
		set refine mesh inside radius factor        = 0.5
		set refine mesh outside radius factor       = 1.5
		set integrate motion                        = false
		set particle nonlinear tolerance            = 1e-6
		set alpha                                   = 1
		set fluid density                           = 1
		subsection gravity
			set Function expression =0;0;0
		end
		subsection particle info 0
			set density    = 1
			subsection position
				set Function expression =0;0;0
			end
			subsection velocity
				set Function expression =0;0;0
			end
		    	subsection omega
		    		set Function expression =0;0;0
		    	end
		    	set inertia    = 1
		    	set pressure x = 0
		    	set pressure y = 0
		    	set pressure z = 0
		    	set radius     = 0.2
	end
	
* The ``number of particles`` Number of particles reprensented by IB.

* The ``stencil order`` parameter control the order of the Lagrange polynomial used to impose the sharp interface immersed boundary condition. The order of the stencil should be equal or higher then the order of the underlying FEM scheme.

* The ``length ratio`` parameter controls the length of the zone used to define the Lagrange polynomial. The length ratio should be kept as small as possible and above 1. A good starting value is two times the average aspect ratio of the elements in the mesh multiplied by the order of the FEM scheme.

* The ``assemble Navier-Stokes inside particles`` parameter determines if the Navier-Stokes equations are solved inside the particles or not. If the Navier Stokes equation is not solved ( the parameter is false ), the solver will solve a Poisson equation for each variable in the problem. This eliminates the need to define a reference value for the pressure. 

* The ``calculate force`` parameter control if the force is evaluated on each particle. 

* The ``ib force output file`` parameter is the file's name where the derivated value associated with each particle is stored. One file will be created for each particle in the simulation.

* The ``ib particles pvd file`` parameter is the file's name that will be created to animate the particles. This file stores all the derivated values calculated for each of the particles. This file is compatible with Paraview 

* The ``initial refinement`` parameter control how many time the refinement zone around each of the particle is applied before the simulation start. Each application of the refinement zone reduce the size of the element by a factor of two.

* The ``refine mesh inside radius factor`` parameter defined the inside radius of the hyper shell that forms the refinement zone around the particle. The radius used is the product between this factor and the particle's radius. 

* The ``refine mesh outside radius factor`` parameter defined the outside radius of the hyper shell that forms the refinement zone around the particle. The radius used is the product between this factor and the particle's radius. 

* The ``integrate motion`` parameter control if the dynamics equations of the particles are calculated. If this parameter is set to false, the dynamics of the particle will not be calculated. The particle's position and velocity will only be defined by the particles' position and velocity function.

* The ``particle nonlinear tolerance`` parameter controls particle dynamics' nonlinear tolerance. The nonlinear solver won't have converged until the residual on the dynamics equations of all the particles is smaller than this threshold.

* The ``alpha`` parameter is the relaxation parameter used when solving the dynamics equation of the particle.

* The ``fluid density`` parameter is the fluid density used in the force calculation of the resolved CFD-DEM. This parameter is redundant with others in the solver and will be removed in the upcoming code modification.

* The subsection ``gravity`` defines the value of the gravity used in the simulation. This gravity can be defined as a function that evolved in time and space. Each component of the ``Function expression`` corresponds respectively to its magnitude in X, Y, and Z.

* The subsection ``particle info 0`` is used to define relevant information that is specific to the particle with id 0. For each particle with the index n, a new subsection name ``particle info n`` should be defined with relevant information.

The following parameter and subsection are all inside the subsection ``particle info 0`` and have to be redefined for all particle separatly.

* The subsection ``position`` defines the initial value of the particle position if the parameter ``integrate motion`` is set to true. Otherwise, it defines the particle's position at all points in time. This position is expressed as a function that can evolve in time. Each component of the ``Function expression`` corresponds to the value of coordinate X, Y, and Z. 

* The subsection ``velocity`` defines the initial value of the particle velocity if the parameter ``integrate motion`` is set to true. Otherwise, it defines the particle's velocity at all points in time. This velocity is expressed as a function that can evolve in time. Each component of the ``Function expression`` corresponds to the value of its component in the X, Y, and Z direction.

* The subsection ``omega`` defines the initial value of the particle rotational velocity if the parameter ``integrate motion`` is set to true. Otherwise, it defines the particle's rotational velocity at all points in time. This rotational velocity is expressed as a function that can evolve in time. Each component of the ``Function expression`` corresponds to the value of its component in the X, Y, and Z direction. It's important to note that even the 2D solver uses the rotational velocity in 3D. In that case, it will only use the Z component of the rotational velocity.

* The ``inertia`` parameter is used to define one of the diagonal elements of the rotational inertia matrix. Since we are defining spherical particles, we assume a uniform distribution of mass, and as such, all the diagonal elements of the rotational inertia matrix are the same.

* The ``pressure x``, ``pressure y``, and ``pressure z`` parameters are used to define the X, Y, and Z coordinate offset of the pressure reference point relative to the center of the particle. These parameters are used when the ``assemble Navier-Stokes inside particles`` parameter is set to true to define the pressure reference point.

* The ``radius`` parameter is used to define the radius of this particle.


