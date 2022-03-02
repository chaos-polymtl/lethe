***********************************************
Resolved CFD-DEM
***********************************************

This subsection contains the parameters related to the resolved CFD-DEM around particles using a sharp interface immersed boundary method. This part of the parameter file concerns the usage of ``gls_sharp_navier_stokes_2d`` or ``gls_sharp_navier_stokes_3d``. These solvers can also be used to simulate the flow around static particles. In that case, using this solver eliminates the need to define a conformal mesh for the fluid between the particles.

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
		set DEM coupling frequency                  = 1000
		set alpha                                   = 1
		set fluid density                           = 1

		
		subsection gravity
			set Function expression =0;0;0
		end
		
		set wall friction coefficient               = 0
		set wall poisson ratio                      = 0.3
		set wall restitution coefficient            = 1
		set wall rolling friction coefficient       = 0
		set wall youngs modulus                     = 100000000
		set lubrication range max		    = 2
		set lubrication range min		    = 0.01
		
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
		    	set friction coefficient         = 0
		    	set poisson ratio                = 0.3
		    	set restitution coefficient      = 1
		    	set rolling friction coefficient = 0
		    	set youngs modulus               = 100000000
		end
	end
	
* The ``number of particles`` is the number of particles simulated by the sharp-edge IB.

* The ``stencil order`` parameter controls the order of the Lagrange polynomial used to impose the sharp interface immersed boundary condition. The order of the stencil should be higher than or equal to the order of the underlying FEM scheme.

* The ``length ratio`` parameter controls the length of the zone used to define the Lagrange polynomial. The length ratio should be kept as small as possible and above 1. A good starting value is twice the average aspect ratio of the elements in the mesh multiplied by the order of the underlying FEM scheme. For example, for Q1 elements with an average aspect ratio of one, the length ratio should be set to 2.

* The ``assemble Navier-Stokes inside particles`` parameter determines if the Navier-Stokes equations are solved inside the particles or not. If the Navier-Stokes equations are not solved (the parameter is false), the solver will solve a Poisson equation for each variable in the problem. This eliminates the need to define a reference value for the pressure. 

* The ``calculate force`` parameter controls if the force is evaluated on each particle. 

* The ``ib force output file`` parameter is the file name where the variables associated with each particle are stored. One file will be created for each particle in the simulation.

* The ``ib particles pvd file`` parameter is the file's name that will be created to animate the particles. This file stores all the variables calculated for each of the particles. This file is compatible with Paraview.

* The ``initial refinement`` parameter controls how many time the refinement zone around each of the particle is applied before the simulation starts. Each application of the refinement zone reduces the size of the elements by a factor two.

* The ``refine mesh inside radius factor`` parameter defines the inside radius of the hyper shell that forms the refinement zone around the particles. The radius used is the product between this factor and the particle's radius. 

* The ``refine mesh outside radius factor`` parameter defines the outside radius of the hyper shell that forms the refinement zone around the particles. The radius used is the product between this factor and the particle's radius. 

* The ``integrate motion`` parameter controls if the dynamics equations of the particles are calculated. If this parameter is set to false, the particles remain static.  If ``ìntegrate motion=true`` the position and the velocity will be defined by the particles' position and velocity function.

* The ``DEM coupling frequency`` parameter controls the number of iterations done on the DEM side for each CFD time step. It's necessary to use a much smaller time step for the particle dynamics than for the fluid in case of contact between the particles. The particle collision happens at a much smaller time-scale than the fluid dynamics.

* The ``particle nonlinear tolerance`` parameter controls particle dynamics' nonlinear tolerance. The nonlinear solver won't have converged until the residual on the dynamics equations of all the particles is smaller than this threshold.

* The ``alpha`` parameter is the relaxation parameter used when solving the dynamics equation of the particle.


* The ``fluid density`` parameter is the fluid density used in the force calculation of the resolved CFD-DEM. This parameter is redundant with others in the solver and will be removed in the upcoming code modification.

* The subsection ``gravity`` defines the value of the gravity used in the simulation. This gravity can be defined as a function that evolves in time and space. Each component of the ``Function expression`` corresponds respectively to its magnitude in X, Y, and Z.

The following properties are used if the particle impact one of the boundaries of the domain. The effective properties used for calculating the impact force are calculated using a harmonic mean of the properties of the wall and the particle.

* The ``wall friction coefficient`` parameter is the coefficient of friction of the wall. This parameter is used to define the effective coefficient of friction between the wall and the particles. At This point in time, all the walls have the same properties.

* The ``wall poisson ratio`` parameter is the Poisson's ratio of the wall's material. This parameter is used to define the nonlinear spring constant used when a particle impacts a wall. At This point in time, all the walls have the same properties.

* The ``wall restitution coefficient`` parameter is the restitution coefficient of the wall's material. This parameter is used to define the effective restitution coefficient for the impact of a particle and the wall. At This point in time, all the walls have the same properties.

* The ``wall rolling friction coefficient`` parameter is the rolling friction coefficient of the wall. This parameter is used to define the effective rolling friction coefficient between the wall and the particles. At This point in time, all the walls have the same properties.

* The ``wall youngs modulus`` parameter is the Young's modulus of the wall's material. This parameter is used to define the nonlinear spring constant used when a particle impacts a wall. At This point in time, all the walls have the same properties.

* The ``lubrication range max`` parameter defines the distance below which the lubrication force between 2 particles or between a particle and a wall is calculated. The range is defined as a multiple of the smallest cell. The lubrication force model is used to model the force between particles when they are too close to each other to model the fluid between them. The distance is expressed as a multiple of the smallest cells. To deactivate this force, simply put this parameter to 0.

* The ``lubrication range min`` parameter defines the minimal distance used in the lubrication force calculation. The range is defined as a multiple of the smallest cell. This limits the force that can be applied on a particle since the lubrification force has a singularity when the distance between 2 particles is 0. We use this parameter to define a lower bound on the distance between 2 particles for the force calculation to avoid this singularity. Physically this distance can be seen as the surface roughness of the particles.

The following parameter and subsection are all inside the subsection ``particle info 0`` and have to be redefined for all particles separatly.

* The subsection ``particle info 0`` is used to define relevant information that is specific to the particle with id 0. For each particle with the index ``n``, a new subsection name ``particle info n`` should be defined with relevant information.



* The subsection ``position`` defines the initial value of the particle position if the parameter ``integrate motion=true``. Otherwise, it defines the particle's position at all points in time. This position is expressed as a function that can evolve in time. Each component of the ``Function expression`` corresponds to the value of coordinate X, Y, and Z. 

* The subsection ``velocity`` defines the initial value of the particle velocity if the parameter ``integrate motion=true``. Otherwise, it defines the particle's velocity at all points in time. This velocity is expressed as a function that can evolve in time. Each component of the ``Function expression`` corresponds to the value of its component in the X, Y, and Z direction.

* The subsection ``omega`` defines the initial value of the particle rotational velocity if the parameter ``integrate motion=true``. Otherwise, it defines the particle's rotational velocity at all times. This rotational velocity is expressed as a function that can evolve in time. Each component of the ``Function expression`` corresponds to the value of its component in the X, Y, and Z direction. It's important to note that even the 2D solver uses the rotational velocity in 3D. In that case, it will only use the Z component of the rotational velocity.

* The ``inertia`` parameter is used to define one of the diagonal elements of the rotational inertia matrix. Since we are defining spherical particles, we assume a uniform distribution of mass, and as such, all the diagonal elements of the rotational inertia matrix are the same.

* The ``pressure x``, ``pressure y``, and ``pressure z`` parameters are used to define the X, Y, and Z coordinate offset of the pressure reference point relative to the center of the particle. These parameters are used when the ``assemble Navier-Stokes inside particles`` parameter is set to true to define the pressure reference point.

* The ``radius`` parameter is used to define the radius of this particle.

The following properties are used if the particle impact one of the boundaries of the domain or another particle. The effective properties used to calculate the impact force are calculated using a harmonic mean of the properties of the particle and the object it impacts.

* The ``friction coefficient`` parameter is the coefficient of friction of the particle. This parameter is used to define the effective coefficient of friction between the wall and the particles.

* The ``poisson ratio`` parameter is the Poisson's ratio of the particle's material. This parameter is used to define the nonlinear spring constant used when a particle impacts a wall.

* The ``restitution coefficient`` parameter is the restitution coefficient of the particles' material. This parameter is used to define the effective restitution coefficient for the impact of a particle and the wall.

* The ``rolling friction coefficient`` parameter is the rolling friction coefficient of the particle. This parameter is used to define the effective rolling friction coefficient between the wall and the particles. The effective coefficient is calculated using a harmonic mean of the properties of the particles and the other objects it impacts.

* The ``youngs modulus`` parameter is the Young's modulus of the particle's material. This parameter is used to define the nonlinear spring constant used when a particle impacts a wall.



