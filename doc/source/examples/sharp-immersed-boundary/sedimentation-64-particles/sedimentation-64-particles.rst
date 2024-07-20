==============================================================================
Sedimentation of 64 Particles
==============================================================================

This example aims to introduce the user to resolved CFD-DEM simulations with a larger number of particles.


.. warning:: 
    * This case is computationally expensive. It can take several days to run on a desktop computer. You can reduce the computational cost by reducing the mesh density or the number of particles.


----------------------------------
Features
----------------------------------

- Solvers: ``lethe-fluid-sharp`` (with Q1-Q1)
- Transient problem
- Displays the capability of the resolved CFD-DEM solver for the flow around multiple particles
- Displays the robustness of the resolved CFD-DEM Solver.

---------------------------
Files Used in This Example
---------------------------

Both files mentioned below are located in the example's folder (``examples/sharp-immersed-boundary-solver/sedimentation-64-particle``).

- Parameter file: ``sedimentation-64-particle.prm``
- Particles file: ``particles.input``


-----------------------
Description of the Case
-----------------------
The case consists of the release of 64 particles (:math:`\rho_p=0.0015 \frac{\text{kg}}{\text{cm}^{3}}`) with a diameter of 0.25 cm arranged in a 4 by 4 by 4 cubic array centered 21 cm above the bottom of the container. The container is a 2 by 2 by 24 cm rectangle. The viscosity of the fluid is :math:`\mu_f=0.0001 \frac{\text{kg}}{\text{s cm}}`. The density of the fluid is :math:`\rho_f=0.001 \frac{\text{kg}}{\text{cm}^{3}}`. The gravity constant is :math:`g= -981 \frac{\text{cm}}{\text{s}^{2}}`. The particles accelerate due to gravity until they hit the bottom of the container, which at this point, we stop the simulation. All the container walls have no-slip boundary conditions except at the top of the container, where we define an open boundary.


---------------
Parameter File
---------------

We explain every part of this parameter file in detail. In each section of the parameter file, we describe relevant parameters. The omitted parameters are only user preference parameters and have no impact on the simulation results. For more information on these parameters, refer to :doc:`../../../parameters/parameters`.
 
Simulation Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: text

    subsection simulation control
      set method             = bdf2
      set bdf startup method = multiple step bdf
      set time step          = 0.0025 
      set time end           = 5      
      set output name        = out    
      set output frequency   = 1      
    end


* The ``method`` is set to  ``bdf2`` to have a second-order time-stepping method. This ensures a low error due to the time discretization in this case.


* The ``time step`` is set to  0.0025. This time step is small enough to prevent large error due to the time discretization. 

* The ``time end`` is set to  4.0. This is slightly longer than the time needed for all the particles to reach the bottom of the container

Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: text

    subsection physical properties
      subsection fluid 0
        set kinematic viscosity = 0.1
        set density             = 0.001
      end
    end


* The ``kinematic viscosity`` is set to  0.1. This value is derived from the case description by dividing :math:`\mu_f` by :math:`\rho_f`. This parameter changes the Reynolds number of the case since it is one of the variables in the evaluation of the Reynolds number but also by changing the terminal velocity of the particle.

* The ``density`` is set to 0.001 according to the description of the problem.

FEM
~~~
.. code-block:: text

    subsection FEM
      set velocity order = 1
      set pressure order = 1
    end

Here we use Q1Q1 elements to reduce the computational cost.

Mesh
~~~~~~
.. code-block:: text

    subsection mesh
        set type                 = dealii
        set grid type            = subdivided_hyper_rectangle
        set grid arguments       = 1,1,12: 0,0,0 : 2 , 2 , 24 : true
        set initial refinement   = 4
    end

The domain is a rectangular box: we can directly use a subdivided hyper rectangle mesh from the deal.II library. In this case, we have oriented the z-direction with gravity. As such, we have the long side of the box along this axis.

* The ``grid arguments`` is set to  ``1,1,12: 0,0,0 : 2,2,24 : true``. This section has 3 subsections. First ``1,1,12`` describes the initial subdivision of the box. This subdivision has been chosen as it is the smallest mesh we can do of the box to have cubic elements. Secondly ``0,0,0 : 2,2,24`` describes the 2 points from which we have derived the rectangular box (0,0,0) and  (2,2,24). Finally, we have ``true``, which is a boolean to activate the coloration of the boundary. This allows us to define separate boundary conditions at each side of the box.

* The ``initial refinement`` is set to 4. This will ensure to have a base mesh that is a bit finer than the particle size.

Mesh Adaptation
~~~~~~~~~~~~~~~
.. code-block:: text

    subsection mesh adaptation
      set fraction coarsening = 0.2
      set fraction refinement = 0.025
      set fraction type = number
      set frequency = 1
      set max number elements = 750000
      set max refinement level = 6
      set min refinement level = 4
      set type = kelly
      set variable = velocity
    end

* The ``fraction coarsening`` is set to 0.2. This limits the accumulation of elements when the particle is moving. It allows for cells far from the particle to be coarsened when the particles get further away.

* The ``fraction refinement`` is set to 0.025. The objective here is to refine elements that become close to the particle when it's moving. This will mostly refine elements around the particle that are not included in the refinement zone around the particle. The refinement zone around the particle will be discussed in more detail in the IB particle section.

* The ``frequency`` is set to 1. Since the particle is moving at each time step, the refinement zone around it should be reevaluated at each time step.

* The ``max refinement level`` is set to 6. This parameter limits how small the elements around the particle can get, limiting the total number of elements in the problem. Here we limit the mesh size to 8 elements per diameter of the particle. This should be sufficient to show the capabilities of the solver. However, the discretization error is not negligible in this case.

* The ``type`` is set to ``kelly``. Since the particle is moving and we do not want a uniform refinement of all the cells, we use the kelly error estimator based on the ``velocity`` variable.

Boundary Conditions
~~~~~~~~~~~~~~~~~~~
.. code-block:: text

    subsection boundary conditions
      set number = 5

      subsection bc 0
        set id   = 0
        set type = noslip
      end
      subsection bc 1
        set id   = 1
        set type = noslip
      end
      subsection bc 2
        set id   = 2
        set type = noslip
      end
      subsection bc 3
        set id   = 3
        set type = noslip
      end
      subsection bc 4
        set id   = 4
        set beta = 10
        set type = function
      subsection u
        set Function expression = 0
      end
      subsection v
        set Function expression = 0
      end
      subsection w
         set Function expression = 0
      end
    end

Here we define the 5 ``no slip`` boundaries for all the box walls and let the 6th boundary free, to represent the top of the box. We refer the reader to the :doc:`../../../parameters/cfd/boundary_conditions_cfd` section on how those boundaries are defined. 

.. note:: 
    The boundary id of dealii rectangular mesh are numbered as such:  :math:`x_{min}=0`, :math:`x_{max}=1`, :math:`y_{min}=2`, :math:`y_{max}=3`, :math:`z_{min}=4`, :math:`z_{max}=5`.


Initial Conditions
~~~~~~~~~~~~~~~~~~
.. code-block:: text

    subsection initial conditions
      set type = nodal
      subsection uvwp
        set Function expression = 0; 0; 0; 0
      end
    end

The initial condition for this case is simple to define. At the start of the simulation, we assume that the particle and the fluid are at rest. From there, we define a uniform velocity field of 0 everywhere. To do that, we used the ``type = nodal`` and then specified a function expression of 0 for all the velocity components.  

Non-linear Solver
~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection non-linear solver
      subsection fluid dynamics
        set verbosity             = verbose
        set tolerance             = 1e-4
        set max iterations        = 10
        set residual precision    = 5
        set force rhs calculation = true
      end
    end

* The ``tolerance`` is set to 1e-4. This is small enough to ensure that the flow field is adequately resolved, since here we expect a velocity of the particle of the order of 10.

* The ``max iterations`` is set to 10. The objective here is to allow enough Newton non-linear steps to ensure the convergence to the tolerance. Also, we should limit the time spent on a single time step if the system is too stiff.  

* The ``force rhs calculation`` is set to ``true``. This is the most important modification for resolved CFD-DEM simulation. By default, the non-linear solver will recalculate the RHS only after the update of the solution. But here, we need to evaluate it before every matrix resolution, and we cannot use the last RHS evaluation that was done after the last newton iteration. The particle position was updated between these two steps, changing the RHS evaluation. This means that for every non-linear step, we evaluate the RHS twice. The non-linear solver follows this sequence of steps for each newton iteration.
    * update the particles positions
    * update the Jacobian matrix
    * update the RHS
    * solve the matrix system
    * reevaluate the RHS to check the convergence.


Linear Solver
~~~~~~~~~~~~~

.. code-block:: text

    subsection linear solver
      subsection fluid dynamics
        set method                                = gmres
        set max iters                             = 1000
        set relative residual                     = 1e-4
        set minimum residual                      = 1e-11
        set preconditioner                        = ilu
        set ilu preconditioner fill               = 0
        set ilu preconditioner absolute tolerance = 1e-6
        set verbosity                             = verbose
        set max krylov vectors                    = 1000
      end
    end


* The ``method`` is set to ``gmres``. This solver is less computationally expensive than the other option, and this case does not require any special preconditioner. This makes the ``gmres`` solver with ``ilu`` preconditioner the best option available.

* The ``max iters`` is set to 1000. This is a lot more steps than how much it should take to solve the system.

* The ``max krylov vectors`` is set to 1000. This is to ensure that we keep the full Arnoldi basis for each new iteration. From experience keeping a maximum of Krylov vector results in a faster resolution for this case than clearing the basis after a lower number of ``gmres`` iterations.

* The ``relative residual`` is set to 1e-4. This is small enough, so we don't under-resolve our matrix and do extra non-linear steps because of it, and at the same, it doesn't require too many ``gmres`` iterations.

* The ``ilu preconditioner fill`` is set to 0. This is the fastest option with the current simulation parameters. In this case, we can use this option without having to do too many ``gmres`` iterations. It requires less computational time to do a few more  ``gmres`` iterations than building the preconditioner and doing fewer ``gmres`` iterations.

* The ``ilu preconditioner absolute tolerance`` is set to 1e-6. This slightly speeds up the first few matrix resolutions. 

IB Particles
~~~~~~~~~~~~~~

.. code-block:: text

    subsection particles
      set assemble Navier-Stokes inside particles = false
      
      subsection extrapolation function
        set length ratio  = 2
        set stencil order = 2
      end
      
      subsection local mesh refinement
        set initial refinement                = 3
        set refine mesh inside radius factor  = 0
        set refine mesh outside radius factor = 2
      end

      subsection DEM
        set DEM coupling frequency            = 1000
        set particle nonlinear tolerance      = 1e-3
        set contact search radius factor      = 1.5
        set enable lubrication force          = true
        set lubrication range max             = 2
        set lubrication range min             = 0.1
        subsection gravity
          set Function expression = 0;0;-981
        end
      end
      
      subsection input file
        set load particles from file = true
        set particles file           = particles.input
      end
    end

In this subsection, we define most of the parameters that are related to the particle.

* The ``stencil order`` is set to 2 since it improves the results in the force evaluation step and does not make the matrix resolution significantly harder.

* The ``refine mesh inside radius factor`` is set to 0. This creates a mesh refinement inside the particle that avoids having hanging nodes in the calculation and helps ensure a small enough mesh around the particle.

* The ``refine mesh outside radius factor`` is set to 2. This creates a mesh refinement around the particle that avoids having hanging nodes in the calculation and helps ensure a small enough mesh around the particle.

* The ``initial refinement`` is set to 3. Here we want to have the mesh as small as possible for the first time step around each of the particles. To achieve this, we refine every element with at least one vertex in the refinement zone around the particle 3 times before the simulation starts. This ensures that all the cells in the refinement zone around the particle are as small as possible.

* The ``integrate motion`` is set to true because we are interested in the dynamic of the particle as it sediments in the rectangular box.

* The ``assemble Navier-Stokes inside particles`` is set to false because we are not interested in the flow inside of the particle.

* The ``length ratio`` has been set to 2. This is small enough so it does not impact the conditioning of the matrix while avoiding interpolation of the immersed boundary stencil in multiple elements.

* The ``contact search radius factor`` is set to 1.5. This parameter is smaller than the default one since the particle motion relative to their size is relatively slow. This enables the use of a smaller search radius which increases the DEM calculation speed.

* The ``particle nonlinear tolerance`` has been set to 1e-3. This is small enough to ensure that the particle dynamics are adequately resolved. We expect a velocity of the particle of the order of 10.

* The ``DEM coupling frequency`` is set to 1000. This is the number of DEM time steps performed per CFD time step. Here 1000 is enough to prevent instability due to particles' contact.

* The ``enable lubrication force`` is set to true since the subgrid lubrication force model is required to capture the lubrication force between the particles when the gap between them is inferior to two times the mesh size.

* The ``lubrication range max`` is set to 2. The subgrid lubrication force model is enabled when the gap between the particles is smaller than two times the mesh size.

* The ``lubrication range min`` is set to 0.1. The subgrid lubrication force model minimal gap considered between the particles is 0.1 times the mesh size.

* The ``load particles from file`` is set to true to enable the particle to be defined using an external file.

* The ``particles file`` is set to ``particles.input``, which is the file where the particles are defined.

* The ``gravity`` ``Function expression`` is set to 0;0;-981 according to the definition of the case. As we choose the long axis of the rectangular box along the Z, we define gravity in this direction. 

.. note:: 
    The number of particles is not defined here since the particles are defined using a file. In this case the number of particles is defined by the number of particles defined in the file.


---------------
Particles File
---------------
The file from which the particles are defined has a header line that goes as follows:

.. code-block:: text

   type; shape_argument; p_x; p_y; p_z; v_x; v_y; v_z; omega_x; omega_y; omega_z; orientation_x; orientation_y; orientation_z; volume ;density; inertia; pressure_x; pressure_y; pressure_z; youngs_modulus; restitution_coefficient; friction_coefficient; poisson_ratio; rolling_friction_coefficient; integrate_motion;


Each line corresponds to a particle and its properties. A space separates each property. For the details on the properties, see the section :doc:`../../../parameters/sharp-immersed-boundary/sharp-immersed-boundary`. Here the particles' Young's moduli are set to 100MPa, the restitution coefficients to 0.9, the Poisson ratios to 0.30, and the friction coefficients to zero.

.. code-block:: text

    type; shape_argument; p_x; p_y; p_z; v_x; v_y; v_z; omega_x; omega_y; omega_z; orientation_x; orientation_y; orientation_z; volume ;density; inertia; pressure_x; pressure_y; pressure_z; youngs_modulus; restitution_coefficient; friction_coefficient; poisson_ratio; rolling_friction_coefficient; integrate_motion;
    sphere; 0.125; 0.25; 0.25; 20.25; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.001953125; 0.0015; 7.6698974609375e-08; 0.0; 0.0; 0.0; 1000000.0; 0.9; 0.0; 0.3; 0.0; 1.0


---------------
Results
---------------
The results are shown in the animation below. We can see the complex motion of the particles and the way they interact with one another. This case demonstrates the stability of the solver for cases with a large number of particle contacts.


.. note:: 
    The results shown in the animation were obtained with a finer mesh and with a finer time-step.

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/Js73OUr08rM" frameborder="0" allowfullscreen></iframe>



