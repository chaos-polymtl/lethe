====================================
Flow past a Backward-Facing Step
====================================

--------
Features
--------

- Solver: ``lethe-fluid`` (with Q1-Q1)
- Steady and pseudo steady-state solution
- Comparison with benchmark solutions
- Mesh refinement and error analysis


----------------------------
Files Used in This Example
----------------------------

All files mentioned below are located in the example's folder (``examples/incompressible-flow/2d-backward-facing-step``).

- Geometry file: ``backward-facing-step.geo``
- Mesh file: ``backward-facing-step.msh``
- Parameter file for the base case (:math:`\mathrm{Re} = 100`): ``Reynolds100_steady.prm``
- Parameter file for the higher-Reynolds case (:math:`\mathrm{Re} = 1000`): ``Reynolds1000_steadybdf.prm``
- Postprocessing Python script for computing the reattachment length: ``reattachment_length.py``
- Postprocessing Python script for computing velocity distributions at outlet: ``velocity_distribution.py``


-----------------------
Description of the Case
-----------------------

In this example, a bidimensional flow goes past a backward-facing step. The flow enters from the left inlet and separates from the bottom wall at the step, and then reattaches to it further downstream at a distance :math:`x_r` from the step.  

.. image:: image/backward-facing-step-description.png

The backward-facing step problem is a classical computational fluid dynamics problem. The fact that it features a non-trivial solution while maintaining simple geometry and boundary conditions makes this problem a good candidate for validation purposes as well as to test the robustness of a given CFD method. First, the basic parameters of the backward-facing step problem will be shown. A solution to two Reynolds numbers (:math:`\mathrm{Re} = 100` to :math:`\mathrm{Re} =1000`) will then be presented and compared to numerical and analytical data.


--------------
Parameter File
--------------

The following subsections show the different parameters used in the simulations. While they all remain more or less the same throughout the various cases, some of them change as the Reynolds number is increased.

Simulation Control
~~~~~~~~~~~~~~~~~~

At :math:`\mathrm{Re} = 100`, the solution is stable enough to be computed as steady state. We do so by setting ``method`` to ``steady``: 

.. code-block:: text

    subsection simulation control
      set method            = steady
      set number mesh adapt = 10
      set output name       = backward_facing_step_output
      set output frequency  = 1
    end
	
A mesh refinement analysis can be done with ``set number mesh adapt = 10``. By starting from a very coarse mesh and by dynamically refining the mesh at least 10 times, asymptotic convergence can be clearly observed.

However, at :math:`\mathrm{Re} = 1000`, convergence can be quite difficult to obtain while doing a steady state simulation. As the Reynolds number increases, the problem becomes progressively stiffer to a point where the ``steady`` solver ultimately fails. With that in mind, the case can be solved as a transient problem until the steady state solution is obtained. This can be achieved with the ``method = steady_bdf`` parameter.

.. code-block:: text

    subsection simulation control
      set method                       = steady_bdf
      set stop tolerance               = 1e-5
      set time step                    = 0.005
      set adapt                        = true
      set max cfl                      = 1e5
      set adaptative time step scaling = 1.2
      set output name                  = backward_facing_step_output
      set output frequency             = 1
    end
  
``stop tolerance``, ``time step``, ``adapt``, ``max cfl`` and ``adaptive time step scaling`` are parameters that control the pseudo steady-state simulation. In this case, choosing ``stop tolerance = 1e-5`` ensures that the simulation reaches steady state while keeping the number of time iterations to a minimum. Moreover, one can notice a very high value for the ``max cfl``; however, since it is used with ``adaptative time step scaling`` (and since *Lethe* is an implicit solver), even a very high value of the CFL does not compromise the results.

Physical Properties
~~~~~~~~~~~~~~~~~~~

In this problem, the Reynolds number is defined as follows: 

.. math::
	Re_{D_h} = \frac{u D_h}{\nu} = \frac{2uh}{\nu}
	
where :math:`h` is the step height, :math:`D_h = 2h` is the characteristic length and :math:`\nu` the kinematic viscosity.

In addition, unit values of :math:`u` and :math:`h` are chosen in the goal of obtaining an adimensional problem.

.. math::
	Re_{D_h} = f(\nu) = \frac{2}{\nu}
	
Consequently, the physical properties are defined as follows : 

.. code-block:: text
	
    subsection physical properties
      set number of fluids = 1
      subsection fluid 0
        set kinematic viscosity = 0.02 # Re_h=2/nu
      end
    end
	
The ``kinematic viscosity`` is the only parameter that changes coherently with :math:`\mathrm{Re}`. For :math:`\mathrm{Re} = 100`, we set the kinematic viscosity to :math:`0.02`, while for :math:`\mathrm{Re} = 1000` we set it to :math:`0.002`.

Mesh
~~~~

.. code-block:: text

    subsection mesh
      set type      = gmsh
      set file name = ../backward-facing-step.msh
    end
	
The mesh features quadrilateral elements as well as unit step and inlet heights (:math:`h_{in}=h=1`). In that direction, the expansion ratio has been set to :math:`\beta=\frac{h_{out}}{h_{in}}=2` throughout the entirety of the simulations. Also, the inlet and outlet must be spaced far enough apart to allow the formation of a fully developed flow. Finally, since a ``gmsh`` mesh file is used, the initial mesh should be as coarse as possible, since these cells cannot be coarsened with the mesh adaptation algorithm.

Mesh Adaptation
~~~~~~~~~~~~~~~

In this example, the mesh adaptation algorithm is based on the Kelly error estimator applied on the velocity variable. This strategy is suitable here, since a fine mesh is required in the vicinity of the step while a coarser mesh is acceptable far way from it.

.. code-block:: text

    subsection mesh adaptation
      set variable            = velocity
      set type                = kelly
      set fraction refinement = 0.2
    end
	
For higher Reynolds numbers with adjoint time stepping, ``frequency = 5`` can be added to the above parameters to obtain a reasonable number of elements throughout the simulation. In this particular case, the mesh would be refined at every fifth time iteration. As an example, the mesh after eight refinement steps for :math:`\mathrm{Re} = 100` looks as follows:

.. image:: image/8th-mesh.png

FEM
~~~

In this example, the interpolation order has been set to one for both velocity and pressure.

.. code-block:: text

    subsection FEM
      set pressure order = 1
      set velocity order = 1
    end

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

As presented in the description of the case (see figure above), three different boundary conditions (or boundary IDs) are necessary to define this particular problem.

.. code-block:: text

    subsection boundary conditions
      set number = 3
      subsection bc 0
        set id   = 0
        set type = noslip
      end
      subsection bc 1
        set id   = 1
        set type = function
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
      subsection bc 2
        set id   = 2
        set type = outlet
      end
    end
	
First, ``subsection bc 0`` sets a Dirichlet boundary condition (or ``noslip``) at each wall where :math:`\mathbf{u}=\mathbf{0}.` The boundary condition at the inlet is a uniform unit flow such that :math:`[u,v,w] = [1,0,0]`. In that case, the parameter ``type = function`` is used in ``subsection bc 1``. With this parameter, :math:`u`, :math:`v` and :math:`w` can be set numerically and independently. The outflow boundary condition is considered a natural boundary condition (also known as the *do nothing* boundary condition) and it is used since we can consider the outlet to be very far from the step. In fact, this condition specifies :math:`p \rightarrow 0` or in other words, that the traction on the fluid equals zero. In *Lethe*, this particular boundary condition is denoted by ``outlet`` and it is specified for the boundary ID :math:`2`.

Non-linear Solver
~~~~~~~~~~~~~~~~~

The ``newton`` non-linear solver is used with a medium ``tolerance``, since convergence can be hard to obtain for high Reynolds number.

.. code-block:: text

    subsection non-linear solver
      subsection fluid dynamics
        set verbosity      = verbose
        set tolerance      = 1e-6
      end
    end

Linear Solver
~~~~~~~~~~~~~

For :math:`\mathrm{Re} = 100`, standard parameters are suitable to achieve convergence:

.. code-block:: text

    subsection linear solver
      subsection fluid dynamics
        set verbosity                             = verbose
        set method                                = gmres
        set max iters                             = 300
        set max krylov vectors                    = 300
        set relative residual                     = 1e-4
        set minimum residual                      = 1e-9
        set preconditioner                        = ilu
        set ilu preconditioner fill               = 2
        set ilu preconditioner absolute tolerance = 1e-12
        set ilu preconditioner relative tolerance = 1.00
      end
    end         
	
For :math:`\mathrm{Re} = 1000`, however, we use an ``amg`` preconditioner with an ILU smoother with ``amg preconditioner ilu fill = 1`` and increase the number of Krylov vectors:

.. code-block:: text

    subsection linear solver
      subsection fluid dynamics
        set verbosity                   = verbose
        set method                      = gmres
        set max iters                   = 500
        set max krylov vectors          = 500
        set relative residual           = 1e-4
        set minimum residual            = 1e-9
        set preconditioner              = amg
        set amg preconditioner ilu fill = 1
      end
    end
	
.. tip::
	It is important to note that the ``minimum residual`` of the linear solver is smaller than the ``tolerance`` of the nonlinear solver. The reader can consult the `Parameters Guide <https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/linear_solver_control.html>`_ for more information.


-----------------------
Running the Simulations
-----------------------

The simulation can be executed using the following command (assuming that the solver's location is in your PATH environment variable and you want to use ``j`` processes for parallel computations):

.. code-block:: text
  :class: copy-button

  mpirun -np j lethe-fluid Reynolds100_steady.prm

For the case where :math:`\textrm{Re}=1000`, replace the name of the parameter file by ``Reynolds1000_steadybdf.prm``.

----------------------
Results and Discussion
----------------------

:math:`\mathrm{Re}=100`
~~~~~~~~~~~~~~~~~~~~~~~

After opening the file ``backward_facing_step_output.pvd`` with Paraview, the following results are observed:

.. image:: image/Reynolds100_profile.png

It is possible to notice a lot of diffusion past the step. This phenomenon is coherent with what is known of the Navier-Stokes equations: the diffusivity term is inversely proportional to the Reynolds number. Most importantly, a small eddy adjacent to the step is clearly observable. It is also visually noticeable that :math:`2.7 \leq x_r \leq 2.9` (:math:`17.7 \leq x \leq 17.9`). With the Python module `PyVista <https://docs.pyvista.org/>`_, raw simulation data can be extracted (from the .vtu files) and this data can be used to compute :math:`x_r` numerically using the following equation:

.. math::
	\left[ \frac{du}{dy} \right]_{y=0} = 0

which can be resolved with a bisection algorithm or with any other appropriate numerical approach. The postprocessing script provided can be used to compute this reattachment length by using the following command: 

.. code-block:: text
  :class: copy-button

  python3 reattachment_length.py -Re 100

The final value of :math:`x_r` is :math:`2.896` with a relative error of :math:`0.8\%`.  The reference value used to compute the error is the one given by Erturk (2008) [#erturk2008]_.


:math:`\mathrm{Re}=1000`
~~~~~~~~~~~~~~~~~~~~~~~~

In a similar way as we did in the last subsection, for :math:`\mathrm{Re} = 1000` the results that can be visualized in Paraview are the following:

.. image:: image/Reynolds1000_profile.png

On the contrary of what we saw in the :math:`\mathrm{Re} = 100` case, it is clear that there is much less diffusion within the flow. This is once more coherent with the theory. The same eddy as mentioned in the previous section is still present but grows as the Reynolds number is increased. Furthermore, a second principal eddy can be seen adjacent to the top wall in the range :math:`x \in [25,37]`. This "oscillating flow" characteristic is expected of a higher Reynolds flow such as this one. Finally, the :math:`x_r` variable is evaluated visually at :math:`x_r \simeq 12.5` (:math:`x \simeq 27.5`). The same Python code as before can be used by setting the Reynolds flag to ``-Re 1000``; we obtain :math:`x_r = 12.602` as a numerical result with a relative error of :math:`3.9\%`.

----------------------
Velocity Distribution
----------------------

To validate the quality of the mesh/geometry as well, it is interesting to compare the obtained velocity distributions with analytical data. The plots are generated using the following command:

.. code-block:: text
  :class: copy-button

  python3 velocity_distribution.py -Re 100

where the Reynolds number is given through the ``-Re`` flag. The figures illustrate the velocity distributions at the outlet (right wall) in comparison to the analytical solution:

.. image:: image/Reynolds100-poiseuille.png
    :width: 49%
.. image:: image/Reynolds1000-poiseuille.png
    :width: 49%

For :math:`\mathrm{Re} = 1000`, an error in the velocity profile is visually noticeable. We can assume that the outlet is not long enough for the flow to be fully developed at its end, meaning that there is still traction on the fluid. Consequently, increasing this length is essential in order to be able to validate cases where :math:`\mathrm{Re} \geq 1000`.

---------------------------
Possibilities for Extension
---------------------------

- **Test the example for other Reynolds numbers**: the parameter file provided for :math:`\mathrm{Re} = 100` should work for all Reynolds numbers below :math:`\mathrm{Re} = 600`, for higher Reynolds numbers use the parameter file provided for :math:`\mathrm{Re} = 1000`. Reference data for other Reynolds numbers can be found in the ``benchmark-data.txt`` file.
- **Validate with a 3D geometry/mesh**: Since experimental data takes into account 3D effects, it would be interesting to compare numerical data to experimental results.
- **Use second order elements for higher Reynolds simulations**: Using second order elements can improve accuracy for more turbulent flows. Also, it can be very powerful in this particular example, since quadratic elements can theoretically interpolate *Poiseuille* flows with genuinely no numerical error. Consequently, the method can yield incredibly precise results while maintaining a very coarse mesh far from the step. 

----------
References
----------

.. [#erturk2008] \E. Erturk, “Numerical solutions of 2-D steady incompressible flow over a backward-facing step, Part I: High Reynolds number solutions,” *Comput. Fluids*, vol. 37, no. 6, pp. 633–655, Jul. 2008, doi: `10.1016/j.compfluid.2007.09.003 <https://doi.org/10.1016/j.compfluid.2007.09.003>`_\.



