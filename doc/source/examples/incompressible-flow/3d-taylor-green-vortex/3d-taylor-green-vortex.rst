==================================
Taylor-Couette Flow
==================================

This example showcases another canonicall fluid mechanics problem, the Taylor-Green vortex.  This examples features both the traditional matrix-based solver within Lethe ( ``lethe-fluid``) and the matrix-free solver  (``lethe-fluid-matrix-free``) which is more computationnaly efficient, especially for high-order elements (Q2 and above). Post-processing capabilities for enstrophy and kinetic energy are also demonstrated.


---------
Features
---------

- Solvers: ``lethe-fluid`` (with Q3-Q3) or  ``lethe-fluid-matrix-free`` (with Q3-Q3)
- Transient problem using ``bdf3``time integrator
- Displays the calculation of enstrophy and total kinetic energy


----------------------------
Files Used in This Example
----------------------------

All files mentioned below are located in the example's folder (``examples/incompressible-flow/3d-taylor-green-vortex``).

- Parameter file: ``tgv-m.prm``
- Postprocessing Python script: ``postprocess_taylor_couette.py``


Description of the Case
-----------------------

The Taylor–Green vortex is an unsteady flow of a decaying vortex, which has an exact closed form solution of the incompressible Navier–Stokes equations in Cartesian coordinates. It is named after the British physicist and mathematician Geoffrey Ingram Taylor and his collaborator A. E. Green [1]. In the present case, we simulate one Taylor-Green vortex at a Reynolds number of 1600 in a domain :math:`\Omega = [-\pi,\pi]\times[-\pi,\pi]\times[-\pi,\pi] ` using periodic boundary conditions.

The three velocity components :math:`[u_x,u_y,u_z]^T` and the pressure :math:`p` are specified at time :math:`t=0` by:

.. math::

  u_{x} = sin(x)*cos(y)*cos(z) \\
  u_{y} = -cos(x)*sin(y)*cos(z)\\
  u_{z} = 0 \\
  p =  \frac{1}{16}*(cos(2x)+cos(2y))*(cos(2z)+2)

In this case, the vortex, which is initially 2D, will decay by generating smaller 3D turbulent structures (vortex tubes, rings and sheets). This decay can be monitored through the total kinetic energy of the system. Since the simulation domain is periodic, it can be demontrated that the time derative of the total kinetic energy :math:`E_k` is directly related to the enstrophy :math:`\mathcal{E}` such that:


with
.. math::

  \frac{\mathrm{d}E_K}{\mathrm{d}t} =  \mathcal{E} \\
  E_k = \frac{1}{\Omega} \int_{\Omega} \frac{\mathbf{u}\cdot \mathbf{u}}{2} \mathrm{d}\Omega \\
  \mathcal{E} = \frac{1}{\Omega} \int_{\Omega} \frac{\mathbf{\omega}\cdot \mathbf{\omega}}{2} \mathrm{d}\Omega

where :math:`\mathbf{\omega}=\nabla \times \mathbf{u}` is the vorticity.


--------------
Parameter File
--------------

Mesh
~~~~~

The ``mesh`` subsection specifies the computational grid:

.. code-block:: text

  subsection mesh
    set type               = dealii
    set grid type          = subdivided_hyper_rectangle
    set grid arguments     = 1, 1, 1: -3.14159265359, -3.14159265359, -3.14159265359 : 3.14159265359, 3.14159265359, 3.14159265359 : true
    set initial refinement = 5 # initial mesh refinement
  end

The ``type`` specifies the mesh format used. We use the ``subdivided_hyper_rectangle`` mesh generated from the deal.II `GridGenerator <https://www.dealii.org/current/doxygen/deal.II/namespaceGridGenerator.html>`_ . We set ``colorize=true``to be able to adequately set-up the periodic bondary conditions.


The last parameter specifies the ``initial refinement`` of the grid. Indicating an ``initial refinement=5`` implies that the initial mesh is refined 5 times. In 3D, each cell is divided by 8 per refinement. Consequently, the final grid is made of 32768 cells.

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

The ``boundary conditions`` subsection establishes the constraints on different parts of the domain:

.. code-block:: text

  subsection boundary conditions
    set number = 3
    subsection bc 0
      set type               = periodic
      set id                 = 0
      set periodic_id        = 1
      set periodic_direction = 0
    end
    subsection bc 1
      set type               = periodic
      set id                 = 2
      set periodic_id        = 3
      set periodic_direction = 1
    end
    subsection bc 2
      set type               = periodic
      set id                 = 4
      set periodic_id        = 5
      set periodic_direction = 2
    end
  end

First, the ``number`` of boundary conditions to be applied must be specified. For each boundary condition, the ``id`` of the boundary as well as its ``type`` must be specified. All boundaries are periodic. The ``x-`` (id=0) is periodic with the ``x+``boundary (id=1), the ``y-`` (id=2) is periodic with the ``y+``boundary (id=3) and so on and so forth. For each periodic boundary condition, the periodic direction must be specified. A periodic direction of ``0``implies that the normal direction of the wall is the :math:`\mathbf{e}_x` vector, ``1``implies that it's the :math:`\mathbf{e}_y`.

Physical Properties
~~~~~~~~~~~~~~~~~~~

The Reynolds number of 1600 is set solely using the kinematic viscosity since the reference velocity is one.
.. code-block:: text

  subsection physical properties
    set number of fluids = 1
    subsection fluid 0
      set kinematic viscosity = 0.000625
    end
  end


FEM Interpolation
~~~~~~~~~~~~~~~~~

The results obtained for the Taylor-Green vortex are highly dependent on the numerical dissipation that occurs within the CFD scheme. Generally, high-order methods outperform traditional second-order accurate methods for this type of flow. In the present case, we will investigate the usage of both second and third degree polynomial.

.. code-block:: text

    subsection FEM
        set velocity order = 2 #3 for Q3
        set pressure order = 2 #3 for Q3
    end

.. note::

Post-processing
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

  subsection post-processing
    set verbosity                = verbose
    set calculate enstrophy      = true
    set calculate kinetic energy = true
  end

To monitor the kinetic energy and the enstrophy, we set both calculation to ``true`` in the post-processing section.

Simulation Control
~~~~~~~~~~~~~~~~~~~~

The ``simulation control`` subsection controls the flow of the simulation. To maximise the temporal accuracy of the simulation, we use a third order ``bdf3`` scheme. Results are written every 2 time-step. To ensure a more adequate visualization of the high-order elements, we set ``subdivision=3``. This will allow paraview to render the high-order solutions with more fidelity.

.. code-block:: text

  subsection simulation control
    set method            = bdf3
    set time step         = 0.05 
    set number mesh adapt = 0    
    set time end          = 20  
    set output frequency  = 2    
    set subdivision       = 3
  end


Non-Linear Solver - Matrix-based
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Linear Solver - Matrix-based
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Linear Solver - Matrix-free
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``lethe-fluid-matrix-free`` has significantly more parameters for its linear solver. The new parameters are all related to the geometric multigrid preconditioner that is used by the matrix free algorithm.

.. code-block:: text

  subsection linear solver
    subsection fluid dynamics
      set method            = gmres
      set max iters         = 100
      set relative residual = 1e-4
      set minimum residual  = 1e-7
      set preconditioner    = gcmg
      set verbosity         = verbos
      
      #MG parameters
      set mg verbosity       = quiet
      set mg min level       = -1
      set mg level min cells = 16

      #smoother
      set mg smoother iterations = 10
      set mg smoother eig estimation = true
      
      # Eigenvalue estimation parameters
      set eig estimation degree          = 3
      set eig estimation smoothing range = 5
      set eig estimation cg n iterations = 20
      set eig estimation verbosity       = verbose

      #coarse-grid solver
      set mg coarse grid max iterations     = 2000
      set mg coarse grid tolerance          = 1e-7
      set mg coarse grid reduce             = 1e-4
      set mg coarse grid max krylov vectors = 30
      set mg coarse grid preconditioner     = ilu
      set ilu preconditioner fill               = 1
      set ilu preconditioner absolute tolerance = 1e-10
      set ilu preconditioner relative tolerance = 1.00
    end
  end

We set ``mg verbosity = quiet`` to prevent logging of the multigrid parameters during the simulation. Setting ``mg min level = -1`` ensures that the ``mg level min cells=16``parameter is used to determine the coarsest level. It is important to ensure that the Taylor-Green vortex has sufficient cell on the coarsest level since periodic boundary conditions are used. Indeed, using a coarsest level with a single cell can lead to a problematic situation where too few degrees of freedom are available on the coarsest level.

The ``smoother``, ``Eigenvalue estimation parameters`` and ``coarse-grid solver`` subsections are explained in the **Theory Guide** (under construction).

Rest of the Subsections
~~~~~~~~~~~~~~~~~~~~~~~~

The ``non-linear solver`` and ``linear solver`` subsections do not contain any new information in this example.


----------------------
Running the Simulation
----------------------
Launching the simulation is as simple as specifying the executable name and the parameter file. Assuming that the ``lethe-fluid`` executable is within your path, the simulation can be launched by typing:

.. code-block:: text
  :class: copy-button

  lethe-fluid taylor-couette.prm

Lethe will generate a number of files. The most important one bears the extension ``.pvd``. It can be read by visualization programs such as `Paraview <https://www.paraview.org/>`_.


----------------------
Results and Discussion
----------------------



----------------------------
Possibilities for Extension
----------------------------

- Calculate the order of convergence for the torque :math:`T_z`.
- It could be very interesting to investigate this flow in 3D at a higher Reynolds number to see the apparition of the Taylor-Couette instability. This, however, would be a major undertaking. 


------------
References
------------

[1] https://en.wikipedia.org/wiki/Main_Page