===============
Mesh Adaptation
===============

This subsection controls the mesh adaptation method, with default values given below.

.. code-block:: text

  subsection mesh adaptation
    # Type of mesh adaptation. Choices are  none, uniform or kelly.
    set type                     = none

    # Variable(s) for kelly estimation. Choices are velocity, pressure, phase or temperature.
    # For multi-variables refinement, separate the different variables with a comma
    set variable                 = velocity

    # Frequency of the mesh refinement
    set frequency                = 1

    # Minimum refinement level
    set min refinement level     = 0

    # Maximum refinement level
    set max refinement level     = 10

    # Fraction of coarsened elements
    # For multi-variables refinement, separate the different fractions with a comma
    set fraction coarsening      = 0.05

    # Fraction of refined elements
    # For multi-variables refinement, separate the different fractions with a comma
    set fraction refinement      = 0.1

    # How the fraction of refinement/coarsening are interpreted
    # Choices are number or fraction
    set fraction type            = number

    # Maximum number of elements
    set max number elements      = 100000000

    # Enabling the mesh adaptation controller to obtain a constant number of elements
    mesh refinement controller   = false

    # Number of initial (pre-solve) refinement steps
    set initial refinement steps = 0

    # Fix the boundary mesh size at the initial refinement level when the kelly mesh refinement estimator is used.
    set fix boundary refinement = false

    # List of boundaries next to which the mesh refinement level should not be modified. The list must contain the boundary ids separated by commas.
    set boundaries fixed = 0, 1
  end


* Two ``type`` of mesh adaptation are available. The ``uniform`` mesh adaptation refines the mesh at every cell, whereas the ``kelly`` uses a `kelly error estimator <https://www.dealii.org/current/doxygen/deal.II/classKellyErrorEstimator.html>`_ to decide which cell are refined, by estimating the error per cell for a given variable. 
* The variable for kelly estimation should be specified with ``set variable``, and can be: velocity, pressure, phase (for multiphase flows), temperature
	* Mesh adaptation can be defined on multiple variables, separated with a coma (e.g. ``set variable = velocity,temperature``, or ``set variable = velocity,phase,pressure`` etc.).

	.. warning::
		The different ``fraction refinement`` and ``fraction coarsening`` must be defined explicitly (see these parameters definition below).

.. important::
	In the case of multiple variable mesh adaptation, the cells are:
		* refined if refinement is necessary for at least one variable
		* coarsened if coarsening is necessary for *all* variables

* The frequency at which the mesh is refined is controlled with the ``frequency`` parameter. If ``set frequency = 1``, the mesh is refined at every iteration. 
	* For transient simulation, this means at every time-step. 
	* For steady-state simulation in which the steady-state problem is solved on successively refined meshes, the user should have ``set frequency = 1``, which is the default value.

* The minimal and maximal refinement level reachable for a cell are controlled respectively with the ``min refinement`` and ``max refinement`` parameters.
   * for ``deal.ii`` meshes, if the ``min refinement level`` is equal to the ``initial refinement`` (see `Mesh paramater <https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/mesh.html>`_), no cell will be coarser than the initial mesh.
   * for ``gmsh`` imported meshes, if ``set min refinement level = 0``, no cell will be coarser than the initial mesh.

.. tip:: 
	For a ``gmsh`` mesh, a cell cannot be coarsened more than it's initial level. Consequently, adaptively refined simulations should start with a mesh as coarse as possible. 

.. tip:: 
	For a good compromise between speed and precision, ``max refinement level`` should be set to ``2`` or ``3`` more than the ``min refinement level``

* The fraction of cell that are refined and coarsened are controlled with the ``fraction refinement`` and ``fraction coarsening`` parameters. 
	* Fractions for mesh adaptation on multiple variables must be separated with a coma (e.g. ``set fraction refinement = 0.2,0.2``, or ``set fraction coarsening = 0.1,0.3`` etc.).

	.. warning::
		The different ``variable`` must be defined explicitly (see this parameter definition above).

.. tip:: 
	For ``set type = kelly``, and ``set variable = velocity`` or ``pressure``, a good first start is achieve with ``set fraction refinement = 0.2`` and ``set fraction coarsening = 0.1``.

	For ``set type = kelly``, and ``set variable = phase``, use ``fraction type = fraction`` (explained below) and ``set fraction refinement = 0.8`` for a good tracking of the entire free surface (see `Multiphysics <file:///home/jeannej/Softwares/lethe/lethe/doc/build/html/parameters/cfd/multiphysics.html>`_).

* The fraction of refinement/coarsening can be interpreted in ``number`` or ``fraction``  depending on the parameter ``fraction type``. At first sight, this is a relatively difficult concept to understand that is inherited from deal.II. 
	* When ``fraction type = number``  the  `refine_and_coarsen_fixed_number <https://www.dealii.org/current/doxygen/deal.II/namespaceGridRefinement.html#a48e5395381ed87155942a61a1edd134d>`_ strategy of deal.II is used. This function provides a strategy to mark cells for refinement and coarsening with the goal of providing predictable growth in the size of the mesh by refining  and coarsening a given fraction of all cells.  
	* When ``fraction type = fraction``,  the `refine_and_coarsen_fixed_fraction <https://www.dealii.org/current/doxygen/deal.II/namespaceGridRefinement.html#ae90dc87c4db158b8d01f6d564ac614e5>`_ strategy is used. This function provides a strategy to mark cells for refinement and coarsening with the goal of controlling the reduction of the error estimate. Also known as the bulk criterion or Dörfler marking, this function computes the thresholds for refinement and coarsening such that the criteria of cells getting flagged for refinement make up for a certain fraction of the total error.


* The maximum number of elements in the entire domain can be controlled with the ``max number elements`` parameter.

* The boolean parameter ``mesh refinement controller`` activates a controller that overrides the value of the ``fraction coarsening`` parameter. If activated, the controller will try to maintain the total number of elements in the domain equal to the value of ``max number elements`` parameter. The control is done using a PID controller.

.. note:: 
    If the ``fraction refinement`` parameter is too high, the controller may not be able to maintain the number of elements constant. If ``fraction type = number``, the maximal ``fraction refinement`` that is stable in 3D is 0.125. In 2D, it is 0.25.

.. tip:: 
	When using the mesh refinement controller, try reducing the ``fraction refinement`` parameter if elements alternate between being refined and coarsened.

.. warning::
	If ``mesh refinement controller`` is set to ``false``, the ``max number elements`` parameter puts a hard limit on the number of cells in the domain, even if the ``fraction refinement`` is increased. 

* The number of initial (before solving) adaptive refinement steps is controlled by the ``initial refinement steps`` parameter. With an ``initial refinement steps`` larger than 0, the triangulation is refined adaptively before the solver starts solving the problem. This enables the user to adapt the initial mesh to the initial condition. For example, if the simulation is a VOF simulation, it is ideal to have an initial mesh that captures the interface between the fluids accurately. This is achieved by refining the mesh using the dynamic mesh adaptation parameters and reapplying the initial condition after each adaptation. This process will be repeated ``initial refinement steps`` times.

* The ``fix boundary refinement`` parameter is used to fix the boundary mesh size at the defined initial refinement level. ``boundaries fixed`` can be used to specify which boundary ids will not be refined or coarsened when the kelly mesh adaptation is used. Cells on the boundary will still be refined when global refinement is active or when neighbouring cells force the refinement of a boundary mesh cell.