Mesh Adaptation Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This subsection controls the mesh adaptation method, with default values given below.

.. code-block:: text

	subsection mesh adaptation
	  # Type of mesh adaptation. Choices are  none, uniform or kelly.
	  set type                 = none

	  # Variable for kelly estimation. Choices are velocity, pressure or phase.
	  set variable             = velocity

	  # Frequency of the mesh refinement
	  set frequency            = 1

	  # Minimum refinement level
	  set min refinement level = 0

	  # Maximum refinement level
	  set max refinement level = 10

	  # Fraction of coarsened elements
	  set fraction coarsening  = 0.05

	  # Fraction of refined elements
	  set fraction refinement  = 0.1

	  # How the fraction of refinement/coarsening are interpreted
	  # Choices are number or fraction 
	  set fraction type        = number

	  # Maximum number of elements
	  set max number elements  = 100000000
	end

* Two ``type`` of mesh adaptation are available. The ``uniform`` mesh adaptation refines the mesh at every cell, whereas the ``kelly`` uses a `kelly error estimator <https://www.dealii.org/current/doxygen/deal.II/classKellyErrorEstimator.html>`_ to decides which cell are refined, by estimating the error per cell for a given variable. The variable for kelly estimation should be specified with ``set variable``, and can be:
	* velocity
	* pressure
	* phase (for multiphase flows)
* The frequency at which the mesh is refined is controlled with the ``frequency`` parameter. If ``set frequency = 1``, the mesh is refined at every time step.

* The minimal and maximal refinement level reachable for a cell are controlled respectively with the ``min refinement`` and ``max refinement`` parameters.
   * for ``deal.ii`` meshes, if the ``min refinement level`` is equal to the ``initial refinement`` (see `Mesh paramater <https://lethe-cfd.github.io/lethe/parameters/cfd/mesh.html>`_), no cell will be coarser than the initial mesh.
   * for ``gmsh`` imported meshes, if ``set min refinement level = 0``, no cell will be coarser than the initial mesh.
.. tip:: 
	For a ``gmsh`` mesh, a cell cannot be coarsened more than it's initial level. Consequently, adaptively refined simulations should start with a mesh as coarse as possible. 
.. tip:: 
	For a good compromise between speed and precision, ``max refinement level`` should be set to ``2`` or ``3`` more than the ``min refinement level``
* The fraction of cell that are refined and coarsened are controlled with the ``fraction refinement`` and ``fraction coarsening`` parameters. 
.. tip:: 
	For ``set type = kelly``, and ``set variable = velocity`` or ``pressure``, a good first start is achieve with ``set fraction refinement = 0.2`` and ``set fraction coarsening = 0.1``, meaning the top 20% of the cells with the highest kelly error estimator will be refined and the top 10% of the cells with the lowest estimator will be coarsened.

	For ``set type = kelly``, and ``set variable = phase``, use ``set fraction refinement = 0.8`` for a good tracking of the entire free surface (see `Multiphysics <file:///home/jeannej/Softwares/lethe/lethe/doc/build/html/parameters/cfd/multiphysics.html>`_).
* The fraction of refinement/coarsening can be interpreted in ``number`` or ``fraction`` in the parameter ``fraction type`` `((explanations needed here))`
* The maximum number of elements in the entire domain can be controlled with the ``max number elements`` parameter.
.. warning::
	The ``max number elements`` parameter puts a hard limit on the number of cells in the domain, even if the ``fraction refinement`` is increased.

