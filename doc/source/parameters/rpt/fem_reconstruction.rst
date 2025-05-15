===================
FEM Reconstruction
===================

In this subsection, the parameters used for ``lethe-rpt-l2-projection-3d`` and ``lethe-rpt-fem-reconstruction-3d`` applications are defined.
Here are the default values:

.. code-block:: text

  subsection fem reconstruction
    set mesh type                           = dealii
    set mesh filename                       = none
    set z subdivisions                      = 1
    set mesh refinement                     = 0
    set l2 projection before reconstruction = false
    set experimental counts file            = none
    set export positions file               = none
    set search extrapolation limit          = 0.005
    set cost function type                  = relative
    set dof handler file                    = none
    set nodal counts file                   = none
    set search type                         = local
    set search cell proximity level         = 1
    set verbose clock                       = false
  end

For both ``lethe-rpt-l2-projection-3d`` and ``lethe-rpt-fem-reconstruction-3d`` applications, we have to define the mesh:

- ``mesh type``: Type of mesh used. Choosing the ``dealii`` option will generate a ``subdivided_cylinder`` grid, for that reason the number of subdivision in the z direction must be specified with the ``z subdivisions`` parameter. For the ``gmsh`` option, only meshes for a cylindrical geometry with tetrahedral elements are accepted at the moment. Furthermore, the z-axis should be the axis of symmetry of the cylinder.
    Options: ``dealii`` or ``gmsh``
- ``mesh filename``: Filename of the imported mesh. Specify only if ``gmsh`` was the selected option for the ``mesh type``.
    Options: Any ``.msh`` file
- ``z subdivisions``: Number of initial subdivisions along the z-axis. Specify only if ``dealii`` was the selected option for the ``mesh type``.
    Options: Any positive integer
- ``mesh refinement``: Number of global mesh refinements. Specify only if ``dealii`` was the selected option for the ``mesh type``.
    Options: Any positive integer

For the ``lethe-rpt-fem-reconstruction-3d`` application only, we have to define the following parameters:

- ``l2 projection before reconstruction``: Enable to run the ``lethe-rpt-l2-projection-3d`` application before the reconstruction.
    Options : ``true`` or ``false``
- ``experimental counts file``: Filename of the file containing the experimental photon counts.
    Options: Any ``.txt`` file
- ``export positions file``: Name of the file in which the calculated positions are exported.
    Options: Any ``.csv`` or ``.dat`` filename
- ``cost function type``: Type of cost function used to find the position.
    Options:``relative`` or ``absolute``
- ``dof handler file``: File containing the saved DoF handler.
    Options: Any ``.dof`` file
- ``nodal counts file``: List of files containing the nodal counts from each detector.
    Options: Any ``.counts`` file
- ``search extrapolation limit``: Tolerance when extrapolating from a cell in the reference space to find a particle's position. The default value is set to :math:`0.005`. However, for a subdivided cylinder geometry from deal.II, the default tolerance is set to :math:`1.15̀` times the ``cell_height``.
    Options: Any positive double
- ``search type``: Type of search algorithm used to find particle positions. The ``local`` option refers to a search algorithm where the next particle position is searched in a scope around the previously found position's cell. The size of that scope is defined by the ``search cell proximity level`` parameter. And the ``global`` option refers to the algorithm  in which we go through every cell of the grid to find every particle position.
    Options: ``local`` or ``global``
- ``search cell proximity level``: Level of proximity of the search scope to find the next particle position. A ``search cell proximity level = 1`` includes in the search scope the previously found position's cell and all its adjacent cells. Specify only if ``local`` was the selected option for the ``search type``.
    Options: Any positive integer
- ``verbose clock``: Enable to show total wallclock time elapsed since start.
    Options: ``true`` or ``false``