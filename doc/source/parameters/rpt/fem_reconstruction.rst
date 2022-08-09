FEM Reconstruction
-------------------

In this subsection are defined the parameters used for the ``rpt_l2_projection_3d`` and ``rpt_fem_reconstruction_3d`` applications.
Here are the default values:

.. code-block:: text

    #---------------------------------------------------
    # FEM reconstruction parameters
    #---------------------------------------------------
    subsection fem reconstruction
        set mesh type                  = dealii
        set mesh filename              = reactor.msh
        set z subdivisions             = 2
        set mesh refinement            = 2
        set experimental counts file   = experimental_counts.txt
        set export positions file      = found_positions.csv
        set cost function type         = relative
        set dof handler file           = temp_dof_handler.dof
        set nodal counts file          = temp_nodal_counts_detector00.counts
        set verbose clock              = false

- ``mesh type``: Type of mesh used.
    Options: ``dealii`` or ``gmsh``
- ``mesh filename``: Filename of the imported mesh. Specify only if ``gmsh`` was the selected option for the ``mesh type``.
    Options: Any ``.msh`` file
- ``z subdivisions``: Number of initial subdivisions along the z-axis. Specify only if ``dealii`` was the selected option for the ``mesh type``.
    Options: Any positive integer
- ``mesh refinement``: Number of global mesh refinements. Specify only if ``dealii`` was the selected option for the ``mesh type``.
    Options: Any positive integer
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
- ``verbose clock``: Enable to show total wallclock time elapsed since start
    Options: ``true`` or ``false``