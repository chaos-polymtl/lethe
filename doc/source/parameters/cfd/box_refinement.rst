..
  SPDX-FileCopyrightText: Copyright (c) 2022-2023, 2026 The Lethe Authors
  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

==============
Box Refinement
==============

The box refinement section allows for a specific region in the grid to be finer before the simulation starts. To do so, a box refinement can be added with the following parameters and their appropriate values:

.. code-block:: text

  subsection box refinement
  set number of refinement boxes = 0
    subsection box 0
      subsection mesh
        set type                   = dealii
        set grid type              = hyper_cube
        set grid arguments         = -1 : 1 : true
        set initial refinement     = 0
        set scale                  = 1
        set initial translation    = 0, 0, 0
        set initial rotation axis  = 1, 0, 0
        set initial rotation angle = 0
      end
      set additional refinement = 0
    end
  end

* The ``number of refinement boxes`` corresponds to the number of refinement boxes we want to introduce.

  .. attention::
    At the moment, the software only supports up to :math:`5` different boxes.

* The ``subsection box 0`` includes information regarding the 1st refinement box.

  .. note::
    When ``number of refinement boxes = 2``, the 2nd refinement box will be described in the ``subsection box 1``. We increment the last digit of the subsection name to describe a new refinement box.

* The ``mesh`` subsection allows to define the region in which the fluid mesh needs to be uniformly refined. A cell will be refined if at least one of its degrees of freedom (DoFs) is located within the outer boundaries of the box specified in the ``grid arguments``. Therefore, in the 3D example above, every cell of the fluid mesh that has at least one of its DoFs located in the hexahedron located between (-1, -1, -1) and (1,1,1) will be refined.

  .. note::
    The used mesh can be of any ``type`` and any ``grid type``. Additionally, box refinement meshes support scaling, translation, and rotation, applied in the listed order.

    .. seealso::
      For more information on meshes, see :doc:`the mesh documentation <./mesh>`.

  .. tip::
    The ``initial refinement`` of the ``hyper_cube`` should be as small as possible, since the initial refinement of the box mesh itself will not have any impact on the definition of the refinement zone.
    However, in the case of a curved shape such as a ``hyper_sphere``, increasing the ``initial refinement`` of the box refinement mesh improves the geometry resolution of the box boundary.

* The ``additional refinement`` parameter in the ``box refinement`` subsection indicates the number of times the cells inside the box will be refined in addition to the ones prescribed in the :doc:`simulation mesh subsection<./mesh>`.

  .. warning::
    This number is not the level of refinement of the cell, but the number of additional refinements applied within the region of the box at the initial condition. For instance, let us consider a simulation mesh with an ``initial refinement`` of :math:`1` and a box refinement mesh with an ``additional refinement`` of :math:`2`. The level of refinement of cells within the box refinement region becomes: :math:`1 + 2 = 3`. However, the mesh resolution within this region, may change dynamically if :doc:`adaptive mesh refinement<./mesh_adaptation_control>` is enabled.
