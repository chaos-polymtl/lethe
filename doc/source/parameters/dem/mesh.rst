..
  SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

====
Mesh
====

DEM and CFD simulations read their triangulation through the same routine, so the ``mesh`` subsection of a DEM simulation is the `CFD one <https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/mesh.html>`_. Every mesh type (``gmsh``, ``dealii`` and the ``lethe`` built-in grids), the mesh transformations (``scale``, ``initial translation``, ``initial rotation axis`` and ``initial rotation angle``), the initial refinement and the refinement in the vicinity of boundaries are all available in ``lethe-particles``.

The :doc:`manifolds <../cfd/manifolds>` and :doc:`box refinement <../cfd/box_refinement>` subsections are also read by ``lethe-particles``, and behave exactly as they do in the CFD solvers.

.. warning::
  Simplex meshes are **not** supported by DEM simulations. Setting ``simplex = true`` in the ``mesh`` subsection raises an error, since the particle-wall contact search, the load balancing and the checkpointing all require a ``parallel::distributed::Triangulation``.

Two parameters of the ``mesh`` subsection are only used by DEM and CFD-DEM simulations.

.. code-block:: text

 subsection mesh
    set check diamond cells                 = false

    set expand particle-wall contact search = false
 end

* ``check diamond cells`` parameter enables searching for diamond-shaped boundary cells and adding them to particle-wall contact search cells. The following image shows a diamond-shaped boundary cell colored in red. A diamond cell is defined as a cell that has a vertex or an edge on the boundary, but no face lying on the boundary.

.. image:: images/diamond_cell.png
    :alt: Schematic
    :align: center
    :width: 400


* ``expand particle-wall contact search`` enables adding the neighboring cells of boundary cells to the particle-wall contact search list. This feature should only be activated in geometries with convex boundaries (for example, particles flowing inside a cylinder or sphere). The following image shows the boundary neighbor cells (colored in teal) of the red boundary cell. In concave geometries, enabling this feature leads to unphysical contacts between particles and the imaginary (unreal) extension of the boundary faces from neighboring cells.

.. image:: images/expand_particle_wall.png
    :alt: Schematic
    :align: center
    :width: 400


.. warning::
     In geometries with convex boundaries, this feature MUST NOT be activated.
