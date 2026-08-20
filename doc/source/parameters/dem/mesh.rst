====
Mesh
====

The mesh subsection of DEM simulations is almost identical to the `CFD <https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/mesh.html>`_ in Lethe. There are two additional parameters mainly used in DEM and CFD-DEM simulations. These parameters are ``check diamond cells`` and ``expand particle-wall contact search``.

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


* ``expand particle-wall contact search`` enables adding the boundary faces of cells neighboring the boundary cell containing the particle, to the particle-wall contact search list. Only neighboring boundary faces that form a convex boundary with the current cell are added; boundary faces forming a non-convex boundary are excluded. The following figure illustrates this concept. When a particle is located in a boundary cell (shown in red), the boundary faces of neighboring cells that form a convex boundary with the current cell (neighboring cells shown in green) are also considered during the particle-wall contact search. This is necessary because a particle whose centroid lies within the red cell can simultaneously contact both the boundary face of the current cell and the boundary face of a neighboring cell when the two faces form a convex boundary. In contrast, boundary faces of neighboring cells that form a concave boundary (neighboring cells shown in orange) are not considered during the particle-wall contact search.

.. image:: images/expand_particle_wall.png
    :alt: Schematic
    :align: center
    :width: 400