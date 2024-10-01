====
Mesh
====

The mesh subsection of DEM simulations is almost identical to the `CFD <https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/mesh.html>`_ in Lethe. There are two additional parameters mainly used in DEM and CFD-DEM simulations. These parameters are ``check diamond cells`` and ``expand particle-wall contact search``.

.. code-block:: text

 subsection mesh
    set check diamond cells                 = false

    set expand particle-wall contact search = false
 end

* ``check diamond cells`` parameter enables searching for diamond-shaped boundary cells and adding them to particle-wall contact search cells. The following image shows a diamond-shaped boundary cell colored in red.

.. image:: images/diamond_cell.png
    :alt: Schematic
    :align: center
    :width: 400


* ``expand particle-wall contact search`` enables adding the boundary neighbor cells of boundary cells to the particle-wall contact search list. This feature should only be activated in geometries with concave boundaries. (For example, for particles flow inside a cylinder or sphere). In convex geometries, enabling this feature leads to unphysical contact between particles and imaginary (unreal) extension of the boundary faces. The following image shows the boundary neighbor cells (colored in teal) of the read boundary cell.

.. image:: images/expand_particle_wall.png
    :alt: Schematic
    :align: center
    :width: 400


.. warning:: 
     In geometries with convex boundaries, this feature MUST NOT be activated.