==============
Box Refinement
==============

The box refinement section allows for a specific region in the grid to be finer before the simulation starts. To do so, a box refinement can be added with the following example parameters:

.. code-block:: text

  subsection box refinement
    subsection mesh
      set type               = dealii
      set grid type          = subdivided_hyper_rectangle
      set grid arguments     = 1,1,1 : -1,-1,-1 : 1,1,1 : true
      set initial refinement = 0
    end
    set initial refinement = 3
  end

* The ``mesh`` subsection allows to define the region in which the fluid mesh needs to be refined. A cell will be refined if at least one of its degrees of freedom (dofs) is located within the outer boundaries of the box specified in the ``grid arguments``. Therefore, in this example, every cell of the fluid mesh that has at least one of its dofs located in the hexahedron located between (-1, -1, -1) and (1,1,1) will be refined. For more information on meshes, see :doc:`../cfd/mesh`. 

.. tip::
  The ``initial refinement`` of the ``subdivided_hyper_rectangle`` should be as small as possible, since the initial refinement of the box mesh itself will not have any impact on the definition of the refinement zone. 

.. note::
  The used mesh can be of any ``type`` and any ``grid type``.

* The ``initial refinement`` parameter in the ``box refinement`` subsection will dictate the number of times the cells inside the box will be refined. 