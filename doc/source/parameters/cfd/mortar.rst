======
Mortar
======

The mortar section is used when simulating rotor-stator geometries, in which the rotor (fluid) mesh is attached to the stator by mortar elements. 

.. code-block:: text

  subsection mortar
    set enable = true
    subsection mesh
      set type                   = dealii
      set grid type              = subdivided_hyper_rectangle
      set grid arguments         = 1,1,1 : -1,-1,-1 : 1,1,1 : true
      set initial refinement     = 0
      set initial rotation angle = 0
    end
    set rotor boundary id  = 4
    set stator boundary id = 2
    set center of rotation = 0, 0
  end

* The mesh parameters in the :doc:`../cfd/mesh` subsection refer to the stator domain. The ``mesh`` subsection herein mentioned contains the parameters of the rotor domain; nonetheless, the input format is the same as in :doc:`../cfd/mesh`.

.. note::
  The initial number of cells at the rotor-stator interface has to be the same; the simulation will be aborted if that is not respected. This restriction will be automatically constrained throughout the simulation if the mesh is refined.
 
* The ``rotor boundary id`` and ``stator boundary id`` refer to the boundary index at the rotor-stator interface.

.. warning::
  In ``dealii`` meshes the boundary IDs are automatically assigned to the geometries using the deal.II `colorization <https://www.dealii.org/current/doxygen/deal.II/DEALGlossary.html#GlossColorization>`_ function, and thus the rotor/stator IDs might be duplicated.
  To circumvent this, the rotor boundary IDs are shifted. The ``rotor boundary id`` entry refers to the shifted ID number, assuming that the enumeration starts sequentially from the last entry of the stator boundary IDs.

* The ``center of rotation`` is the reference point for the prescribed rotation at the rotor domain.