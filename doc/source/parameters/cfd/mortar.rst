======
Mortar
======

The mortar section is used when simulating rotor-stator geometries, in which the rotor mesh is attached to the stator by mortar elements. 

.. code-block:: text

  subsection mortar
    set enable = true
    subsection mesh
      set type                   = dealii
      set grid type              = subdivided_hyper_rectangle
      set grid arguments         = 1,1,1 : -1,-1,-1 : 1,1,1 : true
      set initial refinement     = 0
    end
    set rotor boundary id   = 4
    set stator boundary id  = 2
    set center of rotation  = 0, 0
    set rotation axis       = 0, 0, 1
    subsection rotor rotation angle
      set Function expression = 0
    end
    subsection rotor angular velocity
      set Function expression = 0
    end
    set penalty factor      = 1.0
    set oversampling factor = 2
    set verbosity           = verbose
  end

* The mesh parameters in the :doc:`../cfd/mesh` subsection refer to the stator domain. The ``mesh`` subsection herein mentioned contains the parameters of the rotor domain; nonetheless, the input format is the same as in :doc:`../cfd/mesh`.

.. note::
  The initial number of cells at the rotor-stator interface has to be the same; the simulation will be aborted if that is not respected. This restriction will be automatically constrained throughout the simulation if the mesh is refined.
 
* The ``rotor boundary id`` and ``stator boundary id`` refer to the boundary index at the rotor-stator interface.

.. warning::
  In ``dealii`` meshes the boundary IDs are automatically assigned to the geometries using the deal.II `colorization <https://www.dealii.org/current/doxygen/deal.II/DEALGlossary.html#GlossColorization>`_ function, and thus the rotor/stator IDs might be duplicated.
  To circumvent this, the rotor boundary IDs are shifted. The ``rotor boundary id`` entry refers to the shifted ID number, assuming that the enumeration starts sequentially from the last entry of the stator boundary IDs.

* The ``center of rotation`` is the reference point for the prescribed rotation at the rotor domain.

* The ``rotation axis`` is the unit vector that defines the rotor axis of rotation.

* The ``rotor rotation angle`` subsection allows the imposition of a constant or time-dependent rotation angle for the rotor.

* The ``rotor angular velocity`` subsection allows the imposition of a constant or time-dependent angular velocity for the rotor.

.. warning::
  The ``rotor angular velocity`` expression needs to correspond to the time derivative of the ``rotor rotation angle``.

* The ``penalty factor`` is used for the weak imposition of the mortar coupling at the interface. This parameter is akin to the symmetric interior penalty factor in SIPG (Symmetric Interior penalty Galerkin Method) [#larson2013]_.

* The ``oversampling factor`` is used to increase the number of quadrature points. This feature is used to better approximate the weak imposition of the interface coupling.

* When enabling ``verbosity`` (``set verbosity = verbose``), the rotor rotation information is printed at every iteration.

Reference
---------
.. [#larson2013] \M. G. Larson, F. Bengzon. *The Finite Element Method: Theory, Implementation, and Applications*, vol. 10. Springer Berlin, Heidelberg, 2013. 
