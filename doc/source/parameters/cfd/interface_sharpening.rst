Interface sharpening
--------------------

When running a simulation containing multiple fluids using the volume of fluid method (VOF), the interface between two fluids becomes more blurry after each time step due to diffusion. 
To avoid that, this subsection includes parameters related to the interface sharpening. The default parameters are:

.. code-block:: text

  subsection interface sharpening
   set sharpening threshold   = 0.5
   set interface sharpness    = 2
   set sharpening frequency   = 10
  end

.. warning::
   Do not forget to ``set interface sharpening = true`` in the :doc:multiphysics subsection of the ``.prm``.   
   
.. warning::
   Do not forget to ``set interface sharpening = true`` in the :doc:multiphysics subsection of the ``.prm``.   
   
* The ``sharpening threshold`` is the phase fraction threshold, and therefore takes a value between 0 and 1. It defines the interface sharpening threshold that represents the mass conservation level.

* The ``interface sharpness`` parameter defines the sharpness of the moving interface (parameter :math:`a` in the `interface sharpening model <https://www.researchgate.net/publication/287118331_Development_of_efficient_interface_sharpening_procedure_for_viscous_incompressible_flows>`_).

.. warning::

  This parameter must be larger than 1 for interface sharpening. Choosing values less than 1 leads to interface smoothing instead of sharpening. A good value would be between 1 and 2.

* The ``sharpening frequency`` parameter defines the VOF interface sharpening frequency. It sets at what frequency the interface sharpening calculations should be carried out.

.. seealso::

  The :doc:`../../examples/multiphysics/dam-break-VOF/dam-break-VOF` example using VOF represents well the interface sharpening issue.
