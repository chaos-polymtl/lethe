Interface sharpening
--------------------

When running a multiphase flow simulation using the volume of fluid method (VOF), the interface between two fluids becomes more and more blurry after each time step due to numerical diffusion. 

To avoid that, this subsection includes parameters related to the interface sharpening.

.. code-block:: text

  subsection interface sharpening
   set sharpening threshold   = 0.5
   set interface sharpness    = 2
   set sharpening frequency   = 10
  end

.. warning::
   Do not forget to ``set interface sharpening = true`` in the :doc:`multiphysics` subsection of the ``.prm``.   
  
   
* ``sharpening threshold``: phase fraction threshold, between 0 and 1. It defines the interface sharpening threshold that represents the mass conservation level.

* ``interface sharpness``: sharpness of the moving interface (parameter :math:`a` in the `interface sharpening model <https://www.researchgate.net/publication/287118331_Development_of_efficient_interface_sharpening_procedure_for_viscous_incompressible_flows>`_).

.. tip::

  This parameter must be larger than 1 for interface sharpening. Choosing values less than 1 leads to interface smoothing instead of sharpening. A good value would be between 1 and 2.

* ``sharpening frequency``: sets the frequency (in number of iterations) for the interface sharpening computation.

.. seealso::

  The :doc:`../../examples/multiphysics/dam-break-VOF/dam-break-VOF` example using VOF represents well the interface sharpening issue.
