Interface sharpening
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When running a simulation containing multiple fluids using the volume of fluid method (VOF), the interface between two fluids diffuses. 
To avoid that, this subsection includes parameters related to the interface sharpening. The default parameters are:

.. code-block:: text

  subsection interface sharpening
   set sharpening threshold   = 0.5
   set interface sharpness    = 2
   set sharpening frequency   = 10
  end

* The ``sharpening threshold`` is the phase fraction threshold, and therefore takes a value between 0 and 1. It defines the interface sharpening threshold that represents the mass conservation level.

* The ``interface sharpness`` parameter defines the sharpness of the moving interface (parameter alpha in the `interface sharpening model <https://royalsocietypublishing.org/doi/abs/10.1098/rsta.1952.0006>`_).

.. tip::

  This parameter must be larger than 1 for interface sharpening. Choosing values less than 1 leads to interface smoothing instead of sharpening.

* The ``sharpening frequency`` parameter defines the VOF interface sharpening frequency. It sets at what frequency the interface sharpening calculations should be printed out.

.. note::

  The XX examples well represents que sharpening issue.