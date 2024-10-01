==================================
Rotating Drum with Postprocessing
==================================

This is a mini-example that only adds Lagrangian post-processing features to the rotating drum example (example 3). Hence, we only explain the post-processing subsection in this example.

----------------------------------
Features
----------------------------------
- Solvers: ``lethe-particles``
- Rotational boundary
- Load-balancing
- Lagrangian post-processing


----------------------------
Files Used in This Example
----------------------------

- Parameter file: ``examples/dem/3d-rotating-drum-with-postprocessing/rotating-drum-with-postprocessing.prm``


-----------------------
Description of the Case
-----------------------

This example is identical the rotating drum example. The only difference is that in this example, we use Lagrangian post-processing to obtain granular temperature and average velocity (averaged in cells) distribution.


--------------
Parameter File
--------------

Post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we set the variable ``Lagrangian post-processing`` equal to ``true``. This enables Lagrangian post-processing calculations. Lethe-DEM built-in post-processing capabilities are detailed in the :doc:`../../../parameters/dem/post-processing` section of the parameters guide.

.. code-block:: text

    subsection post-processing
      set Lagrangian post-processing = true
    end

