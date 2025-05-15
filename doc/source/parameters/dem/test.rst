====
Test
====

This subsection deals with simulation testing. If it is set to ``true``, the sorted (using particle id) information of all particles will be printed in `xyz` format at the end of the simulation. Then this information can be compared with predefined values to ensure the reproducibility of the simulations.

.. code-block:: text

  subsection test
    # DEM test
    # Choices are true|false
    set type = false
  end

