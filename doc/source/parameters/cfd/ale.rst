=====================================
Arbitrary Lagrangian-Eulerian
=====================================

Lethe contains a small Arbitrary Lagrangian-Eulerian (ALE) module which enables the user to simulate domains which are moving at time and space-dependent velocities. All physics support the ALE module.


.. code-block:: text

   subsection ALE
    set enable                = false

    subsection uvw
      set Function expression = 0; 0 # In 2D : u;v
        or
      set Function expression = 0; 0; 0 #In 3D u;v;w
    end
   
   end

* The ``enable`` parameter is set to true if the ALE velocity should be substracted from the fluid velocity.

* The ``Function expression`` parameter sets the expression for the ALE velocity field in regards to :math:`u` and :math:`v`  for a 2D simulation and to :math:`u`, :math:`v`, :math:`w` for a 3D simulation.

.. warning:: The current implementation of the ALE module does not support time-dependent or non-divergence velocity fields


