=============
Floating Wall
=============

A floating wall is a temporary flat wall (its start and end times are defined), generally used for holding the particles during the filling and before the discharge stage.

In this subsection, the information on floating walls is defined. First of all, the total ``number of floating walls`` is specified. Then for each floating wall, we should specify its ``normal vector``, a ``point on the wall``, ``start`` and ``end times``.

.. code-block:: text

  subsection floating walls
    # Total number of floating walls
    set number of floating walls = 1

    subsection wall 0
      # Defining a point on the floating wall
      set point on wall = 0., 0., 0.

      # Defining normal vector of the floating wall
      set normal vector = 1., 0., 0.

      # Starting time of the floating wall
      set start time = 0

      # Ending time of the floating wall
      set end time   = 0.2
    end
  end


* ``number of floating walls`` parameter defines the total number of floating walls we wish to insert during the simulation.

* For each floating wall, we need a separate subsection (for instance,	``subsection wall 0``) in which the information of the floating wall (the normal vector, start and end times, and a point on the floating wall) is defined.

* The ``point on wall`` parameter defines a point on the floating wall in `x`, `y`, and `z` directions.

* The ``normal vector`` parameter defines the normal vector of the floating wall in `x`, `y`, and `z` directions.

* The ``start time`` parameter specifies the start time at which the floating wall is applied in the simulation.

* The ``end time`` parameter defines the time at which the floating wall is removed from the simulation. In simulations that the user wishes to keep a floating wall during the entire simulation, ``end time`` should be specified greater than the total simulation time.

