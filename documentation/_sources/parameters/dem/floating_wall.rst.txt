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
      subsection point on wall
        set x = 0
        set y = 0
        set z = -0.06
      end

      # Defining normal vector of the floating wall
      subsection normal vector
        set nx = 0
        set ny = 0.5
        set nz = 1
      end

      # Starting time of the floating wall
      set start time = 0

      # Ending time of the floating wall
      set end time   = 0.2
    end
  end


* ``number of floating walls`` parameter defines the total number of floating walls we wish to insert during the simulation.

* For each floating wall, we need a separate subsection (for instance,	``subsection wall 0``) in which the information of the floating wall (the normal vector, start and end times, and a point on the floating wall) is defined.

* In the subsection ``point on wall``, we define a point (with coordinates ``x``, ``y``, and ``z``) on the floating wall.

* In the subsection ``normal vector of the floating wall``, we define the normal vector (with components of ``nx``, ``ny``, and ``nz`` in `x`, `y`, and `z` directions) of the floating wall.

* The ``start time`` parameter specifies the start time at which the floating wall is applied in the simulation.

* The ``end time`` parameter defines the time at which the floating wall is removed from the simulation. In simulations that the user wishes to keep a floating wall during the entire simulation, ``end time`` should be specified greater than the total simulation time.

