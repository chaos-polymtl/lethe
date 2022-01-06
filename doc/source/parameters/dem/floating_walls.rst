Floating Walls
-------------------
.. note::
    Floating wall is a temporary (its start and end times are defined) flat wall, generally used for holding the particles during the filling and before the discharge stage.

In this subsection, the information on floating walls is defined. First of all, the total ``number of floating walls`` is specified. Then for each floating wall, we should specify its ``normal vector``, a ``point on the wall``, ``start`` and ``end times``.

.. code-block:: text

 subsection floating walls
  # Total number of floating walls
    set number of floating walls         = 1
    	subsection wall 0

  # Defining a point on the floating wall
		subsection point on wall
			set x                        = 0
			set y                        = 0
			set z                        = -0.06
		end

  # Defining normal vector of the floating wall
		subsection normal vector
			set nx                       = 0	
			set ny                       = 0.5
			set nz                       = 1
		end

  # Starting time of the floating wall
		set start time                    = 0

  # Ending time of the floating wall
		set end time                      = 0.2
    	end
 end


