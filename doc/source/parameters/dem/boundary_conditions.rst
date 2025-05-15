===================
Boundary Conditions
===================

In this subsection, the boundary conditions of the DEM simulation are defined. First of all, the ``number of boundary conditions`` is specified. Then for each boundary condition, its information is defined. There are five boundary types: ``fixed_wall``, ``outlet``, ``rotational`` (around the center), ``translational``, and ``periodic``. For ``rotational`` motion, ``rotational speed`` and ``rotational vector`` are required, while for ``translational`` motion, the ``speed`` should be defined in each direction. For ``periodic`` boundaries, ``periodic id 0``, ``periodic id 1`` and ``periodic direction`` are required.

``fixed_wall`` is a static wall, and particles collide with these static walls upon reaching the wall. The only way to move these walls is to move the entire triangulation. If the ``outlet`` condition is chosen for a boundary, particles can leave the simulation domain via this outlet. Using ``rotational`` and ``translational`` boundary conditions, exerts imaginary rotational and translational velocities to that boundary. In other words, the boundary does not move, but the particles that have collisions with these walls receive a rotational or translational velocity from the wall. This feature is used in the rotating drum example.

.. code-block:: text

  subsection DEM boundary conditions
    # Total number of boundary motion
    set number of boundary conditions = 3

    # For each motion, we need a separate subsection
    subsection boundary condition 0
      # ID of boundary
      set boundary id         = 0

      # Boundary type
      # Choices are fixed_wall|outlet|rotational|translational
      set type                = rotational

      # Rotational speed magnitude
      set rotational speed    = 2.5

      # Rotational vector
      set rotational vector   = 1, 0, 0

      # Point on rotational vector
      set point on rotational vector = 0, 0, 0
    end

    # OR for translational motion
    subsection boundary condition 1
      # ID of moving boundary
      set boundary id = 1

      # Motion type
      set type        = translational

      # Speed in each direction
      set speed x     = 0.15
      set speed y     = 0
      set speed z     = 0
    end
  end

* The ``number of boundary conditions`` parameter defines the number of desired boundary conditions to be specified. Note that if a boundary condition ``type`` is not defined explicitly, Lethe defines it as a fixed static wall.

* For each boundary condition, we have to define a separate subsection. In the sample parameter list above, the ``number of boundary conditions`` is equal to 2. Hence, we need to define two subsections (``subsection boundary condition 0`` and ``subsection boundary condition 1``).

* The ``boundary id`` parameter specifies the boundary ID for which the boundary condition should be applied, periodic boundaries are an exception.

* The ``periodic id 0`` and ``periodic id 1`` parameters specify the periodic boundaries ID for which the periodic boundary condition should be applied.

.. note::
        Only periodic boundaries which have co-linear normal vectors which align along one axis of the problem (e.g., x axis) are currently supported.

* The ``type`` parameter specifies the type of the boundary condition. Acceptable types are: ``fixed_wall``, ``outlet``, ``rotational``, ``translational`` and ``periodic``. The default boundary condition type is ``fixed_wall``.

* The ``periodic direction`` parameter specifies the perpendicular axis to the periodic boundaries.

* The ``rotational speed`` parameter defines the rotational speed of the specified boundary.  

* The ``rotational vector`` parameter specifies the rotational vector in `x`, `y`, and `z` directions.

* The ``point on rotational vector`` parameter specifies a point `x, y, z` on the rotating axis.

* The ``speed`` parameter defines the translational speed of the specified boundary.
