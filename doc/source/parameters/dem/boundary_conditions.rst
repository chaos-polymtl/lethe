===================
Boundary Conditions
===================

In this subsection, the boundary conditions of the DEM simulation are defined. First of all, the ``number of boundary conditions`` is specified. Then for each boundary condition, its information is defined. There are five boundary types: ``fixed_wall``, ``outlet``, ``rotational`` (around the center), ``translational``, and ``periodic``. For ``rotational`` motion, ``rotational speed`` and ``rotational vector`` are required, while for ``translational`` motion, the ``speed`` should be defined in each direction. For ``periodic`` boundaries, ``periodic id 0``, ``periodic id 1`` and ``periodic direction`` are required.

``fixed_wall`` is a static wall, and particles collide with these static walls upon reaching the wall. The only way to move these walls is to move the entire triangulation. If the ``outlet`` condition is chosen for a boundary, particles can leave the simulation domain via this outlet. Using ``rotational`` and ``translational`` boundary conditions, exerts imaginary rotational and translational velocities to that boundary. In other words, the boundary does not move, but the particles that have collisions with these walls receive a rotational or translational velocity from the wall. This feature is used in the rotating drum example.

In multiphysic DEM (``solver type = dem_mp`` in the :doc:`model_parameters` subsection), each wall can also be given a ``thermal boundary type``, independently of its motion. A wall is ``adiabatic`` by default: it does not exchange heat with the particles. An ``isothermal`` wall has a temperature imposed by the user, and it exchanges heat with the particles in contact with it. 

.. code-block:: text

  subsection DEM boundary conditions
    # Total number of boundary motion
    set number of boundary conditions = 3

    # For each boundary condition, we need a separate subsection
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

      # Thermal boundary type (multiphysic DEM only)
      # Choices are adiabatic|isothermal
      set thermal boundary type = adiabatic

      # Temperature of the wall
      subsection wall temperature
        set Function expression = 0
      end
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

    # OR for periodic boundaries
    subsection boundary condition 2
      # Boundary type
      set type               = periodic

      # ID of principal boundary
      set periodic id 0      = 3

      # ID of associated periodic boundary
      set periodic id 1      = 2

      # Direction normal to the periodic boundary faces (x=0, y=1, z=2)
      set periodic direction = 1
  end

* The ``number of boundary conditions`` parameter defines the number of desired boundary conditions to be specified. Note that if a boundary condition ``type`` is not defined explicitly, Lethe defines it as a fixed static wall.

* For each boundary condition, we have to define a separate subsection. In the sample parameter list above, the ``number of boundary conditions`` is equal to 3. Hence, we need to define three subsections (``subsection boundary condition 0``, ``subsection boundary condition 1`` and ``subsection boundary condition 2```).

* The ``boundary id`` parameter specifies the boundary ID for which the boundary condition should be applied. Periodic boundaries are an exception.

* The ``type`` parameter specifies the type of the boundary condition. Acceptable types are: ``fixed_wall``, ``outlet``, ``rotational``, ``translational`` and ``periodic``. The default boundary condition type is ``fixed_wall``.

* The ``periodic id 0`` and ``periodic id 1`` parameters specify the periodic boundaries ID for which the periodic boundary condition should be applied. By convention, ``periodic id 0`` should correspond to the boundary ID for which the boundary is further along the coordinate axis (this is the principal boundary).

.. note::
        Only periodic boundaries which have co-linear normal vectors which align along an axis of the problem (e.g., x axis) are currently supported. Multiple simultaneous periodic directions are supported.

* The ``periodic direction`` parameter specifies the perpendicular axis to a pair of periodic boundaries (``0`` is the `x` axis, ``1`` is the `y` and ``2`` is the `z`).

* The ``rotational speed`` parameter defines the rotational speed of the specified boundary.  

* The ``rotational vector`` parameter specifies the rotational vector in `x`, `y`, and `z` directions.

* The ``point on rotational vector`` parameter specifies a point `x, y, z` on the rotating axis.

* The ``speed`` parameter defines the translational speed of the specified boundary.

* The ``thermal boundary type`` parameter defines whether the boundary is ``adiabatic`` or ``isothermal`` in a multiphysic DEM simulation. The default is ``adiabatic``: no heat is exchanged between the particles and the wall. Only ``fixed_wall``, ``translational`` and ``rotational`` boundaries can be ``isothermal``, since particles cannot be in contact with ``outlet`` or ``periodic`` boundaries. An ``isothermal`` boundary can only be used when the ``solver type`` is ``dem_mp``.

* In the subsection ``wall temperature``, we define the temperature of an ``isothermal`` boundary as a function of space and time (:math:`x`, :math:`y`, :math:`z` and :math:`t`). The function is evaluated at the contact point between each particle and the wall, which is the projection of the center of the particle on the wall, so the temperature of the wall can be non-uniform. For example, ``if(y < 0, 80, 20)`` heats the part of the wall below :math:`y = 0` only. 

.. note::
    The temperature field of a ``rotational`` or ``translational`` wall is defined in the fixed frame of reference of the mesh, since the mesh does not move: the wall only transmits its velocity to the particles. For a temperature pattern that moves with the wall, write the function in the frame of the wall. For instance, for a drum rotating at the ``rotational speed`` :math:`\omega` around the :math:`x` axis, use the angle :math:`\operatorname{atan2}(z, y) - \omega t`.

.. note::
    The temperature of an ``isothermal`` wall is imposed, and the wall behaves like a heat reservoir of infinite heat capacity. The heat transfer rate between a particle :math:`i` and the wall :math:`w` is :math:`Q_{iw} = H_{iw} (T_w - T_i)`, where the thermal conductance :math:`H_{iw}` is computed with the particle-wall model described in the `theory guide <../../theory/multiphase/dem/dem.html#particle-wall-resistances>`_. The model is the same as for the ``isothermal`` solid objects, and it uses the wall properties defined in the :doc:`lagrangian_physical_properties` subsection (``thermal conductivity wall``, ``microhardness wall``, etc.).
