Boundary Conditions
-------------------
In this subsection, the boundary conditions of the DEM simulation are defined. First of all, the ``number of boundary conditions`` is specified. Then for each boundary condition, its information is defined. There are four boundary types: ``fixed_wall``, ``outlet``, ``rotational`` (around the center), and ``translational``. For ``rotational`` motion, ``rotational speed`` and ``rotational vector`` are required, while for ``translational`` motion, the ``speed`` should be defined in each direction.

.. note::
    The default boundary condition type is ``fixed_wall``. If the ``outlet`` condition is chosen for a boundary, particles can leave the simulation domain via this outlet.

.. code-block:: text

 subsection DEM boundary conditions
  # Total number of boundary motion
  set number of boundary conditions         		= 1

    # For each motion, we need a separate subsection
    subsection boundary condition 0

        # ID of boundary
	set boundary id					= 0

        # Boundary type
        # Choices are fixed_wall|outlet|rotational|translational
        set type              				= rotational

        # Rotational speed magnitude
	set rotational speed				= 2.5

        # Rotational vector
	set rotational vector x				= 1
	set rotational vector y				= 0
	set rotational vector z				= 0
    end

  # OR for translational motion
    subsection boundary condition 0

  # ID of moving boundary
	set boundary id	 				= 4

  # Motion type
        set type              				= translational

  # Speed in each direction
	set speed x					= 0.15
	set speed y					= 0
	set speed z					= 0
    end
 end

