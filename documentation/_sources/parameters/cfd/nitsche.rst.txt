=========================
Nitsche Immersed Boundary
=========================

These parameters are used for simulations using the Nitsche immersed boundary method. Nitsche immersed boundary method works by forcing the fluid at the location of the Gauss points of the solid triangulation in order to apply the ``noslip`` boundary condition within the solid object.

.. seealso::
	For further understanding about the numerical method used and advanced parameters, the interested reader is referred to `this article <https://doi.org/10.1016/j.jcp.2023.112189>`_.

.. warning::
	The ``lethe-fluid-nitsche`` solver must be used for the Nitsche solid to be accounted for.

.. code-block:: text

  subsection nitsche
    #---------------------------------------------------
    # General parameters
    #---------------------------------------------------
    set verbosity        = quiet
    set number of solids = 0

    #---------------------------------------------------
    # Solid parameters
    #---------------------------------------------------
    subsection nitsche solid 0
      # Solid Mesh
      #=================================================
      subsection mesh
        set type               = dealii
        set grid type          = hyper_ball
        set grid arguments     = 0, 0 : 0.25 : true
        set initial refinement = 5
        set simplex            = false
      end

      # Velocity boundary condition parameters
      #=================================================
      # Penalization term for Nitsche method, on the velocity
      set beta                           = 10

      # Solid velocity
      subsection solid velocity
        # Default values in 2D
        set Function expression = 0 ; 0
        # in 3D: set Function expression = 0 ; 0 ; 0
      end

      # Condition on the motion of particles
      set enable particles motion        = false

      # Number of sub iterations for the motion of the particles
      set particles sub iterations       = 1

      # Heat boundary condition parameters
      #=================================================
      set enable heat boundary condition = false

      # Penalization term for Nitsche method, on the energy equation
      set beta heat                      = 10

      # Solid temperature
      subsection solid temperature
        set Function expression = 0
      end

      # Post-processing parameters
      #=================================================
      # Force and torque calculation on solid
      set calculate force on solid       = false
      set solid force name               = force_solid
      set calculate torque on solid      = false
      set solid torque name              = torque_solid

      # Center of rotation, used for torque calculation
      subsection center of rotation
        set x = 0
        set y = 0
        # in 3D: set z = 0
      end

      # Simulation control parameters
      #=================================================
      # Enable stopping the simulation if Nitsche particles have been lost
      set stop if particles lost         = true

      # Number of Nitsche (quadrature) points to insert in a 1D cell
      set number quadrature points       = 2
    end
  end

* ``verbosity``: controls if Nitsche intermediate results are printed to the terminal.

.. note::
	Even when ``verbosity = false``, Lethe produces additional files corresponding to the Nitsche immersed boundary:

	* the ``<output-name>_solid_triangulation_<id>.pvd``, corresponding to the mesh of the solid with index ``<id>`` ;
	* the ``<output-name>_solid_particles_<id>.pvd``, corresponding to the discrete particles inserted at the Gauss points of the solid triangulation, for the index ``<id>`` .

	The solid particles enable the Nitsche restriction visualization, while the solid triangulation is used for animation purposes.

* ``number of solids``: number of Nitsche solids in the simulation.

.. warning::
	The number of solids must be specified explicitly. This is often a source of error.

.. note::
	Each solid will then correspond to a ``subsection nitsche solid``.

* ``subsection nitsche solid 0``: defines a solid object, with index ``0``, on which the Nitsche immersed boundary is applied. Multiple solids can be added in the same fashion (``subsection nitsche solid 1`` etc.).

* ``subsection mesh``: defines the solid mesh used to apply Nitsche immersed boundary. The syntax is the same as that of the mesh subsection, see :doc:`mesh` for more details.

.. warning::
	If ``set type = gmsh`` and a simplex mesh is given, do not forget to ``set simplex = true`` (default value is ``false``)

.. tip::
	The solid mesh should have a characteristic size of the same order as the fluid dynamics mesh. Using a finer mesh will not cause any problem, but will increase the computational cost without benefits. 

* ``beta``: penalization term, which controls the intensity of the Nitsche method application on the velocity of the fluid region. Higher values of ``beta`` lead to stiffer problems but prevent the fluid from penetrating the solid.

.. tip::
	For flows with Reynolds numbers :math:`Re > 1`, we found that setting ``beta = 10`` (default value) leads to satisfactory results. 

	For ``beta = 0``, the solid has no influence on the flow: this value can be used for debugging purposes.
	
	In case of a static solid, ``beta`` parameter has to be greatly increased, up to ``100`` or ``1000``, to prevent the fluid moving through the solid. For highly viscous flows, even higher values of ``beta`` could be used to compensate for the larger shear stresses acting on the immersed solid.

* ``subsection solid velocity``: defines the velocity of the solid mesh. This velocity is defined by a ``Function expression`` and can depend on both space and time.

.. admonition:: Examples of solid velocity ``Function expression``:

	``set Function expression = 2 ; 0 ; 0``: 3D simulation, the solid is translating along the x-axis, with a norm of :math:`2`.

	``set Function expression = 3 ; -4``: 2D simulation, the solid is translating along a composition of the x and y-axes, with a norm of :math:`\sqrt(3^2+(-4)^2) = 5`.

	``set Function expression = -y ; x``: 2D simulation, the solid is rotating in the anti-clockwise direction around the origin, with a tangential velocity of norm :math:`1`.

.. tip::
	The unit of the solid velocity value depends on the units of the mesh: if the mesh is build with the meter as the base unit, the velocity will be in :math:`m/s`.
	
* ``enable particles motion``: controls if the immersed boundary moves within the fluid domain.

.. tip ::
	For a rotating cylinder, the ``Nitsche solid`` rotates but the boundary location does not change. For such static boundaries, the shape does not have to move within the fluid and this option can be set to ``false``. This saves significant computational time.

.. warning ::
	When the ``solid velocity`` leads to a motion of the solid, use ``enable particles motion = true``.

* ``particles sub iterations``: number of sub iterations for the motion of the particles. 

.. tip ::
	When ``set particles sub iterations = 1`` (default value), there is no sub iteration: the motion of the particle is solved at each ``time step`` (see :doc:`simulation_control`). 

	In case of particle loss, this parameter can be increased (``set particles sub iterations = 5`` is a good start value) to ensure that particles are always located efficiently as they move through the cell. This increases the computational cost, but not as much as lowering the ``time step`` (in :doc:`simulation_control`) would.

	Generally, it is a good practice to have sufficient ``particles sub iterations`` so as to ensure that particles do not move more than half a cell during a particle sub iteration.

* ``enable heat boundary condition``: controls if a heat boundary condition is imposed on the Nitsche immersed boundary. Use to attribute a given temperature to the Nitsche solid (defined in ``subsection solid temperature``).

* ``beta heat``: penalization term on the heat equation, which controls the intensity of the Nitsche method application on the temperature of the fluid region. Higher values of ``beta`` forces the fluid near the solid to have a temperature matching the one of the solid (defined in ``subsection solid temperature``).
	
* ``subsection solid temperature``: defines the temperature of the solid mesh. This temperature is defined by a ``Function expression`` and can depend on both space and time. This parameter is used only if:

  * ``enable heat boundary condition = true``, and
  * ``heat transfert = true`` in :doc:`multiphysics` subsection.

* ``calculate force on solid``: controls if force calculation on the immersed geometry is enabled. If set to ``true``, forces will written in the output file named ``solid force name``, with the solid index automatically added at the end.

* ``calculate torque on solid``: controls if torque calculation on the immersed geometry is enabled. If set to ``true``, torques will be written in the file in the output file named ``solid torque name``, with the solid index automatically added at the end. 

* ``subsection center of rotation``: :math:`(x, y)` coordinates of the center of the rotation, used for torque calculation. Default center of rotation is (0, 0). Add ``set z`` for 3D simulations.

* ``stop if particles lost``: controls if the simulation is stopped when Nitsche particles have been lost. If ``false``, the simulation will continue. 

.. tip ::

	Particle loss can happen when particles move through multiple cells during a time step. This can be caused by a big ``time step`` (see :doc:`simulation_control`), a high fluid ``mesh refinement`` (see :doc:`mesh`), or a high CFL. To prevent particle loss, try increasing the number of ``particles sub iterations``.

* ``number quadrature points``: number of Nitsche (quadrature) points to insert in a 1D cell. The number of inserted points will be higher for higher dimensions. Increasing this number will lead to a higher points density inside the solid.

.. seealso::
	The Nitsche immersed boundary method is used in the examples:
	  * :doc:`../../examples/incompressible-flow/2d-taylor-couette-flow-nitsche/2d-taylor-couette-flow-nitsche`
	  * :doc:`../../examples/incompressible-flow/3d-nitsche-mixer-with-pbt-impeller/nitsche-mixer-with-pbt-impeller`


