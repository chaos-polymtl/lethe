Nitsche
---------

Lethe can also run simulations using the Nitsche immersed boundary method. 

.. code-block:: text

  subsection nitsche
    set beta = 10
    set verbosity = verbose
    set calculate forces on solid = false
    set solid force name = force_solid
    set calculate torques on solid = true
    set solid torque name = torque_solid

    set number of solids = 1

    subsection nitsche solid 0
	  subsection mesh
	  	set type = dealii
	  	set grid type = hyper_ball
	  	set grid arguments = 0, 0 : 0.25 : true
	  	set initial refinement = 5
	  	set simplex = false
	  end
	  subsection solid velocity
	  	set Function expression = -y ; x
	  end
          set particles sub iterations = 5
          set stop if particles lost = true
    end

* ``beta`` is a parameter needed to apply the immersed boundary conditions on the fluid domain, its value is normally between 1 and 1000. ``beta = 0`` (the solid has no influence on the flow) can be used for debugging purposes. A classical value is ``beta = 10``.
* ``verbosity`` enables printing Nitsche intermediate results to the terminal.
* ``calculate forces on solid`` enables forces calculation on the immersed geometry, written in the output file named ``solid force name``. 
* ``calculate torques on solid`` enables torques calculation, written in the file in the output file named ``solid torque name``. 
* ``number of solids`` specifies the number of Nitsche solids in the simulation. Each solid will then correspond to a ``subsection nitsche solid``.
* ``subsection nitsche solid 0`` defines a solid object on which the Nitsche immersed boundary is applied. Multiple solids can be added in the same fashion (add ``subsection nitsche solid 1`` etc.).
* ``subsection mesh`` defines the solid mesh used to apply Nitsche immersed boundary. The syntax is the same as that of the mesh subsection, see :doc:`mesh` for more details.

.. note::
  If ``set type = gmsh`` and a simplex mesh is given, do not forget to ``set simplex = true`` (it is ``false`` by default)

* ``subsection solid velocity`` defines the velocity of the solid mesh. This velocity is defined by a ``Function  expression`` and can depend on both space and time.
* ``enable particles motion`` must be set to ``true`` if the immersed boundary moves. If the boundary is static, for example a rotating cylinder, the shape does not have to move within the fluid and this option can be set to ``false``. This saves significant computational time.
* ``particles sub iterations`` splits the particle time-stepping into ``n`` sub time steps. This enables the particles to move less per iteration and makes the ``sort_particles_into_cells_and_subdomain()`` routine, which is used to locate the cells in which the particles reside, significantly faster. 
* ``stop if particles lost`` enables stopping the simulation if Nitsche particles have been lost. If ``false``, the simulation will continue. To prevent particle loss, try increasing the ``particles sub iterations``.
* ``subsection test`` prints the positions of the particles in the terminal at the end of the simulation. This is mostly used for testing purposes.
