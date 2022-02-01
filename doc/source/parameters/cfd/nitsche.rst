Nitsche
---------

Lethe can also run simulations using the Nitsche immersed boundary method. 

.. code-block:: text

  subsection nitsche
    set beta = 10
    set verbosity = verbose

    set number of solids = 1

    subsection nitsche solid 0
      set calculate forces on solid = false
      set solid force name = force_solid
      set calculate torques on solid = true
      set solid torque name = torque_solid
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
      subsection center of rotation
        set x = 0
        set y = 0
      end
      set number of quadrature points = 2
      set particles sub iterations = 5
      set stop if particles lost = true
      set enable particles motion = false
    end
  end

* ``beta`` is a parameter needed to apply the immersed boundary conditions on the fluid domain, its value is normally between 1 and 1000. ``beta = 0`` (the solid has no influence on the flow) can be used for debugging purposes. A classical value is ``beta = 10``.
* ``verbosity`` enables printing Nitsche intermediate results to the terminal.
* ``number of solids`` specifies the number of Nitsche solids in the simulation. Each solid will then correspond to a ``subsection nitsche solid``.
* ``subsection nitsche solid 0`` defines a solid object on which the Nitsche immersed boundary is applied. Multiple solids can be added in the same fashion (add ``subsection nitsche solid 1`` etc.).
* ``calculate forces on solid`` enables forces calculation on the immersed geometry, written in the output file named ``solid force name``, with the solid index automatically added at the end.
* ``calculate torques on solid`` enables torques calculation, written in the file in the output file named ``solid torque name``, with the solid index automatically added at the end. 
* ``subsection mesh`` defines the solid mesh used to apply Nitsche immersed boundary. The syntax is the same as that of the mesh subsection, see :doc:`mesh` for more details.

.. note::
  If ``set type = gmsh`` and a simplex mesh is given, do not forget to ``set simplex = true`` (it is ``false`` by default)

* ``subsection solid velocity`` defines the velocity of the solid mesh. This velocity is defined by a ``Function expression`` and can depend on both space and time.
* ``subsection center of rotation`` sets coordinates (x, y) of the center of the rotation for torque calculation. Default center of rotation is (0, 0). Add coordinate z for 3D simulations.
* ``number of quadrature points`` sets the number of Nitsche (quadrature) points to insert in a 1D cell. The number of inserted points will be higher for higher dimensions. Increasing this number will lead to a higher points density inside the solid.
* ``particles sub iterations`` splits the particle time-stepping into ``n`` sub time steps. This enables the particles to move less per iteration and makes the ``sort_particles_into_cells_and_subdomain()`` routine, which is used to locate the cells in which the particles reside, significantly faster. 
* ``stop if particles lost`` enables stopping the simulation if Nitsche particles have been lost. If ``false``, the simulation will continue. To prevent particle loss, try increasing the ``particles sub iterations``.
* ``enable particles motion`` must be set to ``true`` if the immersed boundary moves. If the boundary is static, for example a rotating cylinder, the shape does not have to move within the fluid and this option can be set to ``false``. This saves significant computational time.
* ``subsection test`` prints the positions of the particles in the terminal at the end of the simulation. This is mostly used for testing purposes.
