==================================
Bunny Drill
==================================

This example simulates the drilling motion of a bunny within a bed of particles. It illustrates that the DEM module of Lethe can simulate complex moving objects and is a testament to our love of lagomorphs. Do not worry friend, no bunnies were hurt in the making of this example!


----------------------------------
Features
----------------------------------
- Solvers: ``lethe-particles``
- Floating walls
- `GMSH <https://gmsh.info/>`_ grids

----------------------------
Files Used in This Example
----------------------------

All files mentioned below are located in the example's folder (``examples/dem/3d-bunny-drill``).

- GMSH mesh of the bunny: ``bunny-low-poly.msh`` generated using the corresponding STL file ``bunny-low-poly.stl``
- Parameter file used to load the particles: ``bunny-drill-loading.prm``
- Parameter file used to simulate the bunny drill: ``bunny-drill.prm``

-----------------------
Description of the Case
-----------------------

This simulation consists of two stages: filling (0-2 s) and bunny-drilling (2-4.75 s) of particles.

--------------
Parameter File
--------------

Mesh
~~~~~

The mesh is a cylinder generated using the deal.II grid generator.

.. code-block:: text

  subsection mesh
    set type               = cylinder
    set grid type          = balanced
    set grid arguments     = 6 : 0.10 : 0.25
    set initial refinement = 2
  end

Insertion Info
~~~~~~~~~~~~~~~~~~~

An insertion box is defined inside at the top part of the cylinder. The insertion region is sufficiently high to ensure that there are no collisions between the particles bouncing on the bottom wall (or particles accumulated at the bottom) and the ones being inserted.

.. code-block:: text

  subsection insertion info
    set insertion method                               = non_uniform
    set inserted number of particles at each time step = 10000
    set insertion frequency                            = 100000
    set insertion box minimum x                        = 0.1
    set insertion box minimum y                        = -0.07
    set insertion box minimum z                        = -0.07
    set insertion box maximum x                        = 0.24
    set insertion box maximum y                        = 0.07
    set insertion box maximum z                        = 0.07
    set insertion distance threshold                   = 1.10
    set insertion random number range                  = 0.10
    set insertion random number seed                   = 19
    set insertion first direction                      = 0
    set insertion second direction                     = 2
    set insertion third direction                      = 1
  end


Lagrangian Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The total number of particles in this simulation is 8000. All particles have a diameter of 10 mm. The particles are relatively stiff, with a Young Modulus of :math:`E=10^7` to ensure that the overlap between the bunny drill and the particles does not become sufficiently large to allow particles to *jump through* the bunny.

.. code-block:: text

  subsection lagrangian physical properties
    set gx                       = -9.81
    set gy                       = 0
    set gz                       = 0
    set number of particle types = 1
    subsection particle type 0
      set size distribution type            = uniform
      set diameter                          = 0.01
      set number of particles               = 8000
      set density particles                 = 2560
      set young modulus particles           = 1e7
      set poisson ratio particles           = 0.3
      set restitution coefficient particles = 0.9
      set friction coefficient particles    = 0.2
      set rolling friction particles        = 0.3
    end
    set young modulus wall           = 1e7
    set poisson ratio wall           = 0.2
    set restitution coefficient wall = 0.9
    set friction coefficient wall    = 0.5
    set rolling friction wall        = 0.1
  end



Simulation Control (Loading)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The time end of the simulation is 2 seconds after all particles have been loaded.

.. code-block:: text

  subsection simulation control
    set time step         = 2e-6
    set time end          = 2
    set log frequency     = 1000
    set output frequency  = 1000
    set output path       = ./output/
    set output boundaries = true
  end

Simulation Control (Drilling)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The time end of the simulation is 4.75 seconds after which the bunny has done one back-and-forth drilling motion.

.. code-block:: text

  subsection simulation control
    set time step         = 2e-6
    set time end          = 4.75
    set log frequency     = 1000
    set output frequency  = 1000
    set output path       = ./output/
    set output boundaries = true
  end



Solid Objects (Drilling)
~~~~~~~~~~~~~~~~~~~~~~~~~

The bunny is defined using the solid objects feature of Lethe. The surface mesh of the bunny is a GMSH file. The translational of the velocity is defined to have a periodic motion along the axis of the cylinder and the bunny is rotating at a constant angular velocity. This complex drilling motion is fully parametrized from the input file using the function parser functionalities of the translational and the angular velocity of the solid object.

.. code-block:: text

  subsection solid objects
    set number of solids = 1
    subsection solid object 0
      subsection mesh
        set type                   = gmsh
        set file name              = bunny-low-poly.msh
        set simplex                = true
        set initial rotation axis  = 0, 1, 0
        set initial rotation angle = 1.5708 # pi/2
        set initial translation    = 0.05, 0, 0.035
      end
      subsection translational velocity
        set Function expression = -0.27*sin(0.8*3.1416*(t-2)) ; 0 ; 0
      end
      subsection angular velocity
        set Function expression = 31.42 ; 0 ; 0
      end
    end
  end


----------------------
Running the Simulation
----------------------
The loading can be simulated using the following command:

.. code-block:: text
  :class: copy-button

  mpirun -np 8 lethe-particles bunny-drill-loading.prm

Whereas the drilling is launched after the loading using:

.. code-block:: text
  :class: copy-button

  mpirun -np 8 lethe-particles bunny-drill.prm


-------
Results
-------
As seen in the following two animations, the bunny drills into the particles which generates a complex motion within the granular matter. There is not much more to say here, it is a bunny drill.

The first animation displays the drill with the entirety of the particles. It is difficult to see the dynamics of the mighty bunny within these circumstances.

.. raw:: html

    <iframe width="500" height="600" src="https://www.youtube.com/embed/GI_jfsO0ZeM" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

The following animation displays the drill with half of the particles clipped. Here we can clearly see the bunny in action.

.. raw:: html

    <iframe width="500" height="600" src="https://www.youtube.com/embed/VcJ_nt9iNmA" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

----------------------------
Possibilities for Extension
----------------------------

- Use finer particles to see if the drilling dynamics are affected by the particle size.
- Use an STL of an alternative animal. Although we believe lagomorphs are amazing, we are also fans of mustelidae (e.g., otters) and chinchillidae (e.g., chinchillas or, even better, viscachas). Feel free to replace the drill with your favorite animal and to send us your animation to lethe.cfd@gmail.com.


 