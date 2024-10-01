==================================
Granular Mixer
==================================

This example simulates the packing and mixing of particles in a mixer with a pitched-blade impeller. It is recommended to visit `DEM parameters <../../../parameters/dem/dem.html>`_ for more detailed information on the concepts and physical meanings of the parameters in Lethe-DEM.


----------------------------------
Features
----------------------------------

- Solvers: ``lethe-particles``
- Floating mesh
- `GMSH <https://gmsh.info/>`_ grids
- Bidispersed particles (same size and properties, but different types)


----------------------------
Files Used in This Example
----------------------------

- Parameter file: ``/examples/dem/3d-granular-mixer/granular-mixer.prm``


-----------------------
Description of the Case
-----------------------

This simulation consists of two stages: packing (0-0.5 s) and mixing (0.5-5 s) of particles. There are two types of particles in this simulation (bidispersed system), that are inserted on top of each other during the packing stage. The size and properties of the two particle types are the same, we only need to define two particle types to make the visualization easier during post-processing. At :math:`t=0.5` s, the pitched-blade impeller starts rotating with an angular velocity of 6 rad/s and mixes the particles.


--------------
Parameter File
--------------

Mesh
~~~~~

The background mesh (mixer body) is created using deal.II ``subdivided_cylinder``.

.. code-block:: text

    subsection mesh
      set type                = dealii
      set grid type           = subdivided_cylinder
      set grid arguments      = 2 : 0.05 : 0.055
      set initial refinement  = 3
    end

Lagrangian Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As mentioned earlier, there are two types of particles with the same size and properties.

.. code-block:: text

    subsection lagrangian physical properties
      set g                                   = -9.81, 0.0, 0.0
      set number of particle types            = 2
      subsection particle type 0
        set size distribution type              = uniform
        set diameter                            = 0.0015
        set number of particles                 = 23500
        set density particles                   = 1500
        set young modulus particles         	  = 1e6
        set poisson ratio particles             = 0.5
        set restitution coefficient particles   = 0.5
        set friction coefficient particles      = 0.5
      end
      subsection particle type 1
        set size distribution type              = uniform
        set diameter                            = 0.0015
        set number of particles                 = 23500
        set density particles                   = 1500
        set young modulus particles             = 1e6
        set poisson ratio particles             = 0.5
        set restitution coefficient particles   = 0.5
        set friction coefficient particles      = 0.5
      end
      set young modulus wall                  = 1e6
      set poisson ratio wall                  = 0.5
      set restitution coefficient wall        = 0.5
      set friction coefficient wall           = 0.5
    end


Solid Objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this subsection, the floating meshes are defined. We can use deal.II or Gmsh to create the floating meshes. At the moment, solid objects in Lethe have to be defined using triangular (simplex) meshes. Only triangular 2D meshes of 3D surfaces in the ``lethe-particles`` solver are presently supported. Quadrilateral 2D meshes of 3D surfaces and 1D mesh of 2D surfaces are not supported at the moment. For each floating mesh, we need to specify a ``translational velocity``, an ``angular velocity``, and a ``center of rotation``. In this example, we only need an angular motion of the impeller. Note that the ``center of rotation`` of the impeller is at 0, 0, 0.

.. code-block:: text

    subsection solid objects
      subsection solid surfaces
        set number of solids = 1
        subsection solid object 0
          subsection mesh
            set type               = gmsh
            set file name          = pitched-blade-impeller.msh
            set simplex            = true
            set initial refinement = 0
          end
    
          subsection translational velocity
            set Function expression = 0 ; 0 ; 0
          end
          subsection angular velocity
            set Function expression = if(t>0.5,6,0) ; 0 ; 0
          end
          set center of rotation = 0, 0, 0
        end
      end
    end


----------------------
Running the Simulation
----------------------
This simulation can be launched by (in parallel mode on 8 processes):

.. code-block:: text
  :class: copy-button

  mpirun -np 8 lethe-particles granular-mixer.prm

.. warning::
	This example takes approximately 2 hours on 8 cores.


---------
Results
---------

Animation of the granular mixing simulation:

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/ms-gAyZcOXk" frameborder="0" allowfullscreen></iframe>


-----------------------------
Possibility for Extension
-----------------------------

The same simulation can be carried out with particles of different sizes and properties to study segregation.
