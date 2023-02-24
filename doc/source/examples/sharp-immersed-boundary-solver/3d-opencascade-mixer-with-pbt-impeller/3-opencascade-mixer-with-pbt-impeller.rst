.. role:: raw-html(raw)
    :format: html

=======================================================================================
3D Mixer with pitched-blade turbine impeller using OpenCascade Sharp-immersed boundary
=======================================================================================

The mixing of stirred-tanks is a common chemical engineering problem that can be tackled through immersed boundary simulation. This example, a variation of :doc:`../../incompressible-flow/3d-nitsche-mixer-with-pbt-impeller/nitsche-mixer-with-pbt-impeller`, illustrates how the transient flow in a stirred-tank can be simulated by Lethe using the Sharp-Immersed Boundary formulation with a OpenCasacade shape from a step file.

.. seealso::
	The original example using Nitsche immersed boudaries: :doc:`../../incompressible-flow/3d-nitsche-mixer-with-pbt-impeller/nitsche-mixer-with-pbt-impeller`.

:raw-html:`<br />`

Features
----------------------------------
- Solvers: ``gls_sharp_navier_stokes_3d``
- Transient problem
- Rotating complex solid, defined by a step file using the OpenCascade shape, modeled with sharp immersed boundary

:raw-html:`<br />`

Files used in this example
----------------------------

* Parameter file: ``/examples/sharp-immersed-boundary-solver/3d-opencascade-mixer-with-pbt-impeller/mixer.prm``
* Step file: ``/examples/sharp-immersed-boundary-solver/3d-opencascade-mixer-with-pbt-impeller/impeller.step``

:raw-html:`<br />`

Description of the case
-----------------------

In this example, we simulate a mixer using a PBT impeller through the usage of a step file with the OpenCascade shape feature of the sharp immersed boundary solver.

Creation of the step file.
------------------------------------

The step file can be defined using any CAD tool available to the user. The step file must represent a solid. It's favorable to avoid step file that represent shells, composite of solides or compound of objects. 

.. tip::
	Use the union tool at your disposal to avoid issues with step files that are defined by a composite of solids. Most CAD software offers the possibility to define a solid from the union of multiple solids. Similarly, if the step file is only defined by a shell, it is usually possible to define a solid from that shell. If your CAD tool does not allow these operations, the FreeCAD software allows you to do these operations using the part toolbox. 


Definition of the shape and its motion
--------------------------------------

The section defining each parameter for the particles has certains requirements:

1. ``length ratio`` defined the length used to apply the immersed boundaries through interpolation. It should stay as low as possible, but above ``1``.
2. ``type`` and ``shape arguments`` are used to declare that the shape is a ``opencascade`` and that its data is located in ``impeller.step``.
3. ``integrate motion`` is set to ``false``. This way, the solid only moves according to the prescribed `orientation` and angular velocity `omega` (the alternative being the integration of particle movement from forces).

.. code-block:: text

  subsection particles
    set number of particles =1
    set stencil order = 1
    set refine mesh inside radius factor = 0.0
    set refine mesh outside radius factor = 1.1
    set length ratio = 3
    set initial refinement = 1
    set integrate motion = false
    set assemble Navier-Stokes inside particles =false

    subsection particle info 0
      subsection orientation
          set Function expression =-1*2*pi*t;pi/2;0
      end
      subsection omega
          set Function expression =-1*2*pi;0;0
      end
      set type       = opencascade
      set shape arguments = impeller.step
    end
  end

Additionnal information on the ``particles`` parameters can be found on :doc:`../../../parameters/sharp-immersed-boundary-solver/sharp-immersed-boundary-solver`.

Boundary conditions
-----------------------

Because the particles defined for the sharp solver are not divergence-free, it is necessary to have at least one boundary condition that is weakly imposed to ensure the system of equations is well-posed. For this purpose, a ``function weak`` type of boundary is used.
Two aspects need special consideration:

1. ``function weak`` is a variation of ``function``. It is used to weakly imposed a Dirichlet boundary condition, and it is necessary when using ``gls_sharp_navier_stokes_3d``.
2. ``beta`` has to be defined. It is a Nitsche penalization parameter that enforces more strongly the boundary condition when it increases (see :doc:`../../../parameters/cfd/nitsche`).

.. code-block:: text

  subsection boundary conditions
    set number = 3
    subsection bc 0
      set id   = 0
      set type = noslip
    end
    subsection bc 1
      set id   = 1
      set type = noslip
    end
    subsection bc 2
      set id   = 2
      set type              = function weak
          set beta = 1
          subsection u
              set Function expression = 0
          end
          subsection v
              set Function expression = 0
          end
          subsection w
              set Function expression = 0
          end
    end
  end

Results
--------

This example allows to reach similar results as the original example. .....
