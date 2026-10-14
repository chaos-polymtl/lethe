==========================================
Couette Flow with Immersed Particles
==========================================

This example tests Einstein's dilute viscosity on a Couette flow with rigid particles, in two and three dimensions. The wall shear force is compared with :math:`\mu(1 + \tfrac{5}{2}\phi)`. The comparison shows where that formula is sufficient, and where the measured viscosity falls below it.

----------------------------------
Features
----------------------------------

- Solver: ``lethe-fluid-sharp``
- Couette flow in two and three dimensions
- Sharp-edge immersed boundary methods, at fixed centers

- Solid fraction varied by the particle radiu
---------------------------
Files Used in This Example
---------------------------

All files are in ``examples/sharp-immersed-boundary/couette-2d-3d``.

- ``couette-2d-corr.prm``
- ``couette-3d.prm`` 
- ``epsilon_2d.py``
- ``epsilon_3d.py``


-----------------------
Description of the Case
-----------------------

A Newtonian fluid is sheared between a fixed plate and a plate moving at constant speed. Without particles the velocity profile is linear and the plate force is known. Rigid particles are then placed in the gap. Each particle is given the velocity of that linear profile at its center, so the particle follows the fluid and does not slip. The plates feel an extra force. Einstein's theory attributes it to the solid fraction :math:`\phi` [#einstein1906]_, [#einstein1911]_.

The relation is written with the intrinsic viscosity :math:`[\eta]` [#wikiintrinsic]_:

.. math::
    \frac{\mu_{eff}}{\mu} = 1 + [\eta]\,\phi, \qquad [\eta] = \frac{5}{2}.

This is a dilute result: particles interact with the flow, not with each other. It is expected to be accurate at a few percent, and to drift once :math:`\phi` is large. The coefficient :math:`5/2` was derived for spheres. The two-dimensional runs still use it, so that both geometries are judged by the same formula. The force on the moving plate is negative, because the fluid resists the motion. The curves use the absolute value.

Two-Dimensional Channel
~~~~~~~~~~~~~~~~~~~~~~~

The channel is :math:`50 \times 20`. The bottom plate is fixed and the top plate moves at speed 1. The streamwise direction is periodic. With :math:`\mu = 10` the shear rate is :math:`1/20 = 0.05`, and the empty-channel force is

.. math::
    F_{x,0} = 10 \times 0.05 \times 50 = 25.

Seven discs share one set of centers. Only the radius changes, from the reference value 1 up to the value that produces :math:`\phi = 15\,\%`.

.. image:: images/couette-2d-7p.png
    :alt: Two-dimensional Couette flow with seven immersed particles
    :align: center
    :width: 600px

========  ==========
Particle  Position
========  ==========
0         (7, 5)
1         (18, 14)
2         (28, 7)
3         (38, 12)
4         (12, 16)
5         (33, 4)
6         (44, 9)
========  ==========

Three-Dimensional Channel
~~~~~~~~~~~~~~~~~~~~~~~~~

The box is :math:`75 \times 30 \times 20`. The gap is 30 and the top plate moves at 1.5, so the shear rate is again 0.05. The length 75 contains the ten centers, from :math:`x = 7` to :math:`x = 68`. The depth 20 is periodic, and every sphere sits at :math:`z = 10`. The empty-channel force is

.. math::
    F_{x,0} = 10 \times 0.05 \times 75 \times 20 = 750.

Ten spheres of radius 1 fill only 0.093 % of the volume. The same centers cannot exceed about 1 % without overlapping, so the radius sweep is limited to 0.2 %, 0.5 % and 0.8 %. Those three runs use weak imposition of the plates.

----------------------
Parameter File
----------------------

``couette-2d-corr.prm`` is a steady solve at :math:`\mu = 10`. ``epsilon_2d.py`` changes only the radius, the output name and the force name.

Simulation Control
~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection simulation control
      set method           = steady
      set output name      = couette_ex1_b_particles_eps0
      set output path      = out_2d/
    end

Mesh
~~~~

The two-dimensional grid is a rectangle of :math:`50 \times 20`, colored so that each side receives a boundary id.

.. code-block:: text

    subsection mesh
      set type               = dealii
      set grid type          = subdivided_hyper_rectangle
      set grid arguments     = 25,10 : 0,0 : 50,20 : true
      set initial refinement = 2
    end

The three-dimensional grid uses the same construction on :math:`75 \times 30 \times 20`.

.. code-block:: text

    subsection mesh
      set type               = dealii
      set grid type          = subdivided_hyper_rectangle
      set grid arguments     = 25,10,10 : 0,0,0 : 75,30,20 : true
      set initial refinement = 2
    end

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

In two dimensions the streamwise direction is periodic. The plates are imposed weakly, with ``beta = 10``.

.. code-block:: text

    subsection boundary conditions
      set number = 3
      subsection bc 0
        set id                 = 0
        set type               = periodic
        set periodic id        = 1
        set periodic direction = 0
      end
      subsection bc 1
        set id   = 2
        set type = function weak
        subsection u
          set Function expression = 0
        end
        subsection v
          set Function expression = 0
        end
        set beta = 10
      end
      subsection bc 2
        set id   = 3
        set type = function weak
        subsection u
          set Function expression = 1
        end
        subsection v
          set Function expression = 0
        end
        set beta = 10
      end
    end

In three dimensions, ``couette-3d.prm`` adds a periodic pair in :math:`z` (ids 4 and 5, direction 2). The bottom plate (id 2) is fixed and the top plate (id 3) moves at 1.5. That file imposes both plates with ``function``. The runs written by ``epsilon_3d.py`` replace it by ``function weak`` and ``beta = 10``, as in the two-dimensional file.

Particles
~~~~~~~~~

The Navier–Stokes equations are not assembled inside the discs. Motion is prescribed, and each disc receives the linear profile :math:`(y/20,\, 0)`. The six other particles repeat this block at the centers listed above.

.. code-block:: text

    subsection particles
      set assemble Navier-Stokes inside particles = false
      set number of particles                     = 7
      subsection extrapolation function
        set length ratio   = 2
        set stencil degree = 1
      end
      subsection local mesh refinement
        set initial refinement                = 4
        set refine mesh inside radius factor  = 0.7
        set refine mesh outside radius factor = 1.3
      end
      subsection particle info 0
        set type             = sphere
        set shape arguments  = 1
        set integrate motion = false
        subsection position
          set Function expression = 7; 5
        end
        subsection velocity
          set Function expression = y/20; 0
        end
      end
    end

``couette-3d.prm`` contains ten spheres of radius 1. The velocity is the constant :math:`y/20` at the particle center (here :math:`5/20 = 0.25`). The runs written by ``epsilon_3d.py`` set ``assemble Navier-Stokes inside particles = false``.

.. code-block:: text

    subsection particles
      set number of particles = 10
      subsection extrapolation function
        set length ratio   = 2
        set stencil degree = 2
      end
      subsection local mesh refinement
        set initial refinement                = 2
        set refine mesh inside radius factor  = 0.7
        set refine mesh outside radius factor = 1.2
      end
      subsection particle info 0
        set type             = sphere
        set shape arguments  = 1
        set integrate motion = false
        subsection position
          set Function expression = 7;5;10
        end
        subsection velocity
          set Function expression = 0.25;0;0
        end
      end
    end

Mesh Adaptation
~~~~~~~~~~~~~~~

.. code-block:: text

    subsection mesh adaptation
      set type = adaptive
    end

``couette-3d.prm`` sets this type to ``none``.

Forces
~~~~~~

Lethe writes one force file per boundary id. ``.02.dat`` is the bottom plate (id 2) and ``.03.dat`` is the top plate (id 3).

.. code-block:: text

    subsection forces
      set calculate force       = true
      set calculate torque      = false
      set force name            = force_7p_eps0
      set calculation frequency = 1
      set output frequency      = 1
    end

The three-dimensional file uses ``force_10p_3d``.

---------------------------
Running the Simulation
---------------------------

The centers stay fixed. A script reads the base parameter file and, for each solid fraction :math:`\phi`, computes the matching radius and writes one parameter file. The output name, the force name and the output path change with :math:`\phi`.

``epsilon_2d.py`` starts from ``couette-2d-corr.prm``. Seven discs give

.. math::
    \phi = \frac{7\pi R^{2}}{50 \times 20}.

The script builds the cases :math:`\phi = 2\,\%`, :math:`4\,\%`, :math:`5\,\%` and :math:`15\,\%`.

``epsilon_3d.py`` starts from ``couette-3d.prm``. Ten spheres give

.. math::
    \phi = \frac{10 \times \tfrac{4}{3}\pi R^{3}}{75 \times 30 \times 20}.

The same centers overlap above about :math:`1\,\%`, so the script stops at :math:`0.2\,\%`, :math:`0.5\,\%` and :math:`0.8\,\%`. Before those runs, the plates are set to ``function weak`` with ``beta = 10``, and the Navier–Stokes equations are not assembled inside the spheres.

.. code-block:: text

    python3 epsilon_2d.py
    python3 epsilon_3d.py
    mpirun -np 4 lethe-fluid-sharp couette_ref.prm
    mpirun -np 4 lethe-fluid-sharp couette-2d-corr.prm
    python3 calc_ratio.py

Each file written by the script is then run with ``lethe-fluid-sharp`` from its own directory. A three-dimensional case takes about one minute on 4 MPI ranks. The two-dimensional cases take longer because of the local refinement.

-----------------------
Results and Discussion
-----------------------

The two-dimensional curve compares the measured viscosity with :math:`\mu(1 + \tfrac{5}{2}\phi)`.

.. figure:: images/viscosity-2d.png
    :alt: Effective viscosity against solid fraction in two dimensions
    :align: center
    :width: 600px

    Effective viscosity of the two-dimensional sweep. The curve is :math:`\mu(1+\tfrac{5}{2}\phi)`.

The points stay on the curve up to 5 %, where the error is still 0.75 %. At 15 % the simulation falls below the curve, 13.626 against 13.750. Five percent is the highest fraction at which the formula can be used as is.

The three-dimensional curve is the same comparison, on the viscosity :math:`10\,|F_x|/750`. The file ``couette-3d.prm`` is left out: its plate forces, 769 and 724, do not balance.

.. figure:: images/viscosity-3d.png
    :alt: Effective viscosity against solid fraction in three dimensions
    :align: center
    :width: 600px

    Effective viscosity of the three-dimensional sweep. The curve is :math:`\mu(1+\tfrac{5}{2}\phi)`.

The three points remain on the curve, with an error between 0.25 % and 0.39 %. This confirms Einstein's formula in the dilute range. These ten centers overlap before 5 % is reached.

-------------------------
Possibility for Extension
-------------------------

Fit the closure :math:`\mu_{eff}/\mu = 1 + \tfrac{5}{2}\phi + k\phi^{2}` on these runs, and train a regressor to predict :math:`k` from the solid fraction and the minimum gap.

----------
References
----------

.. [#einstein1906] \A. Einstein, "Eine neue Bestimmung der Moleküldimensionen," *Ann. Phys.*, vol. 324, no. 2, pp. 289–306, Feb. 1906, doi: `10.1002/andp.19063240204 <https://doi.org/10.1002/andp.19063240204>`_.

.. [#einstein1911] \A. Einstein, "Berichtigung zu meiner Arbeit: Eine neue Bestimmung der Moleküldimensionen," *Ann. Phys.*, vol. 339, no. 3, pp. 591–592, 1911, doi: `10.1002/andp.19113390313 <https://doi.org/10.1002/andp.19113390313>`_.

.. [#wikiintrinsic] "Intrinsic viscosity," *Wikipedia*. [Online]. Available: `https://en.wikipedia.org/wiki/Intrinsic_viscosity <https://en.wikipedia.org/wiki/Intrinsic_viscosity>`_.https://safirh-portail.polymtl.ca/