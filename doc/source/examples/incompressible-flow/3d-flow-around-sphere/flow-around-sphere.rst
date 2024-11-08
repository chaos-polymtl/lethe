
==================================
Flow around a Sphere
==================================

In this example, a fluid flows past a sphere.


----------------------------------
Features
----------------------------------

- Solvers: ``lethe-fluid`` (with Q1-Q1) or  ``lethe-fluid-block`` (with Q2-Q1)
- Steady-state problem
- Displays the importance of adaptive mesh refinement
- Displays the effect of the Reynolds number on the convergence


----------------------------
Files Used in This Example
----------------------------

All files mentioned below are located in the example's folder (``examples/incompressible-flow/3d-flow-around-sphere``).

- Parameter file for :math:`\mathrm{Re}=0.1`: ``sphere-0.1.prm``
- Parameter file for :math:`\mathrm{Re}=150`: ``sphere-150.prm``
- Parameter file for :math:`\mathrm{Re}=150` using adaptive mesh refinement: ``sphere-adapt.prm``


-----------------------
Description of the Case
-----------------------

The following schematic describes the simulation and its boundary condition indices.

.. image:: images/example-4-setup.png
    :alt: The geometry and boundary conditions
    :align: center
    :name: geometry

Note that the sphere surface has a boundary index of ``0``, the inlet ``1``, the walls parallel to the flow direction have a boundary index of ``2`` and the outlet an index of ``3``. 

As this examples allows for two different Reynolds numbers as well as an adaptive mesh refinement variation, the parameter file section will specify the differences when applicable. 


--------------
Parameter File
--------------

We first establish the parameters.

Mesh
~~~~

The structured mesh is built using `gmsh <https://gmsh.info/#Download>`_. Geometry parameters can be adapted in the "Variables" section of the .geo file, as shown below. 

.. code-block:: text

  // Variables
  gr=4;             // Global refinement
  around=10;        // Refinement around the sphere
  trail=8;          // Refinement of trail of sphere
  near_sphere=1.20; // Progression of cell refinement to sphere surface
  downstream=16;    // Length of domain downstream of sphere
  upstream=-6;      // Length of domain upstream of sphere (must be negative)
  cross_section=10; // Half the cross-sectional width/height
  radius=0.5;       // Radius of sphere

.. note::

  The domain is constructed so that the wake has sufficient distance to develop downstream and there is a sufficient cross-sectional area to avoid any effects from the wall boundary conditions in the flow.

In the parameter file, the mesh is read by specifying its type and file name:

.. code-block:: text

    subsection mesh
      set type      = gmsh
      set file name = sphere.msh
    end

.. note::
  
  The mesh is no further refined for the first static mesh case (:math:`\mathrm{Re}=0.1`) and it has a total of 5823 cells. For the second static mesh case (:math:`\mathrm{Re}=150`), the ``initial refinement`` parameter is used and set to 1, so that it has a number of cells of the same magnitude as the mesh adaptation example. This mesh contains 46,584 cells.

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

The ``boundary conditions`` subsection establishes the constraints on different parts of the domain. They are the same for all three parameter sets.

.. code-block:: text

    subsection boundary conditions
      set number = 4
      subsection bc 0
        set type = noslip
      end
      subsection bc 1
        set type = function
        subsection u
          set Function expression = 1
        end
        subsection v
          set Function expression = 0
        end
        subsection w
          set Function expression = 0
        end
      end
      subsection bc 2
        set type = slip
      end
      subsection bc 3
        set type = outlet
      end
    end

There are four boundary conditions, as shown in the figure above. A ``noslip`` condition is applied on the surface of the sphere, where the velocity should be 0. The inlet velocity is set to `u=1m/s`, and the boundaries of the domain that are parallel to the flow direction have a ``slip`` boundary condition. The outlet boundary condition also needs to be specified.


Physical Properties
~~~~~~~~~~~~~~~~~~~

This is where the parameters differ from the first and the two last examples.

Recall the definition of the Reynolds number for a flow past a sphere:

.. math::
 \mathrm{Re} = \frac{u_{\infty} D}{\nu}

Taking this into account and the fact that the sphere diameter :math:`D` is 1 m, and :math:`u_{\infty}` is set to 1 m/s, we set the ``kinematic viscosity`` parameter to obtain the targeted Reynolds number.

* In ``/examples/incompressible-flow/3d_flow_around_sphere/sphere-0.1.prm`` (Re = 0.1)

.. code-block:: text

    subsection physical properties
      subsection fluid 0
        set kinematic viscosity = 10
      end
    end

* In ``/examples/incompressible-flow/3d_flow_around_sphere/sphere-150.prm`` and ``/examples/incompressible-flow/3d_flow_around_sphere/sphere-adapt.prm`` (:math:`\mathrm{Re}=150`)

.. code-block:: text

    subsection physical properties
      subsection fluid 0
        set kinematic viscosity = 0.006666667
      end
    end

By default, simulations only contain a single fluid which is labeled ``0``.


FEM Interpolation
~~~~~~~~~~~~~~~~~

The default FEM parameters for this example use first order polynomials. They can be easily changed to Q2-Q1 elements.

We specify the interpolation order for both pressure and velocity using the ``FEM`` subsection:

.. code-block:: text

    subsection FEM
      set velocity order = 1
      set pressure order = 1
    end

.. warning:: 

    If you choose to use the ``lethe-fluid-block`` solver; only Q2-Q1 elements are supported. 

Simulation Control
~~~~~~~~~~~~~~~~~~

The parameters also slightly differ from one problem to another the cases have different flow regimes.

In fact, for the :math:`\mathrm{Re} = 0.1`, we have a laminar to creeping flow, meaning that the flow is very steady, and therefore can be solved in steady-state.

.. code-block:: text

    subsection simulation control
      set method            = steady
      set number mesh adapt = 0
      set output path       = sphere-0.1/
      set output frequency  = 1
    end

At :math:`\mathrm{Re} = 150`, the flow has separated, resulting in an unstable wake and recirculation. It is hence more difficult to converge to a steady-state solution. Therefore we use different parameters:

.. code-block:: text

    subsection simulation control
      set method                       = steady_bdf
      set time step                    = 0.1
      set adapt                        = true
      set max cfl                      = 1000
      set stop tolerance               = 1e-5
      set adaptive time step scaling   = 1.2
      set output path                  = sphere-150/
      set output frequency             = 1
    end

The ``steady_bdf`` method solves for a steady-state simulation using adjoint time stepping with a bdf1 scheme. An initial time step is used to complete a transient iteration, and with each iteration, the time step is increased. The simulation is considered to have reached steady-state when the L2 norm of the initial residual is lower than stop tolerance at the start of a non-linear solution step, i.e. until the time step is large enough that a pseudo-steady-state has been reached.

Initial Conditions
~~~~~~~~~~~~~~~~~~

An initial condition is used for the high Reynolds number problem. It is what we call a viscous initial condition that first finds a solution for the problem at :math:`\mathrm{Re} = 10` and then uses it as initial condition to find the solution at :math:`\mathrm{Re}=150`.

.. code-block:: text

    subsection initial conditions
      set type                = viscous
      set kinematic viscosity = 0.1
    end


Mesh Adaptation
~~~~~~~~~~~~~~~

To increase the accuracy of the drag coefficient, the mesh must be refined in areas of interest, such as on the front face of the sphere and in the developing wake. Therefore, adaptive mesh refinement is used to refine the mesh in such regions.

.. code-block:: text

    subsection mesh adaptation
      set type                 = kelly
      set fraction coarsening  = 0.05
      set fraction refinement  = 0.1
      set fraction type        = number
      set max number elements  = 100000
      set min refinement level = 0
      set max refinement level = 4
      set variable             = pressure
      set frequency            = 5
    end

The mesh is dynamically adapted based on an estimate of the error of the solution for the pressure (the Kelly error estimator). The refinement is based on the number of elements. This means that the number of cells refined/coarsened per iteration is based on the fraction of the number of cells, rather than the fraction of the error (where all cells which have the fraction of the error are refined/coarsened).

The ``min refinement level`` refers to the base mesh which has been used in the previous static simulations. The mesh can only become finer than this, not coarser. The ``max refinement level`` is set at 3, giving a maximum possible number of cells of 3 million. However, the ``max number elements`` limits the number of cells to 50,000 to keep the simulation within feasible computational expense.


----------------------
Running the Simulation
----------------------

Launching the simulation is as simple as specifying the executable name and the parameter file. Assuming that the ``lethe-fluid`` executable is within your path, the simulation can be launched by typing:

.. code-block:: text
  :class: copy-button

  lethe-fluid sphere-0.1.prm

or 

.. code-block:: text
  :class: copy-button

  lethe-fluid sphere-150.prm

or

.. code-block:: text
  :class: copy-button

  lethe-fluid sphere-adapt.prm

Lethe will generate a number of files. The most important one bears the extension ``.pvd``. It can be read by popular visualization programs such as `Paraview <https://www.paraview.org/>`_. 


-----------------------
Results and Discussion
-----------------------

First Case Results (:math:`\mathrm{Re}=0.1`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using Paraview, the steady-state velocity profile and the pressure profile can be visualized by creating a *slice* along the xy-plane (z-normal) that cuts in the middle of the sphere (See `documentation <https://forgeanalytics.io/blog/creating-slices-in-paraview/>`_). 

.. image:: images/velocity-0.1.png
    :alt: velocity distribution 0.1
    :align: center

.. image:: images/pressure-0.1.png
    :alt: pressure distribution 0.1
    :align: center

We can appreciate the axisymmetrical behavior of the flow. The drag on the sphere is available in the output file ``force.00.dat`` (the other force files ``force.01.dat``, ``force.02.dat``  and ``force.03.dat`` give the forces on the other boundary conditions 1, 2 and 3, respectively). 

.. note::
  We only perform one iteration, therefore we only have one line in the force file. If several iterations are carried out by further refining the mesh several lines will be obtained. The last line of the file shows the force calculated in the last iteration. Since the flow is in the x-direction, the x-direction force ``f_x`` allows us to calculate the drag force.

.. code-block:: text

  cells      f_x           f_y          f_z        f_xv      f_yv      f_zv       f_xp       f_yp      f_zp  
  5823 98.3705224612 -0.0000000785 0.0000001119 62.270588  0.000000  0.000000   36.099934 -0.000000  0.000000

Given the flow parameters, the calculated drag coefficient is 250.50, using around 6000 cells. At Re = 0.1, an analytical solution of the drag coefficient is known: :math:`C_D = 240` (see `reference <https://kdusling.github.io/teaching/Applied-Fluids/DragCoefficient.html>`_). The deviation from the analytical solution is primarily due to the size of the domain (height of the domain compared to the size of the sphere). The coarseness of the mesh can also have an impact on the result. It would be relevant to carry out a mesh refinement analysis.

Second Case Results (:math:`\mathrm{Re}=150`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We now consider the case at a Reynolds number of 150. At this value of the Reynolds number, the flow has separated, resulting in an unstable wake and recirculation. 

The velocity and pressure are once again visualised as well as the mesh used:

.. image:: images/velocity-150.png
    :alt: velocity distribution 150
    :align: center
.. image:: images/velocity-mesh-150.png
    :alt: velocity distribution mesh 150
    :align: center

.. image:: images/pressure-150.png
    :alt: pressure distribution 150
    :align: center

The drag coefficient at Re = 150 using this example simulation is 0.798, against a predicted coefficient of 0.889 (see `reference <https://kdusling.github.io/teaching/Applied-Fluids/DragCoefficient.html>`_).

Third Case Results (:math:`\mathrm{Re}=150` With an Adaptive Mesh Refinement)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using mesh adaptive refinement, the final mesh contains slightly more than 101,000 cells. The resulting velocity profile is shown without and with the underlying mesh. Refinement around the sphere and wake can be observed:

.. image:: images/adapt-without-mesh.png
    :alt: velocity distribution 150
    :align: center

.. image:: images/adapt-with-mesh.png
    :alt: pressure distribution 150
    :align: center

It is possible to observe that this mesh allows us to obtain a better velocity profile than in the previous example. The resulting drag coefficient of 0.880 is more accurate than the one determined using the static mesh, and does not take much more time to execute than the previous example. 

---------------------------
Possibilities for Extension
---------------------------

- **High-order methods:** Lethe supports higher order interpolation. This can yield much better results with an equal number of degrees of freedom than traditional second-order (Q1-Q1) methods, especially at higher Reynolds numbers. 

- **Mesh size** It would be interesting to increase the height-sphere diameter ratio and see if the drag coefficient obtained is closer to the analytical one for Re = 0.1 A mesh refinement analysis could also be carried out.

- **Dynamic mesh adaptation:** To increase accuracy further, the ``max number elements`` and ``max refinement level`` parameters of the mesh adaption can be increased.
