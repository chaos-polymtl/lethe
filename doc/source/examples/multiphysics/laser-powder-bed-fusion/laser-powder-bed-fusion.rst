===================================
Static Irradiation of a Bare Plate
===================================

This example simulates the static irradiation of a Ti6Al4V bare plate. It is based on the experimental work of Cunningham *et al.* and includes the relevant phenomena involved in the laser powder bed fusion manufacturing process. 

****

--------
Features
--------

- Solver: ``lethe-fluid`` 
- Volume of fluid (VOF) and Heat Transfer (HT)
- Unsteady problem with phase change handled by an adaptive BDF2 time-stepping scheme
- Python scripts for data postprocessing

****

---------------------------
Files Used in This Example
---------------------------

All files mentioned below are located in the example's folder (``examples/multiphysics/static-irradiation``).

- Parameter file: ``static-irradiation.prm``
- Mesh file: ``mesh/bare-plate.msh``
- Postprocessing Python script to extract quatities of interest: ``static-irradiation.py``

****

-----------------------
Description of the Case
-----------------------

Laser powder bed fusion is a manufacturing process using a laser to selectively melt and consolidate, layer-by-layer, a metal powder. Simply, it corresponds to 3D printing with a metal powder as the raw material. The main laser-material interaction takes place at the melt pool scale where the flow dynamics involve multiple driving forces:

- phase change due to laser heating
- surface tension due to the small scale of the melt pool
- evaporative cooling and recoil pressure due to temperature reaching the boiling point

In this example, we are considering the static irradiation of Ti6Al4V bare plate (without powder) to study the melt pool dynamics. We simulate an irradiation of :math:`2 \;\text{ms}` by a laser beam with a diameter of :math:`140\;\mu\text{m}` and a power of :math:`156\;\text{W}`. The following figure shows the case setup, which is based on the experimental work of Cunningham *et al.*

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/new-benchmark.png                                                                             |
|     :align: center                                                                                                |
|     :width: 620                                                                                                   |
|     :name: Case setup                                                                                             |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

The dimensions (:math:`H, \Delta h`, and :math:`L`) and the Dirichlet boundary condition values (:math:`\vec{u}_{\text{in}}, T_\text{in}`, and :math:`T_\text{0}`) are listed bellow.

+---------------------------+---------------------------+----------------------------+-----------------------------+
| Parameter                 | Value                     | Parameter                  | Value                       |
+---------------------------+---------------------------+----------------------------+-----------------------------+
| :math:`H`                 | :math:`170\;\mu\text{m}`  | :math:`\vec{u}_{\text{in}}`| :math:`0.1\;\text{m s}^{-1}`|
+---------------------------+---------------------------+----------------------------+-----------------------------+
| :math:`\Delta h`          | :math:`20\;\mu\text{m}`   | :math:`T_{\text{in}}`      | :math:`298\;\text{K}`       |
+---------------------------+---------------------------+----------------------------+-----------------------------+
| :math:`L`                 | :math:`600\;\mu\text{m}`  | :math:`T_{\text{0}}`       | :math:`298\;\text{K}`       |
+---------------------------+---------------------------+----------------------------+-----------------------------+

The governing equations

--------------
Parameter File
--------------

Simulation Control
~~~~~~~~~~~~~~~~~~

The time integration is handled by a 2nd-order backward differentiation scheme (bdf2) with a maximum time-step of :math:`\Delta t = 5.0 \times 10^{-9} \; \text{s} < \Delta t_\sigma` which corresponds to the capillary time-step constraint (see :doc:`capillary wave example <../capillary-wave/capillary-wave>`). We use adaptive time stepping with a maximum CFL of :math:`0.2` to prevent instability resulting from the explicit coupling of the NS and HT solvers through the recoil pressure. 

.. code-block:: text

    subsection simulation control
      set method           = bdf2
      set time end         = 0.002
      set time step        = 5.0e-9
      set adapt            = true
      set max cfl          = 0.2
      set max time step    = 5.0e-9
      set output name      = static-irradiation
      set output path      = output/
      set output frequency = 100
    end
    
Multiphysics
~~~~~~~~~~~~

In the ``multiphysics`` subsection, we enable both the VOF and HT solvers.

.. code-block:: text

    subsection multiphysics
      set VOF           = true
      set heat transfer = true
    end
    
Mesh and box refinement
~~~~~~~~~~~~~~~~~~~~~~~

The coarse level mesh considered for this example is generated with Pointwise to enable the imposition of the inlet and outlet boundary conditions described in the figure above. It is then uniformly refined :math:`5` times and box refinement is used to insure a well discretized metal-gas interface.

.. code-block:: text

    subsection mesh
      set type               = gmsh
      set file name          = ./mesh/2d-benchmark.msh
      set initial refinement = 5
    end

    subsection box refinement
      subsection mesh
        set type               = dealii
        set grid type          = subdivided_hyper_rectangle
        set grid arguments     = 8,1 : 0,0.3925: 0.6,0.4675: false
        set initial refinement = 0
      end
      set initial refinement = 3
    end

Mesh Adaptation
~~~~~~~~~~~~~~~

As the laser heats the metal-gas interface, the melt pool deepens and the solid-liquid interface reach the bottom boundary of the box refinement. Hence, we dynamically adapt the mesh using the ``temperature`` as the refinement ``variable`` to keep a well discretized melt pool. We choose :math:`8` as the ``min refinement level`` and :math:`5` as the ``max refinement level``. The mesh is adapted each :math:`20` iterations to reduce the computational cost by setting ``frequency = 20``. Note that the ``fraction coarsening`` is set to :math:`0.0` to avoid coarsening in the center of the melt pool, where the temperature gradient is less important than at the solid-liquid interface.

.. code-block:: text

    subsection mesh adaptation
    set type                    = kelly
    set variable                = temperature
    set fraction type           = fraction
    set max refinement level    = 8
    set min refinement level    = 5
    set frequency               = 20
    set fraction refinement     = 0.4
    set fraction coarsening     = 0.0
    end
    
Boundary Conditions
~~~~~~~~~~~~~~~~~~~

In the ``boundary conditions`` subsection, we set the boundary conditions described in the figure above. Note that the ``id`` of each boundary is based on the coarse level mesh 

.. code-block:: text

    subsection boundary conditions
      set number = 6
      subsection bc 0
        set id   = 2
        set type = noslip
      end
      subsection bc 1
        set id   = 5
        set type = noslip
      end
      subsection bc 2
        set id   = 6
        set type = outlet
        set beta = 0
      end
      subsection bc 3
        set id   = 7
        set type = slip
      end
      subsection bc 4
        set id   = 4
        set type = function
        subsection u
          set Function expression = 100.0 
        end
        subsection v
          set Function expression = 0
        end
      end
      subsection bc 5
        set id   = 3
        set type = noslip
      end
    end