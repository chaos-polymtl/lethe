..
    SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
    SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

Fichera Oven
============

This example verifies the implementation of the built-in `discontinuous Petrov-Galerkin (DPG) <https://pdxscholar.library.pdx.edu/mth_fac/426/>`_ residual error estimator and the associated :math:`h`-adaptive mesh refinement for the time-harmonic Maxwell solver. The domain is a "staircase" cavity built by stacking several classical `Fichera corner <https://en.wikipedia.org/wiki/Fichera_corner>`_ reentrant-corner geometries on top of one another, and is fed by a small excitation port at its top, much like the cavity of a microwave oven -- hence the name *Fichera oven*. Because the reentrant edges of the cavity generate electromagnetic field singularities, this benchmark is an ideal stress test for adaptive, as opposed to uniform, mesh refinement.

Features
--------

- Solver: ``lethe-fluid``
- Steady-state problem
- Use of the discontinuous Petrov-Galerkin (DPG) method
- Built-in ``fichera_oven`` grid generator
- Use of the built-in DPG residual error estimator to drive :math:`h`-adaptive mesh refinement
- Postprocessing of PVD/VTU time series with `PyVista <https://pyvista.org/>`_

Files Used in This Example
--------------------------

Both files below are located in the example's folder (``examples/multiphysics/fichera_oven``).

- Parameter file: ``fichera_oven.prm``
- Postprocessing Python script: ``fichera_oven.py``

Description of the Case
-----------------------

Geometry
~~~~~~~~

The domain is generated with Lethe's built-in ``fichera_oven`` grid generator. Starting from a :math:`2\times2\times3` array of unit cubes, four cubes are removed from the top two levels so that the cross-section of the domain shrinks, level by level, as one climbs the :math:`z`-axis:

- bottom level (:math:`z\in[0,1]`): the full :math:`2\times2` array of cubes is kept;
- middle level (:math:`z\in[1,2]`): one of the four cubes is removed, leaving an L-shaped cross-section;
- top level (:math:`z\in[2,3]`): three of the four cubes are removed, leaving a single cube.

This produces the eight-cube "staircase" tower shown below. When ``colorize`` is set to ``true`` in the ``grid arguments``, the small top face of the tower (at :math:`z=3`) is assigned boundary id ``1`` -- this is the excitation port -- while every other face of the cavity is assigned boundary id ``0``.

.. image:: images/geometry.png
    :alt: geometry
    :align: center
    :name: geometry
    :width: 500

The grid is generated within the ``mesh`` subsection using the following arguments:

.. code-block:: text

    subsection mesh
        set type           = lethe
        set grid type      = fichera_oven
        set grid arguments = 0, 0, 0 : 2, 2, 3 : true
    end

where the three groups of arguments are, respectively, the ``bottom_left`` corner, the ``top_right`` corner, and the ``colorize`` flag of the bounding box in which the staircase is inscribed.

Reentrant Edges and the Need for Adaptive Refinement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A cube with an octant removed -- the elementary building block of this geometry -- is commonly referred to in the finite element literature as a *Fichera corner*, and is a classical benchmark for the reentrant-corner singularities that arise in elliptic and electromagnetic boundary value problems. Stacking several such corners, as is done here, produces a domain with five reentrant edges, shown in red in the figure below: one vertical edge where the L-shaped middle level meets the missing corner, and four horizontal edges bounding the two exposed steps at :math:`z=1` and :math:`z=2`.

.. image:: images/reentrant_edges.png
    :alt: reentrant edges
    :align: center
    :name: reentrant-edges
    :width: 500

Near each of these edges, the electromagnetic field develops a corner singularity: its regularity is limited regardless of the polynomial degree used, so that uniform mesh refinement converges at a rate dictated by the (typically low) regularity of the singularity rather than by the polynomial degree of the finite element space. Resolving these singularities efficiently therefore requires concentrating the mesh refinement near the reentrant edges while leaving the rest of the cavity coarse -- precisely what an adaptive refinement strategy driven by a rigorous, method-intrinsic error estimator, such as the DPG residual error estimator, is designed to do. This makes the Fichera oven an ideal benchmark to verify both the DPG error estimator itself and the :math:`h`-adaptive refinement loop of Lethe's time-harmonic Maxwell solver.

Physical Problem
~~~~~~~~~~~~~~~~

.. note::
    As in the :doc:`waveguide example <../waveguide/waveguide>`, Lethe's time-harmonic Maxwell solver solves the dimensionless form of the equations,

    .. math::
        \begin{align*}
        \tag{Faraday's law} \nabla \times \mathbf{E} - i\omega \mu_{\mathrm{r}}\mathbf{H} &= 0 \\
        \tag{Ampere-Maxwell's law} \nabla \times \mathbf{H} + i \omega \varepsilon_{\mathrm{r,eff}} \mathbf{E} &= \mathbf{J}
        \end{align*}

    See the :doc:`DPG formulation for time-harmonic Maxwell problems <../../../theory/multiphysics/electromagnetism/dpg_time_harmonic_maxwell>` theory guide for the full derivation and the DPG variational formulation used by the solver.

The cavity is filled with vacuum (:math:`\varepsilon_\mathrm{r}=\mu_\mathrm{r}=1`) and its walls are perfectly-conducting (``pec``), except for the small top face, on which a Dirichlet ``electric field`` condition of the form :math:`E_x = \sin(\pi y)` is imposed. This profile vanishes on the edges of the port and peaks at its center, which seeds a rich, multi-mode field inside the cavity without requiring a dedicated waveguide-port condition. The excitation frequency is set to :math:`f=\mathrm{238.567258\,MHz}`, corresponding to a free-space wavelength of about :math:`1.26\,\mathrm{m}` -- comparable to the :math:`1\,\mathrm{m}` edge length of the cubes making up the staircase, so that the field structure interacts strongly with the geometric features of the cavity, including its reentrant edges.

Since the cavity has no closed-form analytical solution, this example does not rely on a manufactured or analytical solution to assess convergence, as is done, for instance, in the :doc:`waveguide example <../waveguide/waveguide>`. Instead, it relies on:

- the DPG method's built-in, rigorous a posteriori residual error estimator, which is (up to equivalence constants) also a bound on the true error :math:`\|u_h-u\|_E` in the energy-like norm induced by the DPG test space; and
- a comparison with independently published reference results for this exact benchmark problem, obtained with the ``hp3d`` code by Petrides and Demkowicz [#Petrides2021]_.

Parameter File
--------------

The parameter file follows the structure below. Only the parameters relevant to this example are detailed; refer to the linked documentation pages for the complete description of each subsection.

Simulation Control
~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection simulation control
        set method            = steady
        set output path       = ./output/
        set output frequency  = 1
        set number mesh adapt = 5
    end

The steady-state problem is solved on the initial mesh, then the mesh is adapted and the problem solved again, and so on for ``number mesh adapt = 5`` refinement cycles, for a total of six solves. Each of these solves is written to the ``./output/`` folder as a separate time step of a PVD/VTU time series, with the "time" of each step corresponding to its refinement iteration number (from ``0`` to ``5``).

.. _fichera-oven-mesh-adaptation:

Mesh Adaptation
~~~~~~~~~~~~~~~

This is the core subsection of this example:

.. code-block:: text

    subsection mesh adaptation
        set type                = adaptive
        set variable            = electromagnetic fields
        set error estimator     = dpg
        set fraction refinement = 0.2
        set fraction coarsening = 0.05
    end

- ``set type = adaptive`` requests :math:`h`-adaptive refinement, as opposed to a ``uniform`` refinement of every cell.
- ``set error estimator = dpg`` selects the DPG built-in residual error estimator, which is the only error estimator available for the ``electromagnetic fields`` variable. Unlike the more generic ``kelly`` estimator (a jump-based indicator applicable to every physics), the ``dpg`` estimator is intrinsic to the DPG variational formulation: it is computed from the norm, in the test space, of the local residual representation function of each cell, which is available at essentially no extra cost once the DPG system has been solved on that cell. This provides a mathematically rigorous, physics-aware refinement indicator, without any user-tuned parameter.
- ``set fraction refinement = 0.2`` and ``set fraction coarsening = 0.05`` respectively mark, at every adaptation cycle, the :math:`20\%` of cells with the largest error indicator for refinement and the :math:`5\%` of cells with the smallest error indicator for coarsening (using the deal.II ``refine_and_coarsen_fixed_number`` strategy, since ``fraction type`` is left at its default value of ``number``).

For more details on these parameters, refer to the :doc:`../../../parameters/cfd/mesh_adaptation_control` documentation.

FEM
~~~

.. code-block:: text

    subsection FEM
        set electromagnetics trial degree = 3
        set electromagnetics test degree  = 4
    end

The interior fields are approximated with degree-3 discontinuous elements, while the DPG test space is enriched to degree 4.

.. caution::
    As in the waveguide example, the time-harmonic Maxwell solver requires the test space degree to always be strictly greater than the trial space degree.

Physical Properties
~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection physical properties
        set number of fluids = 1
        subsection fluid 0
            set electric conductivity model = constant
            set electric conductivity       = 0.

            set electric permittivity model     = constant
            set electric permittivity real part = 1.
            set electric permittivity imag part = 0.

            set magnetic permeability model     = constant
            set magnetic permeability real part = 1.
            set magnetic permeability imag part = 0.
        end
    end

The cavity is filled with vacuum: relative permittivity and permeability of 1, with no losses (zero conductivity and zero imaginary parts).

Mesh
~~~~

See the Geometry section above for the description of the ``fichera_oven`` grid generator and its arguments:

.. code-block:: text

    subsection mesh
        set type           = lethe
        set grid type      = fichera_oven
        set grid arguments = 0, 0, 0 : 2, 2, 3 : true
    end

Multiphysics
~~~~~~~~~~~~

.. code-block:: text

    subsection multiphysics
        set fluid dynamics   = false
        set electromagnetics = true
    end

As in the waveguide example, the ``electromagnetics`` physics is enabled and the ``fluid dynamics`` solver, which is enabled by default, is turned off.

Time-Harmonic Maxwell
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection time harmonic maxwell
        set electromagnetic frequency = 2.38567258e8
    end

This is the only parameter required in this subsection since, unlike the waveguide example, no waveguide-port inlet is used here: the excitation is instead applied directly as a Dirichlet ``electric field`` boundary condition (see below).

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

.. attention::
    Although the fluid solver is not used, the ``subsection boundary conditions`` for the fluid dynamics boundaries cannot be removed from the parameter file, otherwise the simulation fails to execute:

    .. code-block:: text

        subsection boundary conditions
            set number = 1
            subsection bc 0
                set id   = 0, 1
                set type = noslip
            end
        end

The electromagnetic boundary conditions are specified in the ``boundary conditions time harmonic maxwell`` subsection:

.. code-block:: text

    subsection boundary conditions time harmonic maxwell
        set number = 2
        subsection bc 0
            set id   = 0
            set type = pec
        end
        subsection bc 1
            set id   = 1
            set type = electric field
            subsection E x real part
                set Function expression = sin(pi*y)
            end
            subsection E x imag part
                set Function expression = 0
            end
            subsection E y real part
                set Function expression = 0
            end
            subsection E y imag part
                set Function expression = 0
            end
            subsection E z real part
                set Function expression = 0
            end
            subsection E z imag part
                set Function expression = 0
            end
        end
    end

- ``bc 0`` (id ``0``) covers every wall of the staircase except the top port and imposes a ``pec`` (perfect electric conductor) condition, :math:`\mathbf{n}\times\mathbf{E}=0`.
- ``bc 1`` (id ``1``) covers the small top face and imposes a Dirichlet ``electric field`` condition. Only the real part of :math:`E_x` is non-zero, and is set to :math:`\sin(\pi y)` -- a profile that vanishes on the edges of the port (:math:`y=0` and :math:`y=1`) and peaks at its center.

Linear Solver Control
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection linear solver
        subsection electromagnetics
            set verbosity         = verbose
            set relative residual = 1e-8
            set minimum residual  = 1e-12
            set preconditioner    = none
        end
    end

For more details on the linear solver parameters, refer to the `Linear Solver <https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/linear_solver_control.html>`_ documentation.

Running the Simulation
----------------------

Call ``lethe-fluid`` by invoking:

.. code-block:: text
    :class: copy-button

    mpirun -np 8 lethe-fluid fichera_oven.prm

to run the simulation using eight CPU cores.

.. warning::
    Make sure to compile Lethe in ``Release`` mode. With five mesh adaptations, this simulation should take on the order of a few minutes on eight processes; because the DPG error estimator concentrates refinement near the reentrant edges rather than refining the entire cavity, the cell count grows much more modestly with each iteration than it would under uniform refinement.

Once the simulation is complete, run the postprocessing script from the same folder:

.. code-block:: text
    :class: copy-button

    python3 fichera_oven.py

The script reads the ``out.pvd`` time series produced in the ``output`` folder (one time step per refinement iteration) using `PyVista <https://pyvista.org/>`_, and for each iteration:

- retrieves the number of cells and, from the surface and edge topology of the mesh, the number of interior and boundary faces and edges;
- estimates the number of degrees of freedom of the interior (degree-3 DG) and skeleton (degree-4 Nedelec) DPG trial spaces from these counts;
- computes the :math:`L^2`- and :math:`L^\infty`-like norms of the DPG error indicator, ``dpg_error_norm``, output by the solver on each cell.

It then produces two convergence plots -- error versus number of cells, and error versus number of degrees of freedom compared against the reference ``hp3d`` results -- as well as, when run with the ``--validate`` flag, a ``solution-fichera-oven.dat`` file used for automated regression testing:

.. code-block:: text
    :class: copy-button

    python3 fichera_oven.py --validate -f /path/to/output/folder

.. tip::
    The ``-f`` (or ``--folder``) argument lets you point the script to an output folder other than the default ``./output``, which is convenient when postprocessing results generated in a different directory.

Results and Discussion
----------------------

The following figure compares the convergence of the DPG error norm obtained with Lethe against the reference results of Petrides and Demkowicz [#Petrides2021]_, obtained with the independent ``hp3d`` code, as a function of the number of degrees of freedom:

.. image:: images/convergence_dofs.png
    :alt: convergence with respect to the number of degrees of freedom
    :align: center
    :name: convergence-dofs
    :width: 500

The following figure shows the same DPG error norm as a function of the number of cells in the mesh, illustrating the reduction of the error achieved over the successive adaptive refinement cycles:

.. image:: images/convergence_cells.png
    :alt: convergence with respect to the number of cells
    :align: center
    :name: convergence-cells
    :width: 500

Finally, the following figure shows the refined mesh obtained after the five adaptive refinement cycles:

.. image:: images/mesh_adaptation.png
    :alt: adaptively refined mesh
    :align: center
    :name: mesh-adaptation
    :width: 500

As expected from the discussion of the geometry above, refinement concentrates near the reentrant edges of the staircase rather than spreading uniformly across the cavity, confirming that the DPG error estimator correctly identifies the regions responsible for the largest contribution to the discretization error.

Possibilities for Extension
---------------------------

- **Compare against uniform refinement:** Set ``mesh adaptation / type = uniform`` and rerun the case. Comparing the resulting convergence curve (error versus number of degrees of freedom) with the adaptive one obtained above directly illustrates the benefit of adaptive refinement for problems with corner singularities: uniform refinement is throttled by the regularity of the reentrant-edge singularities, while adaptive refinement is not.
- **Use the Kelly error estimator:** Set ``mesh adaptation / variable = electric field, magnetic field`` and ``error estimator = kelly, kelly`` to compare the mesh refinement pattern, and the resulting convergence, obtained with a generic jump-based estimator against those obtained with the method-intrinsic DPG estimator.
- **Change the excitation:** Modify the frequency or the profile of the ``electric field`` boundary condition on the port to excite different resonant modes of the cavity and observe how the refinement pattern adapts to the new field distribution.
- **Increase the polynomial degree:** Increase ``electromagnetics trial degree`` (and correspondingly ``electromagnetics test degree``) to study how the polynomial degree affects the convergence rate away from the reentrant edges, where the solution remains smooth.

References
----------

.. [#Petrides2021] \S. Petrides and L. Demkowicz, "An adaptive multigrid solver for DPG methods with applications in linear acoustics and electromagnetics," *Comput. Math. Appl.*, vol. 87, pp. 12-26, 2021, doi: `10.1016/j.camwa.2021.01.017 <https://doi.org/10.1016/j.camwa.2021.01.017>`_\.

.. [#DemkowiczGopalakrishnan2025] \L. Demkowicz and J. Gopalakrishnan, "The discontinuous Petrov-Galerkin method," *Acta Numerica*, vol. 34, pp. 293-384, 2025, doi: `10.1017/S0962492925000021 <https://www.cambridge.org/core/journals/acta-numerica/article/discontinuous-petrovgalerkin-method/71BCF32CDE92B0924051FA31E8F54DC2>`_\.