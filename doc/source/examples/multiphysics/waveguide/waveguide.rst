..
    SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
    SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

=========
Waveguide
=========

This example is used to verify the implementation of the electromagnetic solver. It considers the propagation of an electromagnetic wave in a waveguide, a problem for which an analytical solution is available, thereby allowing us to assess whether the method recovers the expected order of convergence.

--------
Features
--------

- Solvers: ``lethe-fluid`` or ``lethe-fluid-matrix-free``
- Steady-state problem
- Use of discontinous Petrov-Galerkin (DPG) method

.. Note::
    Even though we are using ``lethe-fluid``, we have to enable the electromagnetism solver tin the subsection Multiphysics:

.. code-block:: text

    subsection multiphysics
        set fluid dynamics   = false
        set electromagnetics = true
    end

--------------------------
Files Used In This Example
--------------------------

All files mentioned below are located in the example's folder (``examples/multiphysics/waveguide``)

- Base case parameter file (:math:`\mathrm{f=2.45 GHz}` and :math:`\mathrm{TE}_{10}`): ``waveguide.prm``
- Folders with the results of the simulation for the degree n: ``degree_n``




-----------------------
Description of the Case
-----------------------

Geometry
~~~~~~~~

The simulated waveguide has a square cross section. It is made of four perfectly reflective lateral walls, an inlet port for the wave’s entry, and an outlet for its exit  as shown in the figure below. The wave propagates along the :math:`z`-axis.

.. image:: images/schematic.png
    :alt: schematic
    :align: center
    :name: geometry
    :width: 500


The dimensions chosen here are a width :math:`a = 0.25 \, \mathrm{m}`,  a height :math:`b = 0.25\, \mathrm{m}` and a length :math:`L = 1\, \mathrm{m}`.

Physical Problem
~~~~~~~~~~~~~~~~
.. note::
    Lethe time-harmonic solver always solve the set of equation in a dimensionless form, so that is the convention that will be used in what follows. 


This simulation computes the stationary electromagnetic field in the above waveguide with perfectly conducting metal wall and filled with void. Recall that the time-harmonic maxwell equation in their dimensionless form used in this test case are:

.. math::
    \begin{align*}
    \tag{Farady's law} \nabla \times \mathbf{E} - i\omega \mu_{\mathrm{r}}\mathbf{H} &= 0 \\
    \tag{Ampère-Maxwell's law} \nabla \times \mathbf{H} + i \omega \varepsilon_{\mathrm{r,eff}} \mathbf{E} &= \mathbf{J} 
    \end{align*}

With the parameters of the problem:

- Dielectric caracteristics of the medium : :math:`\mu_{\mathrm{r}}` the relative permeability and :math:`\varepsilon_{\mathrm{r,eff}}` the effective relative permittivity. This caracteristic is "effective" because it takes into account the electric conductivity.   
- Excitation : :math:`\mathbf{J}` the current density and :math:`\omega = \frac{L2\pi f}{c}` the angular frequency of the electromagnetic wave.
- Speed of light in void: :math:`c`

There are 12 unknowns. In fact, one for each of the dimensionless :math:`\mathbf{E}` and :math:`\mathbf{H}` fields, we compute the real and imaginary parts in each of the three spatial directions. 


Therefore we introduce the dimensionless wavenumbers used to describe the EM waves:

- Wave numbers corresponding to a standing wave caused by the walls: :math:`k_\mathrm{x} = \frac{Lm\pi}{a} \ , \ k_\mathrm{y} = \frac{Ln\pi}{b}`
- Norm of the vector :math:`(k_\mathrm{x},k_\mathrm{y})`: :math:`k_\mathrm{c}^2 = (k_\mathrm{x}^2 + k_\mathrm{y}^2)`
- Total wave number: :math:`k^2 = \omega \mu_{\mathrm{r}} \varepsilon_{\mathrm{eff,r}} = \sqrt{k_\mathrm{x}^2 +k_\mathrm{y}^2 + k_\mathrm{z}^2}`
- Wave number corresponding to a wave propagating along the z-axis: :math:`k_\mathrm{z} = \sqrt{\omega^2 \varepsilon_\mathrm{eff,r}\mu_\mathrm{r} - k_\mathrm{c}^2}`

.. note::
    Note that if :math:`k \leq k_\mathrm{c}`, the wave does not propagate because :math:`k_\mathrm{z}^2 \leq 0`. This remark is important for the following part on the frequency and the different modes.

The Different Transverse Modes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`\mathrm{TE}_{mn}` mode refers to "Transverse electric mode". This means that regardless of the values of :math:`m` and :math:`n`, the :math:`z`-component of :math:`\mathbf{E}` is always zero. Therefore, the pair :math:`(m,n)` refers to the :math:`x`-component of the electric field being excited in the :math:`m`-th mode and the :math:`y`-component of the electric field being excited on the :math:`n`-th mode. 

The :math:`\mathrm{TE}_{mn}` mode can be modified at lines 97 and 98, the default settings are:

.. code-block:: text

    subsection waveguide mode
      set mode type    = TE
      set mode order m = 1
      set mode order n = 0
    end


You can also switch to transverse magnetic mode by replacing “TE” with “TM”.

.. caution::
    If :math:`\mathrm{TE}_{mn}` is changed, :math:`k_\mathrm{z} = \sqrt{\omega^2 \varepsilon_{eff}\mu - k_\mathrm{x}^2 - k_\mathrm{y}^2 }`  at line 163, the dimensionless admittance real part :math:`\Re (G)=\frac{1}{Z_s}=\frac{k_\mathrm{z}}{\omega \mu_\mathrm{r}}` at line 147 must be changed because it depends on :math:`(m,n)` through :math:`k_\mathrm{x} = \frac{m \pi}{a}` and :math:`k_\mathrm{y} = \frac{n \pi}{b}`.

.. warning:: 
    If you want to try high modes, you have to adapt the frequency in order for :math:`k_z` to keep a real value. Indeed, the number :math:`k_z^2 = \omega^2 \epsilon_{\mathrm{eff,r}}\mu_\mathrm{r} - k_x^2 - k_y^2` can be negative.

Here are Paraview visuals of different modes. In order to maximise the resolution, we calculated the parameters with :math:`\lambda = L = 1\, \mathrm{m}`, which is equivalent to fix the dimensionless parameter :math:`k_\mathrm{z} = 2\pi`. This represents only one spatial period of the wave. Results of :math:`\mathrm{TE}_{10}`, :math:`\mathrm{TE}_{01}` and :math:`\mathrm{TE}_{21}` are presented below:

.. image:: images/TE_modes.png
    :alt: TE10
    :align: center
    :name: TE10
    :width: 500






Boundary Conditions
~~~~~~~~~~~~~~~~~~~

.. Note::
    Although we are not using the fluid solver, note that we cannot remove the ``subsection boundary conditions`` on line 107 otherwise the simulation fails to execute.

There are three types of boundary conditions in this problem, this explains why the surfaces of the waveguide are sorted in three groups. There is first the inlet :math:`\Gamma_1`, the outlet :math:`\Gamma_2` and finally the metal walls :math:`\Gamma_3`.

- Inlet :math:`\Gamma_1`: waveguide port boundary condition at the inlet of the waveguide to excite the :math:`\mathrm{TE}_{mn}` mode at :math:`z` = 0
- Outlet :math:`\Gamma_2`: impedance matching boundary condition to minimize the reflections (the waveguide is theoretically infinite). The dimensionless admittance in such conditions is:

  .. math::
    \frac{1}{Z_s}=\frac{k_\mathrm{z}}{\omega \mu_\mathrm{r}}

- Metal walls :math:`\Gamma_3`: perfect electric conductor (PEC) boundary conditions

These translate into the three following equations:

.. math::
    \begin{align*}
    \tag{on $\Gamma_1$} \nabla \times \mathbf{H} + \frac{k_\mathrm{z}}{\omega \mu_\mathrm{r}} \mathbf{n} \times (\mathbf{E} \times \mathbf{n}) &=  \nabla \times \mathbf{H}_{\mathrm{TE}_{mn}} + \frac{k_\mathrm{z}}{\omega \mu_\mathrm{r}} \mathbf{n} \times (\mathbf{E}_{\mathrm{TE}_{mn}} \times \mathbf{n})\\
    \tag{on $\Gamma_2$}\nabla \times \mathbf{H} + \frac{k_\mathrm{z}}{\omega \mu_\mathrm{r}} \mathbf{n} \times (\mathbf{E} \times \mathbf{n}) &= 0 \\
    \tag{on $\Gamma_3$}\mathbf{n} \times \mathbf{E} &= 0
    \end{align*}


    
Here are the lines concerning the boundary conditions in the parameter file:

.. code-block:: text

    subsection boundary conditions time harmonic maxwell
        set number = 3
        subsection bc 0
            set id   = 0, 1, 2, 3,
            set type = pec
        end
        subsection bc 1
            set id   = 4
            set type = waveguide port
        end
        subsection bc 2
            set id   = 5
            set type = impedance boundary
            subsection excitation x real part
                set Function expression = 0
            end
            subsection excitation x imag part
                set Function expression = 0
            end
            subsection excitation y real part
                set Function expression = 0
            end
            subsection excitation y imag part
                set Function expression = 0
            end
            subsection excitation z real part
                set Function expression = 0
            end
            subsection excitation z imag part
                set Function expression = 0
            end
            subsection surface admittance real part
                set Function expression = 0.968987646
            end
            subsection surface admittance imag part
                set Function expression = 0.
            end
        end
    end

Here is how each boundary condition is programmed in this file:

- Inlet: ``set type = waveguide port``
- Outlet: The dimensionless admittance :math:`\frac{1}{Z_\mathrm{s}} = \frac{k_\mathrm{z}}{\omega \mu_\mathrm{r}} = \frac{k_\mathrm{z} c}{2\pi f \mathrm{L} \mu_\mathrm{r}}` has to be calculated by the user before launching the simulation. Then, it is programmed here:

.. code-block:: text

    subsection surface admittance real part
        set Function expression = 0.968987646    
    
- Perfect reflective walls: ``set type = pec``


Dimensionless analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The solution for a :math:`\mathrm{TE_mn}` mode, which will be used to assess convergence for this test case, is given by:

.. math::
    \mathbf{E} = i\frac{\omega\mu_r}{k_c^2} \begin{bmatrix}
    -k_\mathrm{y} \cos(k_\mathrm{x} x) \sin(k_\mathrm{y} y) \\
    k_\mathrm{x} \sin(k_\mathrm{x} x) \cos(k_\mathrm{y} y) \\
    0
    \end{bmatrix} e^{ik_\mathrm{z} z}

.. math::
    \mathbf{H} =  \begin{bmatrix}
    -i\frac{k_\mathrm{z} k_\mathrm{x}}{k_c^2} \sin(k_\mathrm{x} x) \cos(k_\mathrm{y} y) \\
    -i\frac{k_\mathrm{z} k_\mathrm{y}}{k_c^2} \cos(k_\mathrm{x} x) \sin(k_\mathrm{y} y) \\
    \cos(k_\mathrm{x} x)\cos(k_\mathrm{y} y)
    \end{bmatrix} e^{ik_\mathrm{z} z}
    




Mesh Adaptation
~~~~~~~~~~~~~~~

We specify, in the subsection ``mesh adaptation``, the ``type``` to be ``uniform``:

.. code-block:: text

    subsection mesh adaptation
      set type = uniform
    end

The first parameter we can change is the the number ``n`` of mesh adaptations in the subsection ``simulation control``:

.. code-block:: text

    subsection simulation control
        set method            = steady
        set output frequency  = 1
        set output path       = ./degree_1/
        set number mesh adapt = n
    end



.. Note:: 
    A typical lab PC (12 cores, 64 Go RAM) can compute up to 3 mesh adaptations. 







Frequency of the Electromagnetic Wave
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The frequency of the inlet is originally set to :math:`f = \mathrm{2.45}\, \mathrm{GHz}` because it the nominal frequency of microwaves (MW) in the industry. This value can still be modified at line 84 : ``set electromagnetic frequency  = 2.45e9``

.. caution::
    If the frequency is changed, :math:`k_\mathrm{z} = \sqrt{\omega^2 \varepsilon_\mathrm{eff,r}\mu_\mathrm{r} - k_\mathrm{x}^2 - k_\mathrm{y}^2}` at line 163, the dimensionless admittance real part :math:`\Re (G) =\frac{1}{Z_s}=\frac{k_\mathrm{z}}{\omega \mu_\mathrm{r}}` at line 147 need to be changed since it depends on :math:`f` through :math:`\omega = \frac{2\pi Lf}{c}`.

.. warning::
    The simulation may not run for low frequencies. Indeed, if :math:`k_\mathrm{z}^2 = \omega^2 \varepsilon_\mathrm{eff,r}\mu_\mathrm{r} - k_\mathrm{x}^2 - k_\mathrm{y}^2` becomes negative, the mathematical solution still exists as evanescent waves, but it does not represent a physical system.



Degree of the FEM Solver
~~~~~~~~~~~~~~~~~~~~~~~~
You can also change the degree of the polynomes used in the Finite Element Method on lines 39 and 40: 

.. code-block:: text

    subsection FEM
        set electromagnetics trial order = 1
        set electromagnetics test order  = 2
    end

.. caution::
    If you change the trial degree, you have to change the test degree on in order to have :math:`test \, degree = trial \, degree + 1`.

Note that a folder is automatically created when the simulation is run. If you try a new order :math:`m`, you cam change the mame of the path on line 22: ``set output path       = ./degree_m/``

Physical Properties
~~~~~~~~~~~~~~~~~~~

You can also change the ``physical properties`` of the medium at the subsection (line 43). These are the default settings for air or void:

.. code-block:: text
    
    subsection physical properties
        set number of fluids = 1
        subsection fluid 0
            set electric conductivity model = constant
            et electric conductivity       = 0.

            set electric permittivity model     = constant
            set electric permittivity real part = 1.
            set electric permittivity imag part = 0.

            set magnetic permeability model     = constant
            set magnetic permeability real part = 1.
            set magnetic permeability imag part = 0.
        end
    end

.. Note::
    Adding an imaginary part to the permeability and the permittivity makes the medium dissipative. 
    
    Here is a comparison of the result in a non-dissipative medium (air) and a disspative one with :math:`\Im(\mu_\mathrm{r})=0.06 \ , \ \Im(\varepsilon_\mathrm{eff,r})=0.06`:

.. image:: images/Dissipation.png
    :alt: TE10V2
    :align: center
    :name: TE10V2
    :width: 500



----------------------
Running the simulation
----------------------

Call ``lethe-fluid`` by invoking:

.. code-block:: text

    mpirun -np 10 lethe-fluid waveguide.prm

to run the simulation using ten CPU cores.

.. warning:: 
    Make sure to compile lethe in `Release` mode and 
    run in parallel using mpirun. With three mesh adaptations, this simulation takes
    :math:`\sim \,6` minutes on :math:`10` processes. 

----------------------
Results and Discussion
----------------------

The following figure shows the :math:`\Re(\mathbf{H}_\mathrm{z})` field of 4 increasing mesh adaptations with he default settings. As expected, with decreasing element size we see better results.

.. image:: images/Mesh_adaptation.png
    :alt: standing-wave-mesh
    :align: center
    :name: standing-wave-mesh
    :width: 500



Here is the error of the default setting simulation in function of the element size (:math:`h`) and the polynomial degree (:math:`p`):

.. image:: images/model_validation.png
    :alt: final graph
    :align: center

As expected, the error is :math:`\mathcal{O}(h^{p+1})` in the :math:`\|L^2\|`  norm for the interior trial space, which validate the implementation of the solver. 

