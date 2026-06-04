=============================
Waveguide
=============================

This example showcases an electromagnetic problem that we can solve analytically, allowing us to test the results of our simulation. It is the first example using the DPG method (Discontinuous Petrov-Gallerkin method) of the time-harmonic section of lethe.

-----------------------------
Features
-----------------------------

- Solvers: ``lethe-fluid``
- Steady-state problem
- Time-harmonic functions



----------------------------
Files Used In This Example
----------------------------

All files mentioned below are located in the example's folder (``examples/multiphysics/waveguide``)

- Base case parameter file (:math:`\mathrm{f=2.45 GHz}`): ``parameters.prm``
- Postprocessing Python script for the :math:`\mathrm{Re}=400` case: ``convergence.py``


---------------------------
Description of the Case
---------------------------

Geometry
~~~~~~~~

The waveguide at stake has a square cross section.It consists of four perfectly reflective walls and two empty sections: one for the wave’s entry and the other for its exit. The wave propagates along the z-axis. It is described on the next figure:

.. image:: images/schematic.png
    :alt: schematic
    :align: center
    :name: geometry
    :width: 500


The dimensions chosen here are a width a = 0.25m,  a height b = 0.25m and a length L = 1m.

Physical problem
~~~~~~~~~~~~~~~~

This simulation calculates the stationnary electromagnetic field in the waveguide. Even though the waveguide is empty, we aim to calculate the fields in matter when we will add the particle flow in our problem. Thus, we start with Maxwell's equations in matter, where :math:`E` and :math:`H` denote the electric and magnetic fields:

:math:`\nabla \times \mathbf(E) = \frac{\partial (\mu \mathbf(H))}{\partial t} \quad \mathrm{(Farady's law)}`

:math:`\nabla \times \mathbf(H) = \mathbf(J) - \frac{\partial (\varepsilon \mathbf(E))}{\partial t} \quad \mathrm{(Ampère-Maxwell's law)}`

Then, we use the time-harmonic ansatz : :math:`\begin{align}
\mathbf{E}(\mathbf{x}, t) &= \Re\left\{\mathbf{E}_{\text{spatial}}(\mathbf{x})\, e^{-i\omega t}\right\} \\
\mathbf{H}(\mathbf{x}, t) &= \Re\left\{\mathbf{H}_{\text{spatial}}(\mathbf{x})\, e^{-i\omega t}\right\}
\end{align}`

This form of solution introduces :math:`\omega`, the angular frequency of the electromagnetic wave, and 12 unknowns. In fact, for each of the :math:`E` and :math:`H` fields, we calculate the real and imaginary parts in each of the three spatial directions. 

Boundary conditions
~~~~~~~~~~~~~~~~~~~

.. warning::
    From now on, the equations are written in a dimensionless form.

There are three types of boundary conditions in this problem, this explains why the surfaces of the waveguide are sorted in three groups. There is first the inlet :math:`\Gamma_1`, the outlet :math:`\Gamma_2` and finally the metal walls :math:`\Gamma_3`.

- Inlet :math:`\Gamma_1`: waveguide port boundary condition at the inlet of the waveguide to excite the :math:`\mathrm{TE}_{mn}` mode at z = 0
- Outlet :math:`\Gamma_2`: impedance matching boundary condition to minimize the reflections (the waveguide is theorically infinite). We calculated the dimensionless admittance in such conditions: :math:`\frac{1}{Z_s}=\frac{k_zc}{\omega \mu_r}`
- Metal walls :math:`\Gamma_3`: perfect electric conductor (PEC) boundary conditions

This results in the three following equations:

.. math::
    \begin{align}
    \nabla \times \mathbf{H} + \frac{k_z}{\omega \mu_r} \mathbf{n} \times (\mathbf{E} \times \mathbf{n}) &=  \nabla \times \mathbf{H}_{\mathrm{TE}_{mn}} + \frac{k_z}{\omega \mu_r} \mathbf{n} \times (\mathbf{E}_{\mathrm{TE}_{mn}} \times \mathbf{n}) \quad (\text{on } \Gamma_1) \\
    \nabla \times \mathbf{H} + \frac{k_z}{\omega \mu_r} \mathbf{n} \times (\mathbf{E} \times \mathbf{n}) = 0 \quad (\text{on } \Gamma_2) \\
    \mathbf{n} \times \mathbf{E} = 0 \quad (\text{on } \Gamma_3)
    \end{align}

.. note:: 
    :math:`\mathrm{TE_mn}` refers to "Tranverse electric mode". This means that regardless of the values of m and n, the z-component of the vector E is always zero. Therefore, the pair (m,n) refers to the fact that the x-component of the electric field is excited in the m-th mode and the y-component of the electric field is excited on the n-th mode. 


Dimensionless analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



Solution of a :math:`\mathrm{TE_mn}` mode that we use for this test case is given by:

.. math::
    \mathbf{E} = i\frac{\omega\mu_r}{k_c^2} \begin{bmatrix}
    -k_y \cos(k_x x) \sin(k_y y) \\
    k_x \sin(k_x x) \cos(k_y y) \\
    0
    \end{bmatrix} e^{ik_z z}

.. math::
    \mathbf{H} =  \begin{bmatrix}
    -i\frac{k_z k_x}{k_c^2} \sin(k_x x) \cos(k_y y) \\
    -i\frac{k_z k_y}{k_c^2} \cos(k_x x) \sin(k_y y) \\
    \cos(k_x x)\cos(k_y y)
    \end{bmatrix} e^{ik_z z}
    
With the dimensionless parameters:
:math:`k_x = \frac{Lm\pi}{a} \ , \ k_y = \frac{Ln\pi}{b} \ , \ k_c^2 = L^2(k_x^2 + k_y^2) \ \text{and} \ k_z = \sqrt{\omega^2 \varepsilon_{eff,r}\mu_r - k_c^2} \ , \ \omega = \frac{L2\pi f}{c}`



--------------
Parameter File
--------------

Mesh adaptation
~~~~~~~~~~~~~~~
The first parameter we can change is the the number n of mesh adaptations:

.. code-block:: text

    subsection simulation control
        set method            = steady
        set output frequency  = 1
        set output path = ./degree_1/
        set number mesh adapt = n
    end

.. Note:: 
    A normal lab PC (12 cores, 64 Go RAM) can only deal with 3 mesh adaptations. 


The following figure represents a visualation of the results on paraview with the default settings and 4 mesh adptations.

.. image:: images/Mesh_adaptation.png
    :alt: standing-wave-mesh
    :align: center
    :name: standing-wave-mesh
    :width: 500


Frequency of the electromagnetic wave
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The frequency of the inlet is originally set up to :math:`f = \mathrm{2.45 GHz}` because it the nominal frequency of the MW in the industry. This value can still be modified at line 84 : ``set electromagnetic frequency  = 2.45e9``

.. Note:: 
    Not only has the frequency to be changed, but also :math:`k_z = \sqrt{\omega^2 \varepsilon_{eff,r}\mu_r - k_x^2 - k_y^2 },\ and \ k = \sqrt{k_x^2 +k_y^2 + k_z^2}` at line 153, and the dimensionless admittance real part :math:`\Re (G)` at line 147 :math:`\Re (G) =\frac{1}{Z_s}=\frac{k_zc}{\omega \mu_r}` because they depend on :math:`f` through :math:`\omega = \frac{2\pi Lf}{c}`. 

.. warning::
    The simulation may not run for low frequencies. Indeed, if :math:`k_z^2 = \omega^2 \varepsilon_{eff,r}\mu_r - k_x^2 - k_y^2` becomes negative, the mathematical solution still exists as evanescent waves, but it does not represent a physical system.

The different Transverse Electric mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :math:`\mathrm{TE_mn}` mode can be modified at lines 97 and 98, the defaultsettings are:

.. code-block:: text

    subsection waveguide mode
      set mode type    = TE
      set mode order m = 1
      set mode order n = 0
    end

.. Note:: 
    Not only have the modes to be changed, but also :math:`k_z = \sqrt{(L \omega)^2 \epsilon_{eff}\mu - k_x^2 - k_y^2 },\ and \ k = \sqrt{k_x^2 +k_y^2 + k_z^2}` at line 153, and the dimensionless admittance real part :math:`\Re (G)` at line 147 :math:`\Re (G) =\frac{1}{Z_s}=\frac{k_z}{\omega \mu_r}` because they depend on :math:`(m,n)` through :math:`k_x = \frac{m \pi}{a} \ , \ k_y = \frac{n \pi}{b}`. 

.. warning:: 
    If you want to try high modes, you have to adapt the frequency in order for :math:`k_z` to keep a real value. Indeed, the number :math:`k_z^2 = (L \omega)^2 \epsilon_{eff}\mu - k_x^2 - k_y^2` can be negative.

Here are paraview visuals with different modes. In order to maximise the resolution, we calculated the parameters with :math:`\lambda = L = 1m`, which is equivalent to fix the dimensionless parameter :math:`k_z = 2\pi`. This represents only one spatial period of the wave. :math:`TE_{10}`, :math:`TE_{01}` and :math:`TE_{21}`:

.. image:: images/TE10.png
    :alt: TE10
    :align: center
    :name: TE10
    :width: 500

.. image:: images/TE01.png
    :alt: TE01
    :align: center
    :name: TE01
    :width: 500

.. image:: images/TE32.png
    :alt: TE32
    :align: center
    :name: TE32
    :width: 500


Order of the FEM solver
~~~~~~~~~~~~~~~~~~~~~~~
You can also change the order of the polynomes used in the Finite Element Method on line 39 : ``set electromagnetics trial order = x``

.. note::
    If you change the trial order, you have to change the test order on line 40 : ``set electromagnetics test order  = y``. in order to have :math:`test \quad order = trial \quad order + 1`



Dielectric caracteristics of the medium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also change the physical properties of the medium on the subsection line 43. These are the default settings for air or void:

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
    
Here is a comparison of the result in a non-dissipative medium (air) and a disspative one with :math:`\Im(\mu_r)=0.06 \ , \ \Im(\varepsilon_r)=0.06`:

.. image:: images/TE10.png
    :alt: TE10V2
    :align: center
    :name: TE10V2
    :width: 500

.. image:: images/Dissipative_medium.png
    :alt: Dissipative medium
    :align: center
    :name: Dissipative medium
    :width: 500

-----------------------
Results and Discussion
-----------------------

Here is the error of the default setting simulation (:math:`\lambda = L`) depending on the mesh adaptation and the order of the simulation :

.. image:: images/model_validation.png
    :alt: final graph
    :align: center

As predicted, the error is on the order :math:`h(o^{k+1})` for a simulation on the order k.

