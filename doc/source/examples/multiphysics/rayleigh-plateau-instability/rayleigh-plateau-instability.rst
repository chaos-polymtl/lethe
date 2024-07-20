================================
Rayleigh-Plateau Instability
================================

This example simulates the transition of a continuous jet to a droplet regime under the influence of a perturbation. The case simulated in this example corresponds to the case J1 in absence of gravity from the work of Denner *et al.* [#denner2022]_ with the Weber number :math:`We = 20` and the Ohnesorge number :math:`Oh = 0.1`.

****

--------
Features
--------

- Solver: ``lethe-fluid`` 
- Volume of fluid (VOF)
- Unsteady problem handled by an adaptive BDF2 time-stepping scheme
- Bash scripts to write, launch, and postprocess multiple cases
- Python scripts for postprocessing data

****

---------------------------
Files Used in This Example
---------------------------

All files mentioned below are located in the example's folder (``examples/multiphysics/rayleigh-plateau-instability``).

- Case generation and simulation launching Bash script: ``rayleigh-plateau-launch.sh``
- Parameter file for case generation: ``rayleigh-plateau-J1-We020-Oh0_10.tpl``
- Parameter file the for 2D case with an excitation amplitude (:math:`\delta_0`) of :math:`0.2`: ``rayleigh-plateau-J1-We020-Oh0_10_delta0_20/rayleigh-plateau-J1-We020-Oh0_10_delta0_20.prm``
- Parameter file for the 3D case with :math:`\delta_0 = 0.3`: ``3D-delta0_30/rayleigh-plateau-J1-3D.prm``
- Postprocessing Python script for breakup lengths extraction: ``rayleigh-plateau-postprocess.py``
- Postprocessing Python script for code to code comparison: ``rayleigh-plateau-compare.py``
- Postprocessing Bash script: ``rayleigh-plateau-postprocess.sh``

****

-----------------------
Description of the Case
-----------------------

Surface tension is renowned for its stabilizing effects, yet it also serves a disruptive role in various applications, such as inkjet printing. There, the surface tension plays a pivotal role in breaking up the continuous inkjet into droplets, governed by the destabilizing mechanism on the interface between the air and the ink known as the Rayleigh-Plateau instability.

In this example, the Rayleigh-Plateau instability is simulated through a continuous glycerol jet that undergoes a perturbation of different excitation amplitudes. The velocity imposed at the jet inlet takes the following form:

.. math::
  u_\mathrm{inlet} = U \left[1+\delta_0 \sin{\left(2 \mathrm{\pi} f t \right)}\right] = U \left[1+\delta_0 \sin{\left(\frac{\kappa U t}{R_\mathrm{inlet}}\right)}\right]

where :math:`U` is the initial uniform velocity of the jet, :math:`\delta_0` is the excitation amplitude, :math:`f = \frac{\kappa U}{2 \mathrm{\pi} R_\mathrm{inlet}}` is the excitation frequency with :math:`\kappa = \frac{2 \mathrm{\pi} R_\mathrm{inlet}}{\lambda}` the dimensionless wavenumber (:math:`\lambda` is the wavelength) and :math:`R_\mathrm{inlet}=1.145 \times 10^{-3} \; \mathrm m` the inlet radius. Lastly, :math:`t` is the time.

.. note::
  It is assumed that at the initial state the air is static :math:`\left(\mathbf{u}_\mathrm{air} = \mathbf{0} \; \mathrm{m\, s^{-1}}\right)`.

The following figure displays the initial state:


+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/initial-state-without-gravity.svg                                                             |
|     :align: center                                                                                                |
|     :width: 1000                                                                                                  |
|     :name: Initial state of the jet                                                                               |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

.. note::
  In this example, the gravity contribution is not considered :math:`\left(Fr = \frac{u}{\sqrt{g R_\mathrm{inlet}}} \rightarrow \infty \right)`.

****

--------------
Parameter File
--------------
The parameter file for the case of :math:`\delta_0 = 0.20` is shown below.

Simulation Control
~~~~~~~~~~~~~~~~~~
The time integration is handled by a 2nd-order backward differentiation scheme (bdf2) with a maximum time-step of :math:`\Delta t = 4.4 \times 10^{-5} \; \text{s} \approx \Delta t_\sigma` which corresponds to the capillary time-step constraint (see :doc:`capillary wave example <../capillary-wave/capillary-wave>`).

.. code-block:: text

    subsection simulation control
      set method           = bdf2
      set time end         = 0.08
      set time step        = 4.4e-5
      set adapt            = true
      set max cfl          = 0.75
      set max time step    = 4.4e-5
      set output name      = rayleigh-plateau
      set output frequency = 5
      set output path      = ./output_delta0_20/
    end

Multiphysics
~~~~~~~~~~~~

The ``multiphysics`` subsection is used to enable the VOF solver.

.. code-block:: text

    subsection multiphysics
      set VOF  = true
    end


Physical Properties
~~~~~~~~~~~~~~~~~~~~

In the ``physical properties`` subsection, we define the jet fluid (``fluid 1``) as presented for case J1 in Denner *et al.* [#denner2022]_ The viscosity is deduced from the imposed Ohnesorge number :math:`\left(Oh=\frac{\mu_1}{\sigma\rho_1 R_\mathrm{inlet}} \right)` value of :math:`0.1`. The ambient fluid (``fluid 0``) is defined such that the density :math:`\left(\frac{\rho_1}{\rho_0} = 10^3 \right)` and dynamic viscosity :math:`\left(\frac{\mu_1}{\mu_0} = 10^2\right)` ratios are respected. A ``fluid-fluid`` type of material interaction is also defined to specify the ``surface tension model``. In this case, it is set to ``constant`` (default value) with the ``surface tension coefficient`` (:math:`\sigma`) set to :math:`0.0674 \; \mathrm{N \, m^{-1}}`.

.. code-block:: text

    subsection physical properties
      set number of fluids = 2
      subsection fluid 0
        set density             = 1.196
        set kinematic viscosity = 2.54e-4
      end
      subsection fluid 1
        set density             = 1196
        set kinematic viscosity = 2.54e-5
      end
      set number of material interactions = 1
      subsection material interaction 0
        set type = fluid-fluid
        subsection fluid-fluid interaction
          set surface tension coefficient = 0.0674
        end
      end
    end

Mesh
~~~~

In the ``mesh`` subsection, we define a subdivided hyper rectangle with appropriate dimensions. The mesh is initially refined :math:`7` times to ensure adequate definition of the interface.

.. code-block:: text

    subsection mesh
      set type               = dealii
      set grid type          = subdivided_hyper_rectangle
      set grid arguments     = 4 , 1 : 0, -0.01145 : 0.0916, 0.01145 : true
      set initial refinement = 7
    end

Mesh Adaptation
~~~~~~~~~~~~~~~~

In the ``mesh adaptation`` subsection, we dynamically adapt the mesh using the ``phase`` as refinement ``variable``. We choose :math:`5` as the ``min refinement level`` and :math:`8` as the ``max refinement level``. We set ``initial refinement steps = 4`` to adapt the mesh to the initial value of the VOF field.

.. code-block:: text

    subsection mesh adaptation
      set type                     = kelly
      set variable                 = phase
      set fraction type            = fraction
      set max refinement level     = 8
      set min refinement level     = 5
      set fraction refinement      = 0.99
      set fraction coarsening      = 0.001
      set initial refinement steps = 4
    end

Initial Conditions
~~~~~~~~~~~~~~~~~~

In the ``initial conditions``, we define the initial condition as presented in the figure above.
The uniform jet velocity :math:`(U = 1.569 \; \mathrm{m \, s^{-1}})` corresponds to :math:`We=\frac{\rho_1 R_\mathrm{inlet} U^2}{\sigma}=20`.

.. code-block:: text

    subsection initial conditions
      set type = nodal
      subsection uvwp
        set Function constants  = U=1.569
        set Function expression = if(y^2 <= 1.3110e-6, U, 0); 0; 0
      end
      subsection VOF
        set Function expression = if(y^2 <= 1.3110e-6, 1, 0)
        subsection projection step
          set enable = true
        end
      end
    end

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

In the ``boundary conditions`` subsection, the inlet velocity perturbation is specified as described in the `description of the case`_ with :math:`\kappa = 0.7`.

.. code-block:: text

    subsection boundary conditions
      set number = 2
      subsection bc 0
        set id   = 0
        set type = function
        subsection u
          set Function constants  = U=1.569, delta=0.2, kappa=0.7, r=1.145e-3
          set Function expression = if (y^2 <= 1.3110e-6, U*(1 + delta*sin(kappa*U*t/r)), 0)
        end
      end
      subsection bc 1
        set id                 = 2
        set type               = periodic
        set periodic_id        = 3
        set periodic_direction = 1
      end
    end

Boundary Conditions VOF
~~~~~~~~~~~~~~~~~~~~~~~

Lasty, in the ``boundary conditions VOF`` subsection we ensure that ``fluid 1`` is at the inlet.

.. code-block:: text

    subsection boundary conditions VOF
      set number = 1
      subsection bc 0
        set id   = 0
        set type = dirichlet
        subsection dirichlet
          set Function expression = if(y^2 <= 1.3110e-6, 1, 0)
        end
      end
    end

****

-----------------------
Running the Simulation
-----------------------

We can call ``lethe-fluid`` for each :math:`\delta_0` value. For :math:`\delta_0 = 0.20`, this can be done by invoking the following command:

.. code-block:: text
  :class: copy-button

  mpirun -np 14 lethe-fluid rayleigh-plateau-J1-We020-Oh0_10_delta0_20.prm

to run the simulation using fourteen CPU cores. Feel free to use more CPU cores.

.. warning:: 
    Make sure to compile Lethe in `Release` mode and run in parallel using mpirun.
    This simulation takes :math:`\sim \, 40` minutes on :math:`14` processes.

.. tip::
  If you want to **generate and launch multiple cases** consecutively, a Bash script (``rayleigh-plateau-launch.sh``) is provided. Make sure that the file has executable permissions before calling it with:

  .. code-block:: text
    :class: copy-button

    ./rayleigh-plateau-launch.sh rayleigh-plateau-J1-We020-Oh0_10.tpl "{0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.5 0.6}"

  where ``"{0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.5 0.6}"`` is the sequence of :math:`\delta_0` values of the different cases.

  .. note::
    An additional ``-ne`` argument can be added at the end before running the script if you do not wish to extract all breakup lengths but only generate the comparison figure.

****

-------
Results
-------

Simulation Results
~~~~~~~~~~~~~~~~~~

The video below displays the results for the case of :math:`\delta_0 = 0.2`.

.. raw:: html

    <iframe width="720" height="428" src="https://www.youtube.com/embed/QA8DEo3-9hA?rel=0&vq=hd720" title="2D Rayleigh-Plateau Instability with an excitation amplitude of 0.20" frameborder="0" allowfullscreen></iframe>

Satellite Droplets
~~~~~~~~~~~~~~~~~~

The video below displays the apparition of satellite droplets (secondary droplets) at at higher excitation amplitudes. Here, :math:`\delta_0 = 0.3`.

.. raw:: html

    <iframe width="720" height="428" src="https://www.youtube.com/embed/gtIBY9FRyvY?rel=0&vq=hd720" title="3D Rayleigh-Plateau Instability with an excitation amplitude of 0.30" frameborder="0" allowfullscreen></iframe>

This 3D simulation was simulated using the ``3D-delta0_30/rayleigh-plateau-J1-3D.prm`` parameter file.

.. note::
  Note that in these simulations, the mass is not perfectly conserved. It can be observed that the satellite droplets are fading away. This will be worked on in future updates.

Code to Code Comparison
~~~~~~~~~~~~~~~~~~~~~~~

We compare the dimensionless breakup length :math:`\left(\frac{L_\mathrm{b}}{R_\mathrm{jet}}\right)` with the simulation results from Denner *et al.* [#denner2022]_ :math:`L_\mathrm{b}` is the breakup length defined as **the shortest distance from the nozzle (inlet) to the tip of the continuous jet**.

The results can be postprocessed using the provided Bash script (``rayleigh-plateau-postprocess.sh``). Make sure that the file has executable permissions before calling it with:

.. code-block:: text
  :class: copy-button

  ./rayleigh-plateau-postprocess.sh denner-et-al-2022-We020.csv "{0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.5 0.6}"

with ``denner-et-al-2022-We020.csv`` being the path to the reference data csv file.

.. important::
  You need to ensure that the ``lethe_pyvista_tools`` is working on your machine. Click :doc:`here <../../../tools/postprocessing/postprocessing_pyvista>` for details.

This script extracts breakup lengths of the cases while excluding the satellite droplets.
The script then calculates an average :math:`L_\mathrm{b}` which is used to evaluate the dimensionless breakup length of the jet.

.. note::
  The script ignores the first 2 breakups of the jet as they as considered not part of the periodical behavior.

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/rayleigh-plateau_comparison_figure.png                                                        |
|     :align: center                                                                                                |
|     :width: 800                                                                                                   |
|     :name: Dimensionless breakup length comparison                                                                |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

As it can be seen above, for :math:`\delta_0 \leq 0.1`, we observe no breakup. The jet stabilizes despite the perturbation. An additional case was studied at :math:`\delta_0 = 0.12` to check the increasing stabilizing tendency of the jet for lower excitation amplitude values.
We also observe that none of the other evaluation points match with the work of Denner *et al.* [#denner2022]_  However, a similar trend in values is observed for :math:`\delta_0 \in [0.2,0.5]`. At :math:`\delta_0 = 0.6`, a huge difference is observed. This is due to the way the satellite droplets are formed. As opposed to previous simulations, the satellite droplets are formed from the broken-off part of the jet, decreasing significantly :math:`L_\mathrm{b}` as displayed in the video below. This might have not been the case in the work of Denner *et al.* [#denner2022]_

.. raw:: html

    <iframe width="720" height="428" src="https://www.youtube.com/embed/p3TXpNErbdc?rel=0&vq=hd720" title="2D Rayleigh-Plateau Instability with an excitation amplitude of 0.60" frameborder="0" allowfullscreen></iframe>

****

----------
References
----------

.. [#denner2022] \F. Denner, F. Evrard, A. A. Castrejón-Pita, J. R. Castrejón-Pita, and B. van Wachem, “Reversal and Inversion of Capillary Jet Breakup at Large Excitation Amplitudes,” *Flow Turbul. Combust.*, vol. 108, no. 3, pp. 843–863, Mar. 2022, doi: `10.1007/s10494-021-00291-w <https://link.springer.com/article/10.1007/s10494-021-00291-w>`_\.
