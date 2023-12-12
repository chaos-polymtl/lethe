================================
Air Bubble Compression
================================

This example simulates the compression of an air bubble by surrounding liquid.
The problem is inspired by the test case of Caltagirone *et al.* `[1] <https://doi.org/10.1016/j.compfluid.2011.06.011>`_


--------
Features
--------

- Solver: ``lethe-fluid`` (with Q1-Q1)
- Volume of fluid (VOF)
- Isothermal compressible fluid
- Unsteady problem handled by an adaptive BDF2 time-stepping scheme
- Usage of a python script for post-processing data


---------------------------
Files Used in This Example
---------------------------

Both files mentioned below are located in the example's folder (``examples/multiphysics/air-bubble-compression``).

- Parameter file: ``air-bubble-compression.prm``
- Postprocessing python script: ``air-bubble-compression-postprocessing.py``


-----------------------
Description of the Case
-----------------------

A circular air bubble of diameter :math:`D=0.06` lies at the center of a square-shaped domain with sides of length :math:`L=0.1`.
On all four sides of the domain, water penetrates with a constant velocity norm of :math:`||\mathbf{v}||=0.0025` causing the compression of the air bubble.
The initial configuration of this example is illustrated below.


+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/air-bubble-initial-configuration.svg                                                          |
|     :align: center                                                                                                |
|     :width: 800                                                                                                   |
|     :name: Air bubble initial configuration                                                                       |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

.. note::
  In this example, gravity and surface tension forces are not considered.

--------------
Parameter File
--------------

Simulation Control
~~~~~~~~~~~~~~~~~~

Time integration is handled by a 2nd-order backward differentiation scheme (``bdf2``) with a variable time step.
The initial time step is set to :math:`0.1 \, \text{s}` and the simulation lasts :math:`2.5 \, \text{s}`.

.. code-block:: text

    subsection simulation control
      set method           = bdf2
      set time end         = 2.5
      set time step        = 0.1
      set adapt            = true
      set max cfl          = 1
      set output name      = air-bubble-compression
      set output frequency = 5
      set output path      = ./output/
    end

Multiphysics
~~~~~~~~~~~~

The ``multiphysics`` subsection is used to enable the VOF solver.

.. code-block:: text

    subsection multiphysics
      set VOF  = true
    end 

VOF
~~~

In the ``VOF`` subsection, the ``compressible``, the ``interface sharpening``, and the ``phase filtration`` features are enabled.
The enabled ``compressible`` parameter allows interface compression by adding the term :math:`\phi (\nabla \cdot \mathbf{u})` to the VOF equation.
The ``interface sharpening`` method and its parameters are explained in the :doc:`../dam-break/dam-break` example.
The ``phase filtration`` filters the phase field used for the calculation of physical properties by stiffening the value of the phase fraction.
We refer the reader to :doc:`../../../../theory/multiphysics/vof` theory guide for further explanation on the ``phase filtration``.

.. code-block:: text

    subsection VOF
      set compressible = true
      subsection interface sharpening
        set enable              = true
        set threshold           = 0.5
        set interface sharpness = 2.2
        set frequency           = 8
      end
      subsection phase filtration
        set type      = tanh
        set beta      = 10
      end
    end

Initial Conditions
~~~~~~~~~~~~~~~~~~

In the ``initial conditions`` subsection, we define the initial air bubble with a radius of :math:`D/2=0.03` surrounded by water.
An initial velocity field is used to avoid discontinuities in the solution.

.. code-block:: text

    subsection initial conditions
      set type = nodal
      subsection uvwp
        set Function expression = 0.0025*-sin(2*pi*x/0.2); 0.0025*-sin(2*pi*y/0.2);0
      end
      subsection VOF
        set Function expression = if (x^2 + y^2 < 0.03^2, 0, 1)
      end
    end

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

On all four sides of the domain, water which is associated with the phase fraction :math:`\phi=1` is injected.
This is done in the simulation by setting the velocities of the fluid in the ``boundary conditions`` subsection and by selecting the correct fluid in the ``boundary conditions VOF`` subsection with a ``dirichlet`` boundary condition on the phase fraction as shown below.

Boundary Conditions - Fluid Dynamics
************************************

.. code-block:: text

    subsection boundary conditions
      set number = 4
      subsection bc 0
        set id   = 0
        set type = function
        subsection u
          set Function expression = 0.0025
        end
      end
      subsection bc 1
        set id   = 1
        set type = function
        subsection u
          set Function expression = -0.0025
        end
      end
      subsection bc 2
        set id   = 2
        set type = function
        subsection v
          set Function expression = 0.0025
        end
      end
      subsection bc 3
        set id   = 3
        set type = function
        subsection v
          set Function expression = -0.0025
        end
      end
    end

Boundary Conditions - VOF
************************************

.. code-block:: text

    subsection boundary conditions VOF
      set number = 4
      subsection bc 0
        set id   = 0
        set type = dirichlet
        subsection dirichlet
          set Function expression = 1
        end
      end
      subsection bc 1
        set id   = 1
        set type = dirichlet
        subsection dirichlet
          set Function expression = 1
        end
      end
      subsection bc 2
        set id   = 2
        set type = dirichlet
        subsection dirichlet
          set Function expression = 1
        end
      end
      subsection bc 3
        set id   = 3
        set type = dirichlet
        subsection dirichlet
          set Function expression = 1
        end
      end
    end

Physical Properties
~~~~~~~~~~~~~~~~~~~~

In the ``physical properties`` subsection, we define the properties of the fluids. For air, represented by ``fluid 0``, the ``isothermal_ideal_gas`` density model is used to account for the fluid's compressibility.
We refer the reader to the `Physical Properties - Density Models <https://lethe-cfd.github.io/lethe/documentation/parameters/cfd/physical_properties.html#density-models>`_ documentation for further explanation on the isothermal compressible density model.
The properties of air and water at :math:`25 \, \text{°C}` are used in this example.

.. code-block:: text

    subsection physical properties
      set number of fluids = 2
      subsection fluid 0
        set density model       = isothermal_ideal_gas
        subsection isothermal_ideal_gas
          set density_ref = 1.18
          set R           = 287.05
          set T           = 298.15
        end
        set kinematic viscosity = 0.0000156
      end
      subsection fluid 1
        set density             = 1000
        set kinematic viscosity = 0.000001
      end
    end

Mesh
~~~~

In the ``mesh`` subsection, we define a hyper cube with appropriate dimensions. The mesh is initially refined :math:`7` times to ensure adequate definition of the interface.

.. code-block:: text

  subsection mesh
    set type               = dealii
    set grid type          = hyper_cube
    set grid arguments     = -0.05 : 0.05 : true
    set initial refinement = 7
  end

Mesh Adaptation
~~~~~~~~~~~~~~~

In the ``mesh adaptation`` subsection, adaptive mesh refinement is defined for the ``phase``. ``min refinement level`` and ``max refinement level`` are :math:`7` and :math:`9`, respectively. Since the size of the bubble changes, we choose a rather large ``fraction refinement`` (:math:`0.99`) and moderate ``fraction coarsening`` (:math:`0.01`).

.. code-block:: text

    subsection mesh adaptation
      set type                     = kelly
      set variable                 = phase
      set fraction type            = fraction
      set max refinement level     = 9
      set min refinement level     = 7
      set frequency                = 1
      set fraction refinement      = 0.99
      set fraction coarsening      = 0.01
      set initial refinement steps = 6
    end


-----------------------
Running the Simulation
-----------------------

We can call ``lethe-fluid`` by invoking the following command:

.. code-block:: text
  :class: copy-button

  mpirun -np 8 lethe-fluid air-bubble-compression.prm

to run the simulation using eight CPU cores. Feel free to use more.

.. warning:: 
    Make sure to compile lethe in `Release` mode and run in parallel using mpirun. This simulation takes :math:`\sim` 1.5 minute on 8 processes.


-------
Results
-------

We compare the density (:math:`\rho_{\text{air}}`) and pressure (:math:`p_{\text{air}}`) in the air bubble with their analytical values. The density is given by:

.. math::

  \rho_{\text{air}}=\frac{\rho_{\text{air,}\;\! \text{initial}}}{1-\frac{4qt}{\pi D^2}}

where :math:`\rho_{\text{air,}\;\! \text{initial}}=1.18` is the initial density of air, :math:`q = 4 \cdot ||\mathbf{v}|| \cdot L = 0.001` is the volumetric flow rate, and :math:`t` is the time.

From the ideal gas law, we obtain the following expression for the pressure:

.. math::

  p_{\text{air}} = (\rho_{\text{air}}-\rho_{\text{air,}\;\! \text{initial}}) \cdot R \cdot T

where :math:`R=287.05` is the specific gas constant of air and :math:`T=298.15` is the temperature of the fluid in Kelvin.

The results can be post-processed by invoking the following command from the folder of the example:

.. code-block:: text
  :class: copy-button

  python3 air-bubble-compression-postprocessing.py . air-bubble-compression.prm

.. important::
    You need to ensure that ``lethe_pyvista_tools`` is working on your machine. Click `here <../../../tools/postprocessing/postprocessing.html>`_ for details.

The following figures present the comparison between the analytical results and the simulation results for the density and pressure evolutions evaluated at the center of the bubble. A pretty good agreement between the simulation and analytical results is observed.


+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/figure-air-bubble-compression-density.svg                                                     |
|     :align: center                                                                                                |
|     :width: 800                                                                                                   |
|     :name: Air bubble density evolution                                                                           |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

|

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/figure-air-bubble-compression-pressure.svg                                                    |
|     :align: center                                                                                                |
|     :width: 800                                                                                                   |
|     :name: Air bubble pressure evolution                                                                          |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+


----------
References
----------

`[1] <https://doi.org/10.1016/j.compfluid.2011.06.011>`_ J.-P. Caltagirone, S. Vincent, and C. Caruyer, “A multiphase compressible model for the simulation of multiphase flows,” *Comput. Fluids*, vol. 50, no. 1, pp. 24–34, Nov. 2011, doi: 10.1016/j.compfluid.2011.06.011.