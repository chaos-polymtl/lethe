================================
Water Injection in a Closed Cell
================================

This example simulates the compression of air in a closed cell by injection of water.
The problem is inspired by the test case of Caltagirone *et al.* [#caltagirone2011]_


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
Both files mentioned below are located in the example's folder (``examples/multiphysics/water-injection-in-a-closed-cell``).

- Parameter file: ``water-injection-in-a-closed-cell.prm``
- Postprocessing python script: ``water-injection-in-a-closed-cell-postprocessing.py``


-----------------------
Description of the Case
-----------------------

A square-shaped cell with sides of length :math:`L=0.1` is initially filled with air.
At :math:`t=0 \, \text{s}`, water starts entering the cell through the bottom side at a constant rate of :math:`\mathbf{v}=[0, 0.1]` causing the air in the cell to compress.
The initial configuration of this example is illustrated below.

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/water-injection-initial-configuration.svg                                                     |
|     :align: center                                                                                                |
|     :width: 700                                                                                                   |
|     :name: Water injection initial configuration                                                                  |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

.. note::
  In this example, gravity force is not considered.

--------------
Parameter File
--------------

Simulation Control
~~~~~~~~~~~~~~~~~~

The time integration is handled by a 2nd-order backward differentiation scheme (``bdf2``) with a variable time step.
The initial time step is set to :math:`0.005 \, \text{s}` and the simulation ends at :math:`t_{end} = 0.49 \, \text{s}`.

.. code-block:: text

    subsection simulation control
      set method           = bdf2
      set time end         = 0.49
      set time step        = 0.005
      set adapt            = true
      set max cfl          = 0.75
      set output name      = water-injection-in-a-closed-cell
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
We refer the reader to :doc:`../../../../theory/multiphase/cfd/vof` theory guide for further explanation on the ``phase filtration``.

.. code-block:: text

    subsection VOF
      set compressible = true
      subsection interface sharpening
        set enable              = true
        set threshold           = 0.5
        set interface sharpness = 1.8
        set frequency           = 25
      end
      subsection phase filtration
        set type      = tanh
        set beta      = 10
      end
    end

Initial Conditions
~~~~~~~~~~~~~~~~~~

In the ``initial conditions`` subsection, we define a cell filled with air (:math:`\phi=0`) at rest.

.. code-block:: text

    subsection initial conditions
      set type = nodal
      subsection uvwp
        set Function expression = 0; 0; 0
      end
      subsection VOF
        set Function expression = 0
      end
    end

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

At the bottom of the domain, water which is associated with the phase fraction :math:`\phi=1` is injected.
This is done in the simulation by setting the velocity of the fluid at the bottom boundary (``id = 2``) in the ``boundary conditions`` subsection and by imposing a ``dirichlet`` condition on the bottom boundary in the ``boundary conditions VOF`` subsection as shown below.

Boundary Conditions - Fluid Dynamics
************************************

.. code-block:: text

    subsection boundary conditions
      set number = 4
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
        set type = function
        subsection v
          set Function expression = 0.1
        end
      end
      subsection bc 3
        set id   = 3
        set type = noslip
      end
    end

Boundary Conditions - VOF
************************************

.. code-block:: text

    subsection boundary conditions VOF
      set number = 1
      subsection bc 0
        set id   = 2
        set type = dirichlet
        subsection dirichlet
          set Function expression = 1
        end
      end
    end

Physical Properties
~~~~~~~~~~~~~~~~~~~~

In the ``physical properties`` subsection, we define the properties of the fluids. For air, represented by ``fluid 0``, the ``isothermal_ideal_gas`` density model is used to account for the fluid's compressibility.
We refer the reader to the `Physical Properties - Density Models <https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/physical_properties.html#density-models>`_ documentation for further explanation on the isothermal compressible density model.
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

In the ``mesh`` subsection, we define a hyper cube with appropriate dimensions. The mesh is initially refined :math:`6` times to ensure adequate definition of the interface.

.. code-block:: text

  subsection mesh
    set type               = dealii
    set grid type          = hyper_cube
    set grid arguments     = -0.05 : 0.05 : true
    set initial refinement = 6
  end

Mesh Adaptation
~~~~~~~~~~~~~~~

In the ``mesh adaptation`` subsection, adaptive mesh refinement is defined for the ``phase``. ``min refinement level`` and ``max refinement level`` are set to :math:`6` and :math:`8`, respectively.

.. code-block:: text

    subsection mesh adaptation
      set type                     = kelly
      set variable                 = phase
      set fraction type            = fraction
      set max refinement level     = 8
      set min refinement level     = 6
      set frequency                = 1
      set fraction refinement      = 0.99
      set fraction coarsening      = 0.01
      set initial refinement steps = 5
    end

-----------------------
Running the Simulation
-----------------------

We can call ``lethe-fluid`` by invoking the following command:

.. code-block:: text
  :class: copy-button

  mpirun -np 8 lethe-fluid water-injection-in-a-closed-cell.prm

to run the simulation using eight CPU cores.

.. warning:: 
    Make sure to compile lethe in `Release` mode and run in parallel using mpirun. This simulation takes approximately half a minute on 8 processes.


-------
Results
-------

We compare the density (:math:`\rho_{\text{air}}`) and pressure (:math:`p_{\text{air}}`) in the air with their analytical values. The density is given by:

.. math::

  \rho_{\text{air}}=\frac{\rho_{\text{air,}\;\! \text{initial}}}{1-\frac{||\mathbf{v}||t}{H_{\text{air,}\;\! \text{initial}}}}

where :math:`\rho_{\text{air,}\;\! \text{initial}}=1.18` is the initial density of air, :math:`t` is the time and :math:`H_{\text{air,}\;\! \text{initial}}=L` is the initial height of the air volume.

From the ideal gas law, we obtain the following expression for the pressure:

.. math::

  p_{\text{air}} = (\rho_{\text{air}}-\rho_{\text{air,}\;\! \text{initial}}) \cdot R \cdot T

where :math:`R=287.05` is the specific gas constant of air and :math:`T=298.15` is the temperature of the fluid in Kelvin.

The results can be post-processed by invoking the following command from the folder of the example:

.. code-block:: text
  :class: copy-button

  python3 water-injection-in-a-closed-cell-postprocessing.py . water-injection-in-a-closed-cell.prm

.. important::

    You need to ensure that the ``lethe_pyvista_tools`` is working on your machine. Click `here <../../../tools/postprocessing/postprocessing.html>`_ for details.

The following figures present the comparison between the analytical results and the simulation results for the density and pressure evolutions evaluated at the center of the cavity in the air. A great agreement between the simulation and analytical results is observed.

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/figure-water-injection-in-a-closed-cell-density.svg                                           |
|     :align: center                                                                                                |
|     :width: 800                                                                                                   |
|     :name: Air density evolution                                                                                  |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

|

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/figure-water-injection-in-a-closed-cell-pressure.svg                                          |
|     :align: center                                                                                                |
|     :width: 800                                                                                                   |
|     :name: Air pressure evolution                                                                                 |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+


----------
References
----------

.. [#caltagirone2011] \J.-P. Caltagirone, S. Vincent, and C. Caruyer, “A multiphase compressible model for the simulation of multiphase flows,” *Comput. Fluids*, vol. 50, no. 1, pp. 24–34, Nov. 2011, doi: `10.1016/j.compfluid.2011.06.011 <https://doi.org/10.1016/j.compfluid.2011.06.011>`_\.
