==========================
Heated Packed Bed
==========================

This example simulates the heating of a packed bed using the discrete element method (DEM) and heat transfer models. It is based on the validation case of Beaulieu [#Beaulieu2020]_ with stainless steel.
More information regarding the DEM parameters is given in the Lethe documentation, i.e. `DEM parameters <../../../parameters/dem/dem.html>`_.


----------------------------------
Features
----------------------------------

- Solvers: ``lethe-particles``
- Multiphysic DEM
- Three-dimensional problem
- Moving solid surfaces
- Heating solid surface
- Post-processing using `Python <https://www.python.org/>`_, `PyVista <https://docs.pyvista.org/>`_, `lethe_pyvista_tools <https://github.com/chaos-polymtl/lethe/tree/master/contrib/postprocessing>`_, and `ParaView <https://www.paraview.org/>`_.


----------------------------
Files Used in This Example
----------------------------

All files mentioned below are located in the example's folder (``examples/dem-mp/3d-heated-packed-bed``).

- Parameter file to load particles: ``load-packed-bed.prm``
- Parameter file for the simulation: ``heat-packed-bed.prm``
- Geometry files for solid surfaces: ``square-top.geo``, ``square-side.geo``, ``square-rake.geo``
- Post-processing Python script: ``bed-postprocessing.py``


-------------------------
Description of the Case
-------------------------

This example is run in three stages. 
First, during the loading stage (:math:`0-6` s), :math:`8849` particles are inserted with a temperature of :math:`20°C` in a rectangular box. In this part of the simulation, :math:`g = 9.81 \times e_x` because we insert the particles from the side of the box. When particles are still enough, we use a solid surface to even out the particles on the side of the packed bed. Then, a second solid surface is placed on the side to close the box and keep the particles in place for the next stage of the simulation.
The second stage (:math:`6-7` s) is where gravity is changed to :math:`g = -9.81 e_z` to be in the same direction as the experiment from Beaulieu [#Beaulieu2020]_ . We let the particles fall into their places for :math:`1` s. 
Finally, during the heating stage (:math:`7-4000` s), the temperature of a solid surface placed on top on the packed bed is set to :math:`53°C` and the particles are heated through this top wall.


--------------
Parameter File
--------------

Mesh
~~~~

The domain we simulate is a rectangular box with dimensions :math:`0.3\times0.1\times0.2` meters and is made using the deal.ii grid generator.

.. code-block:: text

    subsection mesh
      set type               = dealii
      set grid type          = hyper_rectangle
      set grid arguments     = -0.2 , 0.0 , 0.0 : 0.1 , 0.1 , 0.2 : false
      set initial refinement = 2
    end

Insertion Info
~~~~~~~~~~~~~~~~~

In the loading stage, particles are inserted through the side of the box, with a temperature of :math:`20°C`. :math:`3400` particles are inserted at each insertion step (except the last one).

.. code-block:: text

    subsection insertion info
      set insertion method                               = volume
      set inserted number of particles at each time step = 3400
      set insertion frequency                            = 10000
      set insertion box points coordinates               = -0.199, 0.001, 0.001 : -0.03, 0.099, 0.199
      set insertion distance threshold                   = 1.5
      set insertion maximum offset                       = 0.6
      set insertion prn seed                             = 17
      subsection initial temperature function
        set Function expression = 20
      end
    end


Lagrangian Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :math:`8849` particles are mono-dispersed, with a diameter of :math:`6.4` mm.

The physical properties of the steel particles, the walls and the interstitial gas were all chosen to match those used by Beaulieu in her experiment and simulation.

.. code-block:: text

    subsection lagrangian physical properties
      set g                        = 0.0, 0.0 , -9.81
      set number of particle types = 1
      subsection particle type 0
        set size distribution type            = uniform
        set diameter                          = 6.4e-3
        set number of particles               = 8849
        set density particles                 = 7747
        set young modulus particles           = 5e6
        set poisson ratio particles           = 0.29
        set restitution coefficient particles = 0.8
        set friction coefficient particles    = 0.7
        set rolling friction particles        = 0.02
        set real young modulus particles      = 200e9
        set thermal conductivity particles    = 42
        set specific heat particles           = 464
        set microhardness particles           = 3e9
        set surface slope particles           = 0.056
        set surface roughness particles       = 19.e-9
        set thermal accommodation particles   = 0.7
      end
      set young modulus wall           = 5e6
      set poisson ratio wall           = 0.33
      set restitution coefficient wall = 0.8
      set friction coefficient wall    = 0.7
      set rolling friction wall        = 0.02
      set real young modulus wall      = 100e9
      set thermal conductivity wall    = 250
      set microhardness wall           = 1.8e9
      set surface slope wall           = 0.056
      set surface roughness wall       = 0.1e-9
      set thermal accommodation wall   = 0.7
      set thermal conductivity gas     = 0.027
      set specific heat gas            = 1006
      set dynamic viscosity gas        = 1.85e-5
      set specific heats ratio gas     = 1
      set molecular mean free path gas = 68.e-9
    end


Model Parameters
~~~~~~~~~~~~~~~~

For the loading and the stage where gravity is flipped, the model parameters are defined as:

.. code-block:: text

    subsection model parameters
      subsection contact detection
        set contact detection method                = dynamic
        set dynamic contact search size coefficient = 0.9
        set neighborhood threshold                  = 1.3
      end
      subsection load balancing
        set load balance method = frequent
        set frequency           = 100000
      end
      set particle particle contact force method = hertz_mindlin_limit_overlap
      set rolling resistance torque method       = constant_resistance
      set particle wall contact force method     = nonlinear
      set integration method                     = velocity_verlet
      set solver type                            = dem_mp
    end

For the heating of the particles, the parameter ``disable position integration`` is set to ``true`` to freeze the position of the particles. This allows to use a higher time step for the evolution of the temperature. As particles are not moving, ``load balancing`` is no longer necessary.

.. code-block:: text

    subsection model parameters
      subsection contact detection
        set contact detection method                = dynamic
        set dynamic contact search size coefficient = 0.9
        set neighborhood threshold                  = 1.3
      end
      set particle particle contact force method = hertz_mindlin_limit_overlap
      set rolling resistance torque method       = constant_resistance
      set particle wall contact force method     = nonlinear
      set integration method                     = velocity_verlet
      set solver type                            = dem_mp
      set disable position integration           = true
    end


Solid Objects
~~~~~~~~~~~~~~~

Three solid surfaces are used in this example. The first one is the one used to heat the packed bed from :math:`7` s to :math:`4000` s, with a temperature of :math:`53°C`. The second one is used to even the particles on the side of the packed bed. The last one closes the box to maintain the particles within when the direction of gravity is changed.

.. code-block:: text

    subsection solid objects
      subsection solid surfaces
        set number of solids = 3
        subsection solid object 0
          subsection mesh
            set type               = gmsh
            set file name          = square-top.msh
            set simplex            = true
            set initial refinement = 0
          end
          subsection translational velocity
            set Function expression = 0 ; 0 ; 0
          end
          subsection angular velocity
            set Function expression = 0 ; 0 ; 0
          end
          set thermal boundary type = isothermal
          subsection temperature
            set Function expression = if(t>7,53,20)
          end
        end
        subsection solid object 1
          subsection mesh
            set type               = gmsh
            set file name          = square-rake.msh
            set simplex            = true
            set initial refinement = 0
          end
          subsection translational velocity
            set Function expression = if(z<0.19,0,if(t<3.6,-0.5,0)) ; 0 ; if(t>1.6 && z<0.19,0.1,0)
          end
          subsection angular velocity
            set Function expression = 0 ; 0 ; 0
          end
          set center of rotation    = 0 , 0 , 0
          set thermal boundary type = adiabatic
        end
        subsection solid object 2
          subsection mesh
            set type               = gmsh
            set file name          = square-side.msh
            set simplex            = true
            set initial refinement = 0
          end
          subsection translational velocity
            set Function expression = if(t>3.6 && x<-0.001,0.5,0) ; 0 ; 0
          end
          subsection angular velocity
            set Function expression = 0 ; 0 ; 0
          end
          set center of rotation    = -0.2 , 0 , 0
          set thermal boundary type = adiabatic
        end
      end
    end


Simulation Control
~~~~~~~~~~~~~~~~~~

For the loading stage:

.. code-block:: text

    subsection simulation control
      set time step         = 2.5e-5
      set time end          = 6
      set log frequency     = 2000
      set output frequency  = 2000
      set output path       = ./output/
      set output boundaries = true
    end

For the stage where gravity is changed:

.. code-block:: text

    subsection simulation control
      set time step         = 2.5e-5
      set time end          = 7
      set log frequency     = 2000
      set output frequency  = 2000
      set output path       = ./output/
      set output boundaries = true
    end

For the heating stage:

.. code-block:: text

    subsection simulation control
      set time step         = 1
      set time end          = 4000
      set log frequency     = 1
      set output frequency  = 1
      set output path       = ./output/
      set output boundaries = true
    end


-----------------------
Running the Simulation
-----------------------

This simulation is launched in two steps. First the particles are loaded with:

.. code-block:: text
  :class: copy-button

  mpirun -np 4 lethe-particles load-packed-bed.prm

Then we run the simulation to heat the particles:

.. code-block:: text
  :class: copy-button

  mpirun -np 4 lethe-particles heat-packed-bed.prm

.. note::
  In this example, the loading requires approximately 16 minutes, while simulating the temperature evolution requires 3 minutes on 4 cores.


---------------
Post-processing
---------------

A Python post-processing code ``bed-postprocessing.py`` is provided with this example. It is used to compare the temperature of the packed-bed at three different heights :math:`h_1 = 4.0` cm, :math:`h_2 = 6.0` cm and :math:`h_3 = 7.3` cm (:math:`h = 0.0` cm corresponds to the top wall), with the results obtained by Beaulieu for stainless steel.

.. figure:: images/heights.png
    :height: 400
    :align: center

The post-processing code can be run with the following line. The argument are the folder which contains the ``.prm`` file, the vtu id at which the loading ends (corresponding to :math:`7` s) and the height of the top wall when it has stopped.

.. code-block:: text
  :class: copy-button

    python3 bed-postprocessing.py  --folder ./ --start 139 --htop 0.198

.. note::

  The post-processing code may take a bit of time to run.

.. important::

    You need to ensure that ``lethe_pyvista_tools`` is working on your machine. Click `here <../../../tools/postprocessing/postprocessing_pyvista.html>`_ for details.


-------
Results
-------

The simulation can be visualised using Paraview as seen below.

.. figure:: images/heated-steel.png
    :width: 500
    :align: center

    Temperatures at the end of the simulation

The following figure compares the temperature of the packed-bed at three different heights :math:`h_1 = 4.0` cm, :math:`h_2 = 6.0` cm and :math:`h_3 = 7.3` cm, with the results obtained by Beaulieu for stainless steel.

.. figure:: images/mean-temperatures.png
    :width: 500
    :align: center

The results show good agreement with the experimental data. However, since the heat transfer is very sensitive to the overlap, pushing the top wall more or less against the particles affects the results a lot and particularly the speed at which heat transfer propagates from the top wall. So the position of the top wall could be set more precisely to get even better results.


---------
Reference
---------

.. [#Beaulieu2020] \C. Beaulieu, “Impact de la ségrégation granulaire sur le transfert de chaleur dans un lit rotatif,” (Order No. 28990310), Ph.D. thesis, Polytechnique Montréal, 2020. Available: `<https://www.proquest.com/dissertations-thèses/impact-de-la-ségrégation-granulaire-sur-le/docview/2626891455/se-2>`_\.

