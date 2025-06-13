==========================
Heating of a Packed Bed
==========================

This example simulates the heating of a packed bed using the discrete element method (DEM) and heat transfer models. It is based on the validation case of Beaulieu [#Beaulieu2020]_ with stainless steel.
More information regarding the DEM parameters are given in the Lethe documentation, i.e. `DEM parameters <../../../parameters/dem/dem.html>`_.


----------------------------------
Features
----------------------------------

- Solvers: ``lethe-particles``
- Multiphysic DEM
- Three-dimensional problem
- Moving solid surface
- Heating solid surface
- Post-processing using `Python <https://www.python.org/>`_, `PyVista <https://docs.pyvista.org/>`_, `lethe_pyvista_tools <https://github.com/chaos-polymtl/lethe/tree/master/contrib/postprocessing>`_, and `ParaView <https://www.paraview.org/>`_.


----------------------------
Files Used in This Example
----------------------------

All files mentioned below are located in the example's folder (``examples/dem-mp/3d-heated-packed-bed``).

- Parameter file to load particles: ``load-packed-bed.prm``
- Parameter file for the simulation: ``heat-packed-bed.prm``
- Geometry file solid surface: ``square-top.geo``
- Post-processing Python script: ``bed-postprocessing.py``


-------------------------
Description of the Case
-------------------------

During the loading stage (:math:`0-7` s), :math:`8849` particles are inserted with a temperature of :math:`20°C` in a rectangular box to form a packed bed. Then, a moving solid surface is used to even out the height of the particles on top. After several back and forth, the top wall stops at a height of :math:`0.198` m. The height of the top wall is set so that it pushes just enough on the particles to enable heat transfer. At :math:`7` s, the temperature of the top wall is set to :math:`53°C` and the particles are heated through the top wall until :math:`4000` s, to match the experiment led by Beaulieu [#Beaulieu2020]_.


--------------
Parameter File
--------------

Mesh
~~~~

The domain we simulate is a rectangular box which is :math:`0.1\times0.1\times0.4` meters and is made using the deal.ii grid generator.

.. code-block:: text

    subsection mesh
      set type               = dealii
      set grid type          = hyper_rectangle
      set grid arguments     = 0.0 , 0.0 , 0.0 : 0.1 , 0.1 , 0.4 : false
      set initial refinement = 2
    end

Insertion Info
~~~~~~~~~~~~~~~~~

In the loading stage, particles are inserted in an insertion box in the upper part of the domain, with a temperature of :math:`20°C`. :math:`900` particles are inserted at each insertion step.

.. code-block:: text

    subsection insertion info
      set insertion method                               = volume
      set inserted number of particles at each time step = 900
      set insertion frequency                            = 10000
      set insertion box points coordinates               = 0.001, 0.001, 0.25 : 0.099, 0.099, 0.35
      set insertion distance threshold                   = 1.5
      set insertion maximum offset                       = 0.1
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
      set surface slope wall           = 0.117
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

For the loading of the particles, the model parameters are defined as:

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


Solid object
~~~~~~~~~~~~~~~

From :math:`3` to :math:`7` s, the moving solid surface is used to even out the height of the particles on top of the packed bed and it then stops at a height of :math:`0.198` m. At :math:`7` s, the temperature of the solid object is changed from :math:`20°C`  to :math:`53°C` to heat the packed bed.

.. code-block:: text

    subsection solid objects
      subsection solid surfaces
        set number of solids = 1
        subsection solid object 0
          subsection mesh
            set type               = gmsh
            set file name          = square_top.msh
            set simplex            = true
            set initial refinement = 0
          end
          subsection translational velocity
            set Function expression = 0 ; 0 ; if( (t>3 && t<7), if(t<3.359,-0.5, if(t<3.73,0.5, if (z>0.201,-0.5,if(t>5.3 && z>0.198,-0.5,0)))), 0 )
          end
          subsection angular velocity
            set Function expression = 0 ; 0 ; 0
          end
          set center of rotation    = 0 , 0 , 0.38
          set thermal boundary type = isothermal
          subsection temperature
            set Function expression = if(t>7,53,20)
          end
        end
      end
    end


Simulation Control
~~~~~~~~~~~~~~~~~~

The ``simulation control`` subsection is the main difference between the two parameter files. While the loading is :math:`7` s long with a time step of :math:`2.5\times10^{-5}`, the heating lasts until :math:`4000` s with a time step of :math:`1` s.

.. code-block:: text

    subsection simulation control
      set time step         = 2.5e-5
      set time end          = 7
      set log frequency     = 2000
      set output frequency  = 2000
      set output path       = ./output/
      set output boundaries = true
    end

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

