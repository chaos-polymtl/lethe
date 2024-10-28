==================================
Bubble detachment in shear flow
==================================

This example simulates the detachment of a growing bubble in a liquid shear flow.


----------------------------------
Features
----------------------------------

- Solver: ``lethe-fluid`` (Q1-Q1)
- Two phase flow handled by the Cahn-Hilliard-Navier-Stokes (CHNS) approach
- Dimensionality of the length
- Setting the mobility coefficient in problems where the lenghthscale of interest is below the critical radius in CHNS simulations
- Parametric sweep on the shear rate of the fluid
- Post-processing of the detachment time and volume of the bubble, plotting the contour of the bubble when the detachment occurs


--------------------------
Files Used in This Example
--------------------------

All files mentioned below are located in the example's folder (``examples/multiphysics/bubble-detachment-shear-flow``).

- Parameter file: ``bubble-detachment-shear-flow.prm``
- Postprocessing Python script: ``multiple_folders_bubble_detachment_post_processing.py`` (using the functions of ``functions_bubble_detachment_post_processing.py``)


-------------------------
Description of the Case
-------------------------

**INSERT A LITTLE INTRODUCTION ON THE SUBJECT HERE**

The problem consists of a rectangular box of length :math:`L = 12 \ \text{mm}`, height :math:`H = 5 \ \text{mm}` and thickness :math:`l = 5 \ \text{mm}`. The liquid phase flows from left to right following a Couette profile, that is determined with a given **shear rate**, which varies between two simulations. At the bottom wall is a circular air inlet of radius :math:`R_0 = 0.5 \ \text{mm}`. Air is injected from this inlet following a Poiseuille profile. As the time passes, the bubble, which is initialized as a semi-sphere of radius :math:`R_0`, grows, until the shear force from the surrounding liquid and the buoyancy force detach it.

The computational domain with relevant boundary conditions is described in the following figure (not to scale):

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/bubble-detachment-case.svg                                                                    |
|     :alt: Computational domain definition and initial conditions of the bubble detachment in a shear flow.        |
|     :align: center                                                                                                |
|     :name: Computational domain definition and initial conditions of the bubble detachment in a shear flow.       |
|                                                                                                                   |
|     Representation of the initial conditions of the bubble detachment in a shear flow                             |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

The quantity of interest of this problem is the detachment time :math:`t_\text{det}`. It is defined as the last time where the number of closed contour of the phase order field is equal to 1. From this we derive the detachment volume :math:`V_\text{det}` which is the bubble volume at :math:`t_\text{det}`. We perform numerous simulations by changing the shear rate and compute the detachment time and volumes using the python scripts provided. Those results are then compared to the results from Mirsandi *et al.* [#mirsandi2020]_

-----------------
Parameter File
-----------------

Simulation Control
~~~~~~~~~~~~~~~~~~

Time integration is handled by a 2nd order backward differentiation scheme (`bdf2`), for a :math:`0.5 \ \text{s}` simulation time with an initial time step of :math:`0.0005 \ \text{s}`. Time-step adaptation is enabled using ``adapt=true`` and the max CFL is :math:`0.9`. ``output boundaries`` is set to ``true`` to get a ``.vtu`` file containing the indices of the boundaries of the domain. The maximum time step is computed using the capillary time step condition.

.. code-block:: text

    subsection simulation control
      set method            = bdf2
      set output name       = article-3d
      set output frequency  = 15
      set output path       = ./outputs/output/
      set time end          = 0.5
      set adapt             = true
      set max cfl           = 0.9
      set time step         = 5e-7
      set max time step     = 7.9784e-06 # Capillary timestep condition
      set output boundaries = true
    end


Multiphysics
~~~~~~~~~~~~

The ``multiphysics`` subsection is used to enable the ``cahn hilliard`` solver.
Note that the fluid dynamics are solved by default.

.. code-block:: text

    subsection multiphysics
      set cahn hilliard = true
    end
    
Dimensionality
~~~~~~~~~~~~~~

The ``dimensionality`` subsection is used to define the unit length as :math:`0.001 \text{m} = 1 \ \text{mm}`. This setting helps with the convergence of the solver.	

.. Note:: When using the dimensionality parameters, the problem and the physical properties are rescaled using the new units specified by the user. This means that physical properties can be given their value in SI units and will automatically be rescaled. The resulting fields (velocity and pressure for instance) will also be rescaled accordingly. One exception to this is the source terms and the initial conditions, which need to be specified in rescaled units.


.. code-block:: text

    subsection dimensionality
      set length = 0.001 # meter
    end
    
Mesh
~~~~

In the ``mesh`` subsection, we specify the mesh used in this example. The grid arguments are specified so the origin of the coordinate system is located at one sixth of :math:`L`, in the middle of the transversal direction and on the bottom wall. The subdivisions are such that the cells are close to a square shape.

.. code-block:: text

    subsection mesh
      set type               = dealii
      set grid type          = subdivided_hyper_rectangle
      set grid arguments     = 9,3,3 : -2,0,-2.5 : 10,5,2.5 : true
      set initial refinement = 3
    end
    
Mesh Adaptation
~~~~~~~~~~~~~~~

The ``mesh adaptation`` section controls the dynamic mesh adaptation. Here, we choose ``phase_cahn_hilliard`` as the refinement ``variable``. The maximum and minimum refinement levels are respectively set to :math:`6` and :math:`3` with the number of ``initial refinement steps`` set to :math:`4`. This ensures the physics close to the interface to be well resolved, while keeping a coearse cell size far from the interface.

.. code-block:: text

    subsection mesh adaptation
      set type                     = kelly
      set variable                 = phase_cahn_hilliard
      set fraction type            = fraction
      set max refinement level     = 6
      set min refinement level     = 3
      set frequency                = 3
      set fraction refinement      = 0.995
      set fraction coarsening      = 0.005
      set initial refinement steps = 4
    end
    
Physical Properties
~~~~~~~~~~~~~~~~~~~

The ``physical properties`` subsection defines the physical properties of the fluids. In this example, we need first to define the properties of the surrounding liquid as that of water, hence the choice of :math:`\rho_0 = 2000 \ \text{kg}\cdot\text{m}^{-3}` and :math:`\nu_0 = 10^{-4} \ \text{m}^2\cdot\text{s}^{-1}`. The gas in the bubble is air, whose physical properties are :math:`\rho_1 = 1.23 \ \text{kg}\cdot\text{m}^{-3}` and :math:`\nu_1 = 1.455^{-5} \ \text{m}^2\cdot\text{s}^{-1}` . Since we have a water-air interface, the surface tension coefficient is : :math:`\sigma = 0.073 \ \text{N}\cdot\text{m}^{-1}`. 

In this problem, the radius of the bubble is below the critical radius (see Yue *et al.* [#yue2007]_), which means the air bubble will diffuse in the liquid phase over the timescale of the problem is the mobility coefficient is not set adequately. Yue *et al.* derive a criterion for setting the mobility coefficient :math:`D` that depends on the parameters of the problem. It is given below:

.. math::
    D = \frac{(S + S_a)R_0\varepsilon}{\sigma}
    
where :math:`S` is the shear rate related to the liquid flow, :math:`S_a` is the shear rate related to the air flow, :math:`\varepsilon` is the interface thickness and :math:`\sigma` is the surface tension coefficient. :math:`S_a` is the only unknown, it is estimated as follows:

.. math::
    S_a = \frac{v_\text{in,a}}{2R_0}
    


.. code-block:: text

    subsection physical properties
      set number of fluids = 2
      subsection fluid 0 # Water (phase = 1)
        set kinematic viscosity = 1.0016e-06
        set density             = 1000
      end
      subsection fluid 1 # Air (phase = -1)
        set kinematic viscosity = 14.55e-6
        set density             = 1.23
      end
      set number of material interactions = 1
      subsection material interaction 0
        subsection fluid-fluid interaction
          set surface tension coefficient     = 0.073
          set cahn hilliard mobility model    = constant
          set cahn hilliard mobility constant = 2.8177e-08 # Non-diffusion on problem time-scale condition
        end
      end
    end
    
Cahn-Hilliard
~~~~~~~~~~~~~

In the ``cahn hilliard`` subsection, we set the ``potential smoothing coefficient`` (soon to be deprecated) to :math:`0`. The interface thickness is set to be determined automatically based on the mesh size in the ``epsilon`` subsection. We also output the interface thickness for each time-step by setting the ``verbosity`` to ``verbose`` to know its exact value for the initial conditions.

.. code-block:: text

    subsection cahn hilliard
      set potential smoothing coefficient = 0
      subsection epsilon
        set method    = automatic
        set verbosity = verbose
      end
    end
    
Initial Conditions
~~~~~~~~~~~~~~~~~~

In the ``initial conditions`` subsection, we initialize both the phase field in the ``cahn hilliard`` subsection and the fluid velocity in the ``uvwp`` subsection. 

First the whole domain velocity is initialized to that of a Couette flow of a given shear rate. The velocity profile of a Couette flow and the associated shear rate :math:`S` are related as:

.. math::
    \mathbf{u}(y) = S\cdot y\mathbf{e}_x
    

The chemical potential field is set to :math:`0` uniformly. The interface is initialized with the equilibrium interface thickness, which requires to know the value of :math:`\epsilon` that is outputted at every iteration. Here, :math:`\epsilon = 0.0418546`.

.. code-block:: text

    subsection initial conditions
      subsection uvwp
        set Function expression = 1000*2.25*(y/5);0;0;0
      end

      subsection cahn hilliard
        set Function expression = -tanh((5e-1 - sqrt(y*y + x*x + z*z))/(1.41*0.0418546));0
      end
    end
    
Boundary Conditions
~~~~~~~~~~~~~~~~~~~

We need to set boundary conditions both for the fluid dynamics solver and the Cahn-Hilliard solver. For the latter, we constraint the angle of contact between the left side of the plate and the fluid using the ``angle_of_contact`` boundary condition of the Cahn-Hilliard solver.

.. code-block:: text

    subsection boundary conditions cahn hilliard
    set number = 1
    	subsection bc 0
    		set id          = 2 
    		set type        = angle_of_contact
    		set angle value = 50
    	end
    end
    
Then, a ``slip`` boundary condition is applied everywhere, except for the upper boundary, where it is set as ``none``.

.. code-block:: text

    subsection boundary conditions
      set number = 4
      subsection bc 0
        set id   = 2 # inner wall
        set type = slip 
      end
      subsection bc 1
        set id   = 5 # remaining walls
        set type = slip
      end
      subsection bc 2
        set id = 4 # upper surface
        set type = none
      end
      subsection bc 3
        set id   = 3 # middle
        set type = slip
      end
    end
    
Source Term
~~~~~~~~~~~

In the ``source term`` subsection, we define the gravitational acceleration. Since the unit length is the millimeter, the usual value of :math:`g` needs to be multiplied by :math:`1000`.

.. code-block:: text

    subsection source term
      subsection fluid dynamics
        set Function expression = 0; 0; -9810; 0
      end
    end
    
-----------------------
Running the Simulation
-----------------------

We call ``lethe-fluid`` by invoking:

.. code-block:: text
  :class: copy-button
  
   mpirun -np 10 lethe-fluid jurins-law-2d.prm
   
to run the simulation using ten CPU cores. Feel free to use more CPU cores.

.. warning::
    Make sure to compile Lethe in `Release` mode and run in parallel using ``mpirun``. The simulation should take 3-4 minutes for 10 processors.

-----------------
Results
-----------------

The detachment time and volume are computed and compared to the results of Mirsandi *et al.* in the following figure, which shows an excellent agreement.

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/bubble-detachment_volume.png                                                                  |
|     :alt: Plot of the detachment volume (our results and literature results) with respect to shear rate.          |
|     :align: center                                                                                                |
|     :name: Detachment volumes                                                                                     |
|                                                                                                                   |
|     Plot of the detachment volume (our results and literature results) with respect to shear rate.                |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/bubble-detachment_time.png                                                                    |
|     :alt: Plot of the detachment time (our results and literature results) with respect to shear rate.            |
|     :align: center                                                                                                |
|     :name: Detachment times                                                                                       |
|                                                                                                                   |
|     Plot of the detachment times (our results and literature results) with respect to shear rate.                 |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+

Below are the plots of the contour of the bubble in the plane :math:`z = 0` when detachment occurs for different values of the shear rate.

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/bubble-detachment-contour.png                                                                 |
|     :alt: Cut of the contour of the bubble at detachment time in the plane :math:`z = 0` for different shear      |
|      rates.                                                                                                       |
|     :align: center                                                                                                |
|     :name: Contour cuts at detachment                                                                             |
|                                                                                                                   |
|     Cut of the contour of the bubble at detachment time in the plane :math:`z = 0` for different shear rates.     |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+


**make the comments on the figure better based on the content of the article and the memoire**

---------------------------
Possibilities for Extension
---------------------------

- **Non-newtonian case**: extend the case to a non-newtonian liquid. This should yield fairly different results because of all the shear around the bubble surface. 

- **FIND ANOTHER POSSIBLITY FOR EXTENSION** : blablablablabla

-----------
References
-----------


.. [#mirsandi2020] \H. Mirsandi, W. J. Smit, G. Kong, M. W. Baltussen, E. A. J. F. Peters, and J. A. M. Kuipers, ‘Bubble formation from an orifice in liquid cross-flow’, Chemical Engineering Journal, vol. 386, p. 120902, Apr. 2020, doi: 10.1016/j.cej.2019.01.181.

.. [#yue2007] \P. Yue, C. Zhou, and J. J. Feng, ‘Spontaneous shrinkage of drops and mass conservation in phase-field simulations’, Journal of Computational Physics, vol. 223, no. 1, pp. 1–9, Apr. 2007, doi: 10.1016/j.jcp.2006.11.020.

