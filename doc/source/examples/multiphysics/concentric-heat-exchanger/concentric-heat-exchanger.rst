====================================
Concentric heat exchanger
====================================

This example simulates an heat exchanger which is made of two concentric pipes in which a hot and a cold fluid circulate in a counter-current fashion. In this problem, we simulate the fluid flow within two fluid regions as well as the heat transfer within the entire domain (solid and fluid).

----------------------------------
Features
----------------------------------
- Solver: ``gls_navier_stokes_3d`` 
- Heat transfer pĥysics
- Conjugated heat transfer

Files used in this example
---------------------------
Prm file: ``examples/multiphysics/concentric-heat-exchanger/concentric-heat-exchanger.prm``
GMSH file: ``examples/multiphysics/concentric-heat-exchanger/concentric-cylinders.geo``
GMSH mesh: ``examples/multiphysics/concentric-heat-exchanger/concentric-cylinders.msh``


Description of the case
-------------------------

Heat exchangers are a common unit operations used in many type of industries to transfer energy from one fluid to another. In this case, we simulate the most simple heat exchanger geometry, which is a cocentric tube in which the hot fluid circulates within the outer region and a cold fluid circulates within the inner region. We model the full heat transfer by simulating the motion of the fluid in both regions and the heat transfer within the entire domain. 

We consider steel concentric tubes with radii of :math:`R_0=2\text{mm} ,R_1=5\text{mm},R_2=7.5\text{mm}` in which water circulates. We consider a counter-current flow with an inner tube velocity of :math:`u_i=10\text{mm/s}` and a outer tube velocity of :math:`u_o=2\text{mm/s}`. We will not formulate the problem in SI units, but instead we will express the fundamental length in mm. This will ensure that most variables are of order 1 and will lead to a system matrix that will have a better condition number.

.. 
  image:: images/stefan-problem-illustration.png
    :alt: problem_illustration
    :align: center



We will compare the results we obtain with the CFD simulations with results that would be obtained using a correlation. Using correlations, we can estimate the heat transfer coefficients and the temperature profile to be:

The analytical solution for the temperature in the liquid is given by:

Lethe has the capability to solve 


--------------
Parameter file
--------------

We first define the geometry in which the simulation is carried out using the mesh subsection:

.. code-block:: text

    #---------------------------------------------------
    # Mesh
    #---------------------------------------------------
    subsection mesh
      set type               = dealii
      set grid type          = subdivided_hyper_rectangle
      set grid arguments     = 100, 1 : 0, 0 : 1, 0.1 : true
      set initial refinement = 0
    end

We use the ``dealii`` GridGenerator to generate a ``subdivided_hyper_rectangle``. This rectangle contains 100 cells in the x direction and 1 in the y direction. It is created from two points, :math:`(0,0)` and :math:`(1,0.1)`. Finally, we give a different id to each boundary of the domain, hence the colorize option is set to true.

The next step is establishing the boundary conditions:

.. code-block:: text

    # --------------------------------------------------
    # Boundary Conditions
    #---------------------------------------------------
    subsection boundary conditions heat transfer
      set number = 1
      subsection bc 0
        set id    = 0
        set type  = temperature
        set value = 1
      end
    end

Note that we only set one boundary condition for the temperature, which is a constant temperature on the wall which bears the ID 0 (the left wall). By default, boundaries on which boundary conditions are not specified are no-flux Neumann boundary conditions which, for a heat transfer problem results in:

.. math::
  \nabla T \cdot \mathbf{n} = 0

Next, we define the physical properties:

.. code-block:: text

    subsection physical properties
      set number of fluids = 1
      subsection fluid 0
        set thermal conductivity model = constant
        set thermal conductivity       = 1
    
        set specific heat model = phase_change
        subsection phase change
          # Enthalpy of the phase change
          set latent enthalpy = 100
    
          # Temperature of the liquidus
          set liquidus temperature = 0.02
    
          # Temperature of the solidus
          set solidus temperature = 0
    
          # Specific heat of the liquid phase
          set specific heat liquid = 1
    
          # Specific heat of the solid phase
          set specific heat solid = 1
        end
      end
    end

This subsection defines the various parameters of the specific heat model for phase change. Key parameters to note are the solidus and liquidus temperatures. These parameters define the phase change interval, that is the temperature interval over which the phase change occurs. For pure substance, this interval should, in theory, be infinitely small. However, this leads to a numerically unstable solution. Consequently, we set a finite value which should be relatively small, but not too small as to lead to numerical instabilities. In the present case, we set this interval to 0.02C, which is sufficient to guarantee a high degree of accuracy while maintaining numerical stability. The impact of this parameter on the stability and the accuracy of the model has been studied in depth by `Blais & Ilinca (2018)`_.

Finally, the only remaining section is the simulation control, which controls the flow of the simulation. We simulate until a :math:`t=5s` using a time step of :math:`\Delta t=0.02s` using a BDF1 (implicit Euler) time integration scheme and we output the solution at every iteration.

.. code-block:: text

    # --------------------------------------------------
    # Simulation Control
    #---------------------------------------------------
    subsection simulation control
      set method           = bdf1
      set output frequency = 1
      set output name      = stefan
      set output path      = ./output/
      set time end         = 5
      set time step        = 0.02
    end



-------
Results
-------

The following image compares the results obtained with Lethe with the analytical solution for the Stefan problem at :math:`t=5`. This data is extracted through the use of a python script available in the folder of the example. We see that a quasi perfect agreement can be obtained with the analytical solution of the Stefan problem. 

.. image:: images/lethe-stefan-comparison.png
    :alt: comparison_analytical_solution
    :align: center

Refining the mesh, decreasing the time step and decreasing the phase change interval (by decreasing ``liquidus temperature``) would increase the accuracy of the solution since the analytical solution of the Stefan problem is defined for a pure fluid (for which the liquid and the solidus temperatures are equal).


Possibilities for extension
----------------------------

- **Consider different Stefan numbers:** The solver in Lethe is sufficiently robust to simulate a large range of Stefan numbers. You can try to simulate the problem with different Stefan number and see how the value of the Stefan number affects the solution.

- **Simulate a more complex geometry:** The phase change model can be readily used in any sort of geometry using, for example, a simplex mesh. An easy extension of this problem is to consider any 2D or 3D geometry.

----------------------------
References
----------------------------

`[1] <https://doi.org/10.1016/j.applthermaleng.2007.01.008>`_ aus der Wiesche, Stefan. "Numerical heat transfer and thermal engineering of AdBlue (SCR) tanks for combustion engine emission reduction." Applied Thermal Engineering 27.11-12 (2007): 1790-1798.
