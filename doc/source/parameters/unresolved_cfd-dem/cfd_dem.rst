***********************************************
CFD-DEM
***********************************************
This subsection includes parameters related to multiphase flow simulations using the both the gls_vans solver and the cfd-dem_coupling solver within Lethe.

.. code-block:: text

   subsection cfd-dem
      set grad div = false
      set void fraction time derivative = false
      set drag force = true
      set buoyancy force = true
      set shear force = true
      set pressure force = true
      set drag model = dallavalle
      set post processing = true
      set coupling frequency = 200
      set inlet boundary id = 2
      set outlet boundary id = 3
   end


* The ``grad div`` parameter allows the enabling of the grad div stabilization for the Volume Averaged Navier Stokes equations [1]. This allows for a much better mass conservation of the system.
* The ``void fraction time derivative`` parameters allows us to choose whether or not we want to account for the time derivative of the void fraction or take it equal to zero.
* The ``drag force``, ``buoyancy force``, ``shear force``, and ``pressure force`` parameters allow us to enable or disable the respective forces in a cfd-dem simulation.
* The ``drag model`` parameter allows one to choose the type of drag model to be implemented for the calculation of the drag force between the particles and the fluids. Available drag models at the time of writing are: Difelice [2], Rong [3], Dallavalle [4], and Koch and Hill [5].
* The ``post processing`` parameter, when enabled, allows the calculation of the pressure drop, void fraction in the packed region, and the mass conservation in a packed bed at each time step.
* The ``coupling frequency`` parameter is only applicable for the cfd-dem solver and it determines the number of DEM iterations per 1 CFD iteration.

.. note::
   The ``coupling frequency`` parameter is used to calculate the dem time step as it is not explicitly determined in the parameter file. It is calculated as: 

   .. math::
      \Delta t_{DEM} = \frac{\Delta t_{CFD}}{coupling frequency}

* The ``inlet boundary`` and ``outlet boundary`` parameters allow us to specify the ID of the inlet and outlet for pressure drop calculations.

[1] B. Blais, L. Barbeau, V. Bibeau, S. Gauvin, T. E. Geitani, S. Golshan, R. Kamble, G. Mirakhori, J. Chaouki, Lethe: An open-source parallel high- order adaptative cfd solver for incompressible flows, SoftwareX 12 (2020) 100579.

[2] R. Di Felice, The voidage function for fluid-particle interaction systems, International journal of multiphase flow 20 (1) (1994) 153–159.

[3] L. Rong, K. Dong, A. Yu, Lattice-boltzmann simulation of fluid flow through packed beds of uniform spheres: Effect of porosity, Chemical engineering science 99 (2013) 44–58.

[4] Sobieski, Wojciech. (2011). Drag Coefficient in Solid–Fluid System Modeling with the Eulerian Multiphase Model. Drying Technology. 29. 111-125. 10.1080/07373937.2010.482714. 

[5]  D. Jajcevic, E. Siegmann, C. Radeke, J. G. Khinast, Large-scale cfd–dem simulations of fluidized granular systems, Chemical Engineering Science 98 (2013) 298–310
