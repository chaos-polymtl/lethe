***********************************************
CFD-DEM
***********************************************
This subsection includes parameters related to multiphase flow simulations using the both the gls_vans solver and the cfd-dem solver within Lethe.

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

* The `grad div` parameter allows the enabling of the grad div stabilization for the Volume Averaged Navier Stokes equations. This allows for a much better mass conservation of the system.
* The `void fraction time derivative` parameters allows us to choose whether or not we want to account for the time derivative of the void fraction or take it equal to zero.
* The `drag force`, `buoyancy force`, `shear force`, and `pressure force` parameters allows us to enable or disable the respective forces in a cfd-dem simulation.
* The `drag model` parameter allows to choose the type of drag model to be implemented for the calculation of the drag force between the particles and the fluids. Available drag models at the time of writing are: Difelice, Rong, and Dallavalle.
* The `post processing` parameter, when enabled allows the calculation of the pressure drop, void fraction in the packed region, and the mass conservation in a packed bed at each time step.
* The `coupling frequency` parameter is only applicable for the cfd-dem solver and it determines the number of DEM iterations per 1 CFD iteration.
* The `inlet boundary` and `outlet boundary` parameters allows us to specify the ID of the inlet and outlet for pressure drop calculations.
