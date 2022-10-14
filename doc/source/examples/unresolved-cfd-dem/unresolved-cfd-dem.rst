****************************
Unresolved CFD-DEM
****************************

This section includes examples related to multiphase fluid-solid flows. We organize the examples from the single phase flow in porous media (Packed Bed Example) to multiphase flows (Fluidized Solid-Gas Beds). The packed bed example uses the gls_vans_3d solver which solves the Volume Average Navier Stokes (VANS) equations. The fluidized bed, the spouted bed and the Boycott effect example use the cfd_dem_coupling_3d solver which solves the VANS equations for the fluid phase coupled with the DEM equations for the solid phase.

.. toctree::
    :maxdepth: 1
    :glob:
    :numbered:


    cylindrical-packed-bed/cylindrical-packed-bed
    gas-solid-fluidized-bed/gas-solid-fluidized-bed
    gas-solid-spouted-bed/gas-solid-spouted-bed
    boycott-effect/boycott-effect
