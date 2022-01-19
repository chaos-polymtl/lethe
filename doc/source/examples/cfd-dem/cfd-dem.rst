****************************
CFD-DEM
****************************

This section includes examples related to multiphase flows. We organise the examples from the single phase flow in porous media (Packed Bed Example) to multiphase flows (Fluidized Solid-Gas Beds). The packed bed example uses the gls_vans_3d solver which solves the Volume Average Navier Stokes (VANS) equations. The fluidized bed example uses the cfd_dem_coupling_3d solver which solves the VANS equations for the fluid phase and the DEM equations for the solid phase. 

.. toctree::
    :maxdepth: 1
    :glob:
    :numbered:

    cylindrical-packed-bed/cylindrical-packed-bed
    square-fluidized-bed/square-fluidized-bed
