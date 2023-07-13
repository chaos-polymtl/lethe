****************************
Unresolved CFD-DEM
****************************

This section includes examples related to multiphase fluid-solid flows. We organize the examples from the single phase flow in porous media (Packed Bed Example) to multiphase flows (Fluidized Solid-Gas Beds). The packed bed example uses the gls_vans_3d solver which solves the Volume Average Navier Stokes (VANS) equations. The fluidized bed, the spouted bed and the Boycott effect example use the cfd_dem_coupling_3d solver which solves the VANS equations for the fluid phase coupled with the DEM equations for the solid phase.

.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:

    cylindrical-packed-bed/cylindrical-packed-bed
    gas-solid-fluidized-bed/gas-solid-fluidized-bed
    gas-solid-spouted-bed/gas-solid-spouted-bed
    boycott-effect/boycott-effect

.. graphviz:: 

    digraph unresolved_cfd_dem_diagram {
      graph [bgcolor="transparent", align=true, ranksep=1.5];
      node [fontname=Arial, fontsize=18, shape=box, fontcolor=royalblue, color=royalblue];
      edge [color=royalblue];
      rankdir="LR";
      size = "9,9";

      unresolved_cfd_dem [label="Unresolved \nCFD-DEM", href="https://lethe-cfd.github.io/lethe/examples/unresolved-cfd-dem/unresolved-cfd-dem.html"];

      cfd_dem_1 [label="Cylindrical packed bed", href="https://lethe-cfd.github.io/lethe/examples/unresolved-cfd-dem/cylindrical-packed-bed/cylindrical-packed-bed.html"]; 

      cfd_dem_2 [label="Gas-solid fluidized bed", href="https://lethe-cfd.github.io/lethe/examples/unresolved-cfd-dem/gas-solid-fluidized-bed/gas-solid-fluidized-bed.html"]; 

      cfd_dem_3 [label="Gas-solid spouted bed", href="https://lethe-cfd.github.io/lethe/examples/unresolved-cfd-dem/gas-solid-spouted-bed/gas-solid-spouted-bed.html"]; 

      cfd_dem_4 [label="Boycott effect", href="https://lethe-cfd.github.io/lethe/examples/unresolved-cfd-dem/boycott-effect/boycott-effect.html"]; 
      
      unresolved_cfd_dem -> cfd_dem_1:w;
      unresolved_cfd_dem -> cfd_dem_2:w;
      unresolved_cfd_dem -> cfd_dem_3:w;
      unresolved_cfd_dem -> cfd_dem_4:w;
    }