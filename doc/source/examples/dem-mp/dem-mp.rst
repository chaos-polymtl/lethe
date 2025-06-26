******************
DEM-Multiphysics
******************

We organize the DEM-Multiphysics examples from the simplest to the most complicated example: 

.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:

    heat-transfer-in-aligned-particles/heat-transfer-in-aligned-particles
    heated-packed-bed/heated-packed-bed

.. graphviz:: 

    digraph dem_mp_diagram {
      graph [bgcolor="transparent", align=true, ranksep=1.5];
      node [fontname=Arial, fontsize=15, shape=box, fontcolor=royalblue, color=royalblue];
      edge [color=royalblue];
      rankdir="LR";
      size = "9,9";
      
      dem_mp [label="DEM-Multiphysics", href="https://chaos-polymtl.github.io/lethe/documentation/examples/dem-mp/dem-mp.html"];

      dem_mp_1 [label="Heat Transfer in Aligned Particles", href="https://chaos-polymtl.github.io/lethe/documentation/examples/dem-mp/heat-transfer-in-aligned-particles/heat-transfer-in-aligned-particles.html"];

      dem_mp_2 [label="Heated Packed Bed", href="https://chaos-polymtl.github.io/lethe/documentation/examples/dem-mp/heated-packed-bed/heated-packed-bed.html"];

      dem_mp -> dem_mp_1:w;
      dem_mp -> dem_mp_2:w;
   
    }