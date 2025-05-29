******************
DEM-Multiphysics
******************

We organize the DEM-Multiphysics examples from the simplest to the most complicated example: 

.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:

    thermal-lines/thermal-lines

.. graphviz:: 

    digraph dem_mp_diagram {
      graph [bgcolor="transparent", align=true, ranksep=1.5];
      node [fontname=Arial, fontsize=15, shape=box, fontcolor=royalblue, color=royalblue];
      edge [color=royalblue];
      rankdir="LR";
      size = "9,9";
      
      dem_mp [label="DEM-Multiphysics", href="https://chaos-polymtl.github.io/lethe/documentation/examples/dem_mp/dem_mp.html"];

      dem_mp_1 [label="Thermal Lines", href="https://chaos-polymtl.github.io/lethe/documentation/examples/dem_mp/thermal-lines/thermal-lines.html"];

      dem_mp -> dem_mp_1:w;
   
    }