*************
Multiphysics
*************

.. toctree::
    :hidden:
    :maxdepth: 2
    :glob:

    vof
    heat-transfer
    tracer

.. graphviz:: 

    digraph multiphysics_diagram {
      graph [bgcolor="transparent", align=true, ranksep=1.5];
      node [fontname=Arial, fontsize=20, shape=box, fontcolor=royalblue, color=royalblue, height=1];
      edge [color=royalblue];
      rankdir="LR";
      size = "9,9";

      multiphysics [label="Multiphysics", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/multiphysics.html"];

      multiphysics_1 [label="VOF", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/vof.html"];  

      multiphysics_1_1 [label="Dam-Break", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/dam-break/dam-break.html"];

      multiphysics_1_2 [label="Rayleigh-Taylor Instability", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/rayleigh-taylor-instability/rayleigh-taylor-instability.html"];

      multiphysics_1_3 [label="Static Bubble", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/static-bubble/static-bubble.html"];

      multiphysics_1_4 [label="Rising Bubble", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/rising-bubble/rising-bubble.html"];

      multiphysics_1_5 [label="Sloshing in a Rectangular Tank", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/sloshing-in-rectangular-tank/sloshing-in-rectangular-tank.html"];

      multiphysics_1_6 [label="Capillary Wave", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/capillary-wave/capillary-wave.html"];

      multiphysics_1_7 [label="3D Dam-Break with an Obstacle", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/3d-dam-break/3d-dam-break.html"];

      multiphysics_1_8 [label="Water Injection in a Closed Cell", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/water-injection-in-a-closed-cell/water-injection-in-a-closed-cell.html"];

      multiphysics_1_9 [label="Air Bubble Compression", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/air-bubble-compression/air-bubble-compression.html"];

      multiphysics_2 [label="Heat Transfer", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/heat-transfer.html"];

      multiphysics_2_1 [label="Rayleigh-Bénard Convection", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/rayleigh-benard-convection/rayleigh-benard-convection.html"];

      multiphysics_2_2 [label="Stefan Problem: Melting of a Solid", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/stefan-problem/stefan-problem.html"];

      multiphysics_2_3 [label="Warming up a Viscous Fluid", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/warming-up-a-viscous-fluid/warming-up-a-viscous-fluid.html"];

      multiphysics_2_4 [label="Melting Cavity", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/melting-cavity/melting-cavity.html"];

      multiphysics_2_5 [label="Laser Heating", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/laser-heating/laser-heating.html"];

      multiphysics_2_6 [label="Laser Melt Pool", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/laser-melt-pool/laser-melt-pool.html"];

      multiphysics_2_7 [label="Concentric Heat Exchanger", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/concentric-heat-exchanger/concentric-heat-exchanger.html"];

      multiphysics_3 [label="Tracer", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/tracer.html"];  

      multiphysics_3_1 [label="Tracer through CAD Junction in Simplex", href="https://lethe-cfd.github.io/lethe/documentation/examples/multiphysics/tracer-through-cad-junction/tracer-through-cad-junction.html"];

      multiphysics -> multiphysics_1:w;
      multiphysics -> multiphysics_2:w;
      multiphysics -> multiphysics_3:w;

      multiphysics_1 -> multiphysics_1_1:w;
      multiphysics_1 -> multiphysics_1_2:w;
      multiphysics_1 -> multiphysics_1_3:w;
      multiphysics_1 -> multiphysics_1_4:w;
      multiphysics_1 -> multiphysics_1_5:w;
      multiphysics_1 -> multiphysics_1_6:w;
      multiphysics_1 -> multiphysics_1_7:w;
      multiphysics_1 -> multiphysics_1_8:w;
      multiphysics_1 -> multiphysics_1_9:w;

      multiphysics_2 -> multiphysics_2_1:w;
      multiphysics_2 -> multiphysics_2_2:w;
      multiphysics_2 -> multiphysics_2_3:w;
      multiphysics_2 -> multiphysics_2_4:w;
      multiphysics_2 -> multiphysics_2_5:w;
      multiphysics_2 -> multiphysics_2_6:w;
      multiphysics_2 -> multiphysics_2_7:w;

      multiphysics_3 -> multiphysics_3_1:w;
    }

