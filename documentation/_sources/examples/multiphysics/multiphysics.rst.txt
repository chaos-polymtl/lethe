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
    chns

.. graphviz:: 

    digraph multiphysics_diagram {
      graph [bgcolor="transparent", align=true, ranksep=1.5];
      node [fontname=Arial, fontsize=20, shape=box, fontcolor=royalblue, color=royalblue, height=1];
      edge [color=royalblue];
      rankdir="LR";
      size = "9,9";

      multiphysics [label="Multiphysics", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/multiphysics.html"];

      multiphysics_1 [label="VOF", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/vof.html"];  

      multiphysics_1_1 [label="Dam-Break", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/dam-break/dam-break.html"];

      multiphysics_1_2 [label="Rayleigh-Taylor Instability", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/rayleigh-taylor-instability/rayleigh-taylor-instability.html"];

      multiphysics_1_3 [label="Static Bubble", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/static-bubble/static-bubble.html"];

      multiphysics_1_4 [label="Rising Bubble", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/rising-bubble/rising-bubble.html"];

      multiphysics_1_5 [label="Sloshing in a Rectangular Tank", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/sloshing-in-rectangular-tank/sloshing-in-rectangular-tank.html"];

      multiphysics_1_6 [label="Capillary Wave", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/capillary-wave/capillary-wave.html"];

      multiphysics_1_7 [label="Rayleigh-Plateau Instability", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/rayleigh-plateau-instability/rayleigh-plateau-instability.html"];

      multiphysics_1_8 [label="3D Dam-Break with an Obstacle", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/3d-dam-break/3d-dam-break.html"];

      multiphysics_1_9 [label="Water Injection in a Closed Cell", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/water-injection-in-a-closed-cell/water-injection-in-a-closed-cell.html"];

      multiphysics_1_10 [label="Air Bubble Compression", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/air-bubble-compression/air-bubble-compression.html"];

      multiphysics_2 [label="Heat Transfer", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/heat-transfer.html"];

      multiphysics_2_1 [label="Cooling fin", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/cooling-fin/cooling-fin.html"];

      multiphysics_2_2 [label="Rayleigh-Bénard Convection", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/rayleigh-benard-convection/rayleigh-benard-convection.html"];

      multiphysics_2_3 [label="Stefan Problem: Melting of a Solid", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/stefan-problem/stefan-problem.html"];

      multiphysics_2_4 [label="Warming up a Viscous Fluid", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/warming-up-a-viscous-fluid/warming-up-a-viscous-fluid.html"];

      multiphysics_2_5 [label="Melting Cavity", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/melting-cavity/melting-cavity.html"];
      
      multiphysics_2_6 [label="Static Irradiation of a Bare Plate", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/static-irradiation/static-irradiation.html"];

      multiphysics_2_7 [label="Concentric Heat Exchanger", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/concentric-heat-exchanger/concentric-heat-exchanger.html"];

      multiphysics_2_8 [label="Cooling of a cylinder", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/cylinder-cooling/cylinder-cooling.html"];

      multiphysics_3 [label="Tracer", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/tracer.html"];  

      multiphysics_3_1 [label="Tracer in Static Mixer", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/tracer-in-static-mixer/tracer-in-static-mixer.html"];
      
      multiphysics_4 [label="CHNS", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/chns.html"];

      multiphysics_4_1 [label="Jurin's Law", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/jurins-law/jurins-law.html"];

      multiphysics_4_2 [label="Shear Bubble Detachment", href="https://chaos-polymtl.github.io/lethe/documentation/examples/multiphysics/bubble-detachment-shear-flow/bubble-detachment-shear-flow.html"];


      multiphysics -> multiphysics_1:w;
      multiphysics -> multiphysics_2:w;
      multiphysics -> multiphysics_3:w;
      multiphysics -> multiphysics_4:w;

      multiphysics_1 -> multiphysics_1_1:w;
      multiphysics_1 -> multiphysics_1_2:w;
      multiphysics_1 -> multiphysics_1_3:w;
      multiphysics_1 -> multiphysics_1_4:w;
      multiphysics_1 -> multiphysics_1_5:w;
      multiphysics_1 -> multiphysics_1_6:w;
      multiphysics_1 -> multiphysics_1_7:w;
      multiphysics_1 -> multiphysics_1_8:w;
      multiphysics_1 -> multiphysics_1_9:w;
      multiphysics_1 -> multiphysics_1_10:w;

      multiphysics_2 -> multiphysics_2_1:w;
      multiphysics_2 -> multiphysics_2_2:w;
      multiphysics_2 -> multiphysics_2_3:w;
      multiphysics_2 -> multiphysics_2_4:w;
      multiphysics_2 -> multiphysics_2_5:w;
      multiphysics_2 -> multiphysics_2_6:w;
      multiphysics_2 -> multiphysics_2_7:w;
      multiphysics_2 -> multiphysics_2_8:w;

      multiphysics_3 -> multiphysics_3_1:w;

      multiphysics_4 -> multiphysics_4_1:w;
      multiphysics_4 -> multiphysics_4_2:w;
    }

