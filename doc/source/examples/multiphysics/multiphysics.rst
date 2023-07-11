*********************
Multiphysics Examples
*********************

.. toctree::
    :maxdepth: 2
    :glob:
    :numbered:

    vof
    heat-transfer
    tracer

.. graphviz:: 

    digraph multiphysics_diagram {
      graph [align=true];
      node [fontname=Arial, fontsize=18];
      center=true;
      fontname="Comic Sans MS";
      fontsize = "20";
      rankdir="LR";
      size = "9,9";
      
      multiphysics [label="Multiphysics", href="https://lethe-cfd.github.io/lethe/examples/multiphysics/multiphysics.html"];

      multiphysics_1 [label="1. VOF",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/vof.html"];  

      multiphysics_1_1 [label="1.1. Dam-break",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/dam-break/dam-break.html"];  

      multiphysics_1_2 [label="1.2. Rayleigh-Taylor instability",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/rayleigh-taylor-instability/rayleigh-taylor-instability.html"];  

      multiphysics_1_3 [label="1.3. Static bubble",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/static-bubble/static-bubble.html"];  

      multiphysics_1_4 [label="1.4. Rising bubble",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/rising-bubble/rising-bubble.html"];  

      multiphysics_1_5 [label="1.5. Sloshing in a rectangular tank",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/sloshing-in-rectangular-tank/sloshing-in-rectangular-tank.html"];  

      multiphysics_1_6 [label="1.6. 3D Dam Break with an Obstacle",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/3d-dam-break/3d-dam-break.html"];  

      multiphysics_2 [label="2. Heat transfer",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/heat-transfer.html"];  

      multiphysics_2_1 [label="2.1. Rayleigh-Bénard Convection",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/rayleigh-benard-convection/rayleigh-benard-convection.html"]; 

      multiphysics_2_2 [label="2.2. Stefan problem: melting of a solid",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/stefan-problem/stefan-problem.html"]; 

      multiphysics_2_3 [label="2.3. Warming Up a Viscous Fluid",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/warming-up-a-viscous-fluid/warming-up-a-viscous-fluid.html"]; 

      multiphysics_2_4 [label="2.4. Melting Cavity",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/melting-cavity/melting-cavity.html"]; 

      multiphysics_2_5 [label="2.5. Laser heating",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/laser-heating/laser-heating.html"]; 

      multiphysics_2_6 [label="2.6. Laser meltpool",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/laser-meltpool/laser-meltpool.html"]; 

      multiphysics_2_7 [label="2.7. Concentric heat exchanger",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/concentric-heat-exchanger/concentric-heat-exchanger.html"]; 

      multiphysics_3 [label="3. Tracer",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/tracer.html"];  

      multiphysics_3_1 [label="3.1. Tracer through CAD Junction in Simplex",shape=polygon,sides=4, href="https://lethe-cfd.github.io/lethe/examples/multiphysics/tracer-through-cad-junction/tracer-through-cad-junction.html"]; 

      multiphysics -> multiphysics_1;
      multiphysics -> multiphysics_2;
      multiphysics -> multiphysics_3;

      multiphysics_1 -> multiphysics_1_1;
      multiphysics_1 -> multiphysics_1_2;
      multiphysics_1 -> multiphysics_1_3;
      multiphysics_1 -> multiphysics_1_4;
      multiphysics_1 -> multiphysics_1_5;
      multiphysics_1 -> multiphysics_1_6;

      multiphysics_2 -> multiphysics_2_1;
      multiphysics_2 -> multiphysics_2_2;
      multiphysics_2 -> multiphysics_2_3;
      multiphysics_2 -> multiphysics_2_4;
      multiphysics_2 -> multiphysics_2_5;
      multiphysics_2 -> multiphysics_2_6;
      multiphysics_2 -> multiphysics_2_7;

      multiphysics_3 -> multiphysics_3_1;
    }

