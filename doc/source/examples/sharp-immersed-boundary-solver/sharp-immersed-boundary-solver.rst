******************************
Sharp Immersed Boundary Solver
******************************

.. toctree::
    :hidden:
    :maxdepth: 2
    :glob:
    :numbered:

    resolved-cfd-dem
    incompressible-flow

.. graphviz:: 

    digraph sharp_diagram {
      graph [align=true];
      node [fontname=Arial, fontsize=18, shape=box];
      center=true;
      fontname="Comic Sans MS";
      fontsize = "20";
      rankdir="LR";
      size = "9,9";
      
      sharp_immersed_boundary_solver [label="Sharp Immersed \nBoundary Solver", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/sharp-immersed-boundary-solver.html"];

      sharp_1 [label="1. Resolved CFD-DEM", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/resolved-cfd-dem.html"];   

      sharp_1_1 [label="1.1. Sedimentation of one particle", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/sedimentation-1-particle/sedimentation-1-particle.html"]; 

      sharp_1_2 [label="1.2. Sedimentation of 64 particles", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/sedimentation-64-particles/sedimentation-64-particles.html"]; 

      sharp_2 [label="2. Incompressible flow", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/incompressible-flow.html"];  

      sharp_2_1 [label="2.1. Flow around a cylinder using \nthe sharp interface method", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/cylinder-with-sharp-interface/cylinder-with-sharp-interface.html"]; 

      sharp_2_2 [label="2.2. Non-Newtonian flow \npast a sphere", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/sphere-carreau-with-sharp-interface/sphere-carreau-with-sharp-interface.html"]; 

      sharp_2_3 [label="2.3. 3D Mixer with pitched-blade \nturbine impeller using Composite \nSharp-immersed boundary", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/3d-composite-mixer-with-pbt-impeller/3d-composite-mixer-with-pbt-impeller.html"]; 

      sharp_2_4 [label="2.4. 3D Mixer with pitched-blade turbine \nimpeller using OpenCascade \nSharp-immersed boundary", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/3d-opencascade-mixer-with-pbt-impeller/3d-opencascade-mixer-with-pbt-impeller.html"]; 

      sharp_immersed_boundary_solver -> sharp_1
      sharp_immersed_boundary_solver -> sharp_2

      sharp_1 -> sharp_1_1:w;
      sharp_1 -> sharp_1_2:w;
      sharp_2 -> sharp_2_1:w;
      sharp_2 -> sharp_2_2:w;
      sharp_2 -> sharp_2_3:w;
      sharp_2 -> sharp_2_4:w;
    }