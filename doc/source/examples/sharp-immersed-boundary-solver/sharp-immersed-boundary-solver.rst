******************************
Sharp immersed boundary solver
******************************

.. toctree::
    :hidden:
    :maxdepth: 2
    :glob:

    resolved-cfd-dem
    incompressible-flow

.. graphviz:: 

    digraph sharp_diagram {
      graph [bgcolor="transparent", align=true, ranksep=1.5];
      node [fontname=Arial, fontsize=18, shape=box, fontcolor=royalblue, color=royalblue];
      edge [color=royalblue];
      rankdir="LR";
      size = "9,9";
      
      sharp_immersed_boundary_solver [label="Sharp Immersed \nBoundary Solver", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/sharp-immersed-boundary-solver.html"];

      sharp_1 [label="Resolved CFD-DEM", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/resolved-cfd-dem.html"];   

      sharp_1_1 [label="Sedimentation of one particle", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/sedimentation-1-particle/sedimentation-1-particle.html"]; 

      sharp_1_2 [label="Sedimentation of 64 particles", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/sedimentation-64-particles/sedimentation-64-particles.html"]; 

      sharp_2 [label="Incompressible flow", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/incompressible-flow.html"];  

      sharp_2_1 [label="Flow around a cylinder using \nthe sharp interface method", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/cylinder-with-sharp-interface/cylinder-with-sharp-interface.html", tooltip="Flow around a cylinder using the sharp interface method"]; 

      sharp_2_2 [label="Non-Newtonian flow \npast a sphere", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/sphere-carreau-with-sharp-interface/sphere-carreau-with-sharp-interface.html", tooltip="Non-Newtonian flow past a sphere"]; 

      sharp_2_3 [label="3D Mixer with pitched-blade \nturbine impeller using Composite \nSharp-immersed boundary", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/3d-composite-mixer-with-pbt-impeller/3d-composite-mixer-with-pbt-impeller.html", tooltip="3D Mixer with pitched-blade turbine impeller using Composite Sharp-immersed boundary"]; 

      sharp_2_4 [label="3D Mixer with pitched-blade turbine \nimpeller using OpenCascade \nSharp-immersed boundary", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/3d-opencascade-mixer-with-pbt-impeller/3d-opencascade-mixer-with-pbt-impeller.html", tooltip="3D Mixer with pitched-blade turbine impeller using OpenCascade Sharp-immersed boundary"]; 

      sharp_immersed_boundary_solver -> sharp_1:w;
      sharp_immersed_boundary_solver -> sharp_2:w;

      sharp_1 -> sharp_1_1:w;
      sharp_1 -> sharp_1_2:w;
      sharp_2 -> sharp_2_1:w;
      sharp_2 -> sharp_2_2:w;
      sharp_2 -> sharp_2_3:w;
      sharp_2 -> sharp_2_4:w;
    }