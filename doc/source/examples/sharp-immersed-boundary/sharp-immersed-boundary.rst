******************************
Sharp Immersed Boundary Solver
******************************

.. toctree::
    :hidden:
    :maxdepth: 2
    :glob:

    resolved-cfd-dem
    incompressible-flow
    geometry-definition

.. graphviz:: 

    digraph sharp_diagram {
      graph [bgcolor="transparent", align=true];
      node [fontname=Arial, fontsize=20, shape=box, fontcolor=royalblue, color=royalblue, height=1];
      edge [color=royalblue];
      rankdir="LR";
      size = "9,9";
      
      sharp_immersed_boundary_solver [label="Sharp Immersed \nBoundary", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/sharp-immersed-boundary.html"];

      sharp_1 [label="Resolved CFD-DEM", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/resolved-cfd-dem.html"];

      sharp_1_1 [label="Sedimentation of One Particle", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/sedimentation-1-particle/sedimentation-1-particle.html"];

      sharp_1_2 [label="Sedimentation of 64 Particles", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/sedimentation-64-particles/sedimentation-64-particles.html"];

      sharp_2 [label="Incompressible Flow", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/incompressible-flow.html"];

      sharp_2_1 [label="Flow around a Cylinder Using \nthe Sharp Interface Method", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/cylinder-with-sharp-interface/cylinder-with-sharp-interface.html", tooltip="Flow around a cylinder using the sharp interface method"];

      sharp_2_2 [label="Non-Newtonian Flow \npast a Sphere", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/sphere-carreau-with-sharp-interface/sphere-carreau-with-sharp-interface.html", tooltip="Non-Newtonian flow past a sphere"];

      sharp_2_3 [label="3D Mixer with Pitched-Blade \nTurbine Impeller Using Composite \nSharp-Immersed Boundary", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/3d-composite-mixer-with-pbt-impeller/3d-composite-mixer-with-pbt-impeller.html", tooltip="3D Mixer with pitched-blade turbine impeller using Composite Sharp-immersed boundary"];

      sharp_2_4 [label="3D Mixer with Pitched-Blade Turbine \nImpeller Using OpenCascade \nSharp-Immersed Boundary", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/3d-opencascade-mixer-with-pbt-impeller/3d-opencascade-mixer-with-pbt-impeller.html", tooltip="3D Mixer with pitched-blade turbine impeller using OpenCascade Sharp-immersed boundary"];
      
      sharp_2_5 [label="Inclined 3D Mixer with Pitched-Blade \nTurbine Impeller Using Composite \nSharp-Immersed Boundary", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/inclined-3d-composite-mixer-with-pbt-impeller/inclined-3d-composite-mixer-with-pbt-impeller.html", tooltip="Inclined 3D mixer with pitched-blade turbine impeller using composite sharp-immersed boundary"];
      
      sharp_3 [label="Geometry Definition", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/geometry-definition.html"];
      
      sharp_3_1 [label="Simple Plane Model From Composite", href="https://lethe-cfd.github.io/lethe/documentation/examples/sharp-immersed-boundary/sharp-immersed-boundary/simple-plane-model-from-composite.html", tooltip="Simple Plane Model From Composite"];

      sharp_immersed_boundary -> sharp_1:w;
      sharp_immersed_boundary -> sharp_2:w;
      sharp_immersed_boundary -> sharp_3:w;

      sharp_1 -> sharp_1_1:w;
      sharp_1 -> sharp_1_2:w;
      sharp_2 -> sharp_2_1:w;
      sharp_2 -> sharp_2_2:w;
      sharp_2 -> sharp_2_3:w;
      sharp_2 -> sharp_2_4:w;
      sharp_2 -> sharp_2_5:w;
      sharp_3 -> sharp_3_1:w;
    }
