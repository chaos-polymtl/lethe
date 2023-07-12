****************************
Incompressible flow
****************************

.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:

    2d-lid‐driven-cavity-flow/lid‐driven-cavity-flow
    2d-flow-around-cylinder/2d-flow-around-cylinder
    2d-taylor-couette-flow/2d-taylor-couette-flow
    2d-backward-facing-step/2d-backward-facing-step
    2d-transient-around-ahmed-body/2d-transient-around-ahmed-body
    2d-transient-flow-around-cylinder/2d-transient-flow-around-cylinder
    2d-naca0012-low-reynolds/2d-naca0012-low-reynolds
    2d-taylor-couette-flow-nitsche/2d-taylor-couette-flow-nitsche
    3d-flow-around-sphere/flow-around-sphere
    3d-flow-over-periodic-hills/3d-flow-over-periodic-hills
    3d-mixer-using-single-rotating-frame/3d-mixer-using-single-rotating-frame
    3d-nitsche-mixer-with-pbt-impeller/nitsche-mixer-with-pbt-impeller

.. graphviz:: 

    digraph incompressible_diagram {
      graph [bgcolor="transparent", align=true, ranksep=1.5];
      node [fontname=Arial, fontsize=18, shape=box, fontcolor=royalblue, color=royalblue];
      edge [color=royalblue];
      rankdir="LR";
      size = "9,9";

      {
      node [color=none];
      image_1 [image="./map_images/3d-mixer-using-single-rotating-frame.png", label="", fixedsize=true, width=2, height=3]
      }

      incompressible_flow [label="Incompressible \nflow", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/incompressible-flow.html"];

      {
      node [color=none];
      image_2 [image="./map_images/lid-driven-cavity.png", label="", fixedsize=true, width=3, height=3]
      }
      
      incompressible_1 [label="2D"]; 

      incompressible_1_1 [label="Lid-driven cavity flow",href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-lid%E2%80%90driven-cavity-flow/lid%E2%80%90driven-cavity-flow.html"]; 

      incompressible_1_2 [label="Flow around a cylinder", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-flow-around-cylinder/2d-flow-around-cylinder.html"];

      incompressible_1_3 [label="Taylor-Couette flow", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-taylor-couette-flow/2d-taylor-couette-flow.html"];

      incompressible_1_4 [label="Flow past a backward-facing step", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-backward-facing-step/2d-backward-facing-step.html"];

      incompressible_1_5 [label="Transient flow around an Ahmed body", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-transient-around-ahmed-body/2d-transient-around-ahmed-body.html"];

      incompressible_1_6 [label="Transient flow around a cylinder", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-transient-flow-around-cylinder/2d-transient-flow-around-cylinder.html"];

      incompressible_1_7 [label="Flow around NACA0012 \nat low Reynolds number", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-naca0012-low-reynolds/2d-naca0012-low-reynolds.html", tooltip="Flow around NACA0012 at low Reynolds number"];

      incompressible_1_8 [label="Taylor-Couette flow using \nNitsche immersed boundary", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-taylor-couette-flow-nitsche/2d-taylor-couette-flow-nitsche.html", tooltip="Taylor-Couette flow using Nitsche immersed boundary"];

      incompressible_2 [label="3D",shape=polygon,sides=4]; 

      incompressible_2_1 [label="Flow around a sphere", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-flow-around-sphere/flow-around-sphere.html"];
      
      incompressible_2_2 [label="Flow over periodic hills", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-flow-over-periodic-hills/3d-flow-over-periodic-hills.html"];

      incompressible_2_3 [label="Ribbon mixer using a single \n rotating reference frame", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-mixer-using-single-rotating-frame/3d-mixer-using-single-rotating-frame.html", tooltip="Ribbon mixer using a single rotating reference frame"];

      incompressible_2_4 [label="Mixer with pitched-blade turbine impeller \nusing Nitsche immersed boundary", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-nitsche-mixer-with-pbt-impeller/nitsche-mixer-with-pbt-impeller.html", tooltip="Mixer with pitched-blade turbine impeller using Nitsche immersed boundary"];

      incompressible_flow:e -> incompressible_1:w;
      incompressible_flow:e -> incompressible_2:w;

      incompressible_1 -> incompressible_1_1:w;
      incompressible_1 -> incompressible_1_2:w;
      incompressible_1 -> incompressible_1_3:w;
      incompressible_1 -> incompressible_1_4:w;
      incompressible_1 -> incompressible_1_5:w;
      incompressible_1 -> incompressible_1_6:w;
      incompressible_1 -> incompressible_1_7:w;
      incompressible_1 -> incompressible_1_8:w;
      incompressible_2 -> incompressible_2_1:w;
      incompressible_2 -> incompressible_2_2:w;
      incompressible_2 -> incompressible_2_3:w;
      incompressible_2 -> incompressible_2_4:w;

    }