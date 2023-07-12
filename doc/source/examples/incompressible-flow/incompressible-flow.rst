****************************
Incompressible Flow
****************************

.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:
    :numbered:

    2d-lid‐driven-cavity-flow/lid‐driven-cavity-flow
    2d-flow-around-cylinder/2d-flow-around-cylinder
    2d-taylor-couette-flow/2d-taylor-couette-flow
    2d-backward-facing-step/2d-backward-facing-step
    2d-transient-around-ahmed-body/2d-transient-around-ahmed-body
    2d-transient-flow-around-cylinder/2d-transient-flow-around-cylinder
    2d-naca0012-low-reynolds/2d-naca0012-low-reynolds
    3d-flow-around-sphere/flow-around-sphere
    3d-flow-over-periodic-hills/3d-flow-over-periodic-hills
    2d-taylor-couette-flow-nitsche/2d-taylor-couette-flow-nitsche
    3d-mixer-using-single-rotating-frame/3d-mixer-using-single-rotating-frame
    3d-nitsche-mixer-with-pbt-impeller/nitsche-mixer-with-pbt-impeller

.. graphviz:: 

    digraph incompressible_diagram {
      graph [bgcolor="transparent", align=true];
      node [fontname=Arial, fontsize=18, shape=box, fontcolor=royalblue, color=royalblue];
      edge [color=royalblue];
      rankdir="LR";
      size = "9,9";

      incompressible_flow [label="Incompressible \nFlow", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/incompressible-flow.html"];

      incompressible_1 [label="2D"]; 

      incompressible_1_1 [label="1. Lid-driven cavity flow",href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-lid%E2%80%90driven-cavity-flow/lid%E2%80%90driven-cavity-flow.html"]; 

      incompressible_1_2 [label="2. Flow around a cylinder", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-flow-around-cylinder/2d-flow-around-cylinder.html"];

      incompressible_1_3 [label="3. Taylor-Couette flow", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-taylor-couette-flow/2d-taylor-couette-flow.html"];

      incompressible_1_4 [label="4. Flow past a backward-facing step", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-backward-facing-step/2d-backward-facing-step.html"];

      incompressible_1_5 [label="5. Transient flow around an Ahmed body", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-transient-around-ahmed-body/2d-transient-around-ahmed-body.html"];

      incompressible_1_6 [label="6. Transient flow around a cylinder", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-transient-flow-around-cylinder/2d-transient-flow-around-cylinder.html"];

      incompressible_1_7 [label="7. Flow around NACA0012 \nat low Reynolds number", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-naca0012-low-reynolds/2d-naca0012-low-reynolds.html"];

      incompressible_1_8 [label="10. Taylor-Couette flow using \nNitsche immersed boundary", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-taylor-couette-flow-nitsche/2d-taylor-couette-flow-nitsche.html"];

      incompressible_2 [label="3D",shape=polygon,sides=4]; 

      incompressible_2_1 [label="8. Flow around a sphere", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-flow-around-sphere/flow-around-sphere.html"];
      
      incompressible_2_2 [label="9. Flow over periodic hills", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-flow-over-periodic-hills/3d-flow-over-periodic-hills.html"];

      incompressible_2_3 [label="11. Ribbon mixer using a single \n rotating reference frame", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-mixer-using-single-rotating-frame/3d-mixer-using-single-rotating-frame.html"];

      incompressible_2_4 [label="12. Mixer with pitched-blade turbine impeller \nusing Nitsche immersed boundary", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-nitsche-mixer-with-pbt-impeller/nitsche-mixer-with-pbt-impeller.html"];

      incompressible_flow -> incompressible_1:w;
      incompressible_flow -> incompressible_2:w;
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