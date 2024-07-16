****************************
Incompressible Flow
****************************

.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:

    2d-lid-driven-cavity-flow/lid-driven-cavity-flow
    2d-flow-around-cylinder/2d-flow-around-cylinder
    2d-taylor-couette-flow/2d-taylor-couette-flow
    2d-backward-facing-step/2d-backward-facing-step
    2d-transient-flow-around-ahmed-body/2d-transient-flow-around-ahmed-body
    2d-transient-flow-around-cylinder/2d-transient-flow-around-cylinder
    2d-naca0012-low-reynolds/2d-naca0012-low-reynolds
    2d-taylor-couette-flow-nitsche/2d-taylor-couette-flow-nitsche
    3d-flow-around-sphere/flow-around-sphere
    3d-taylor-green-vortex/3d-taylor-green-vortex
    3d-turbulent-taylor-couette/3d-turbulent-taylor-couette
    3d-flow-over-periodic-hills/3d-flow-over-periodic-hills
    3d-mixer-using-single-rotating-frame/3d-mixer-using-single-rotating-frame
    3d-nitsche-mixer-with-pbt-impeller/nitsche-mixer-with-pbt-impeller

.. graphviz:: 

    digraph incompressible_diagram {
      graph [bgcolor="transparent", align=true, ranksep=1.5];
      node [fontname=Arial, fontsize=20, shape=box, fontcolor=royalblue, color=royalblue,height=1];
      edge [color=royalblue];
      rankdir="LR";
      size = "9,9";

      {
      node [color=none];
      image_1 [image="./map_images/3d-mixer-using-single-rotating-frame.png", label="", fixedsize=true, width=2, height=3]
      }

      incompressible_flow [label="Incompressible \nFlow", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/incompressible-flow.html"];

      {
      node [color=none];
      image_2 [image="./map_images/lid-driven-cavity.png", label="", fixedsize=true, width=3, height=3]
      }
      
      incompressible_1 [label="2D"]; 

      incompressible_1_1 [label="Lid-Driven Cavity Flow",href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/2d-lid-driven-cavity-flow/lid-driven-cavity-flow.html"];

      incompressible_1_2 [label="Flow around a Cylinder", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/2d-flow-around-cylinder/2d-flow-around-cylinder.html"];

      incompressible_1_3 [label="Taylor-Couette Flow", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/2d-taylor-couette-flow/2d-taylor-couette-flow.html"];

      incompressible_1_4 [label="Flow past a Backward-Facing Step", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/2d-backward-facing-step/2d-backward-facing-step.html"];

      incompressible_1_5 [label="Transient Flow around an Ahmed Body", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/2d-transient-around-ahmed-body/2d-transient-around-ahmed-body.html"];

      incompressible_1_6 [label="Transient Flow around a Cylinder", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/2d-transient-flow-around-cylinder/2d-transient-flow-around-cylinder.html"];

      incompressible_1_7 [label="Flow around NACA0012 \nat low Reynolds Number", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/2d-naca0012-low-reynolds/2d-naca0012-low-reynolds.html", tooltip="Flow around NACA0012 at low Reynolds number"];

      incompressible_1_8 [label="Taylor-Couette Flow Using \nNitsche Immersed Boundary", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/2d-taylor-couette-flow-nitsche/2d-taylor-couette-flow-nitsche.html", tooltip="Taylor-Couette flow using Nitsche immersed boundary"];

      incompressible_2 [label="3D",shape=polygon,sides=4]; 

      incompressible_2_1 [label="Flow around a Sphere", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/3d-flow-around-sphere/flow-around-sphere.html"];
      
      incompressible_2_2 [label="Taylor-Green Vortex", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/3d-taylor-green-vortex/3d-taylor-green-vortex.html"];

      incompressible_2_3 [label="Turbulent Taylor-Couette", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/3d-turbulent-taylor-couette/3d-turbulent-taylor-couette.html"];

      incompressible_2_4 [label="Flow over Periodic Hills", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/3d-flow-over-periodic-hills/3d-flow-over-periodic-hills.html"];

      incompressible_2_5 [label="Ribbon Mixer Using a Single \n Rotating Reference Frame", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/3d-mixer-using-single-rotating-frame/3d-mixer-using-single-rotating-frame.html", tooltip="Ribbon mixer using a single rotating reference frame"];

      incompressible_2_6 [label="Mixer with Pitched-Blade Turbine Impeller \nUsing Nitsche Immersed Boundary", href="https://chaos-polymtl.github.io/lethe/documentation/examples/incompressible-flow/3d-nitsche-mixer-with-pbt-impeller/nitsche-mixer-with-pbt-impeller.html", tooltip="Mixer with pitched-blade turbine impeller using Nitsche immersed boundary"];

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
      incompressible_2 -> incompressible_2_5:w;
      incompressible_2 -> incompressible_2_6:w;
    }