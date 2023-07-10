########
Examples
########

.. toctree::
    :maxdepth: 2
    :glob:
    :titlesonly:

    incompressible-flow/incompressible-flow
    multiphysics/multiphysics
    dem/dem
    unresolved-cfd-dem/unresolved-cfd-dem
    sharp-immersed-boundary-solver/sharp-immersed-boundary-solver
    rpt/rpt

.. graphviz:: 

    digraph incompressible_diagram {
      graph [align=true];
      node [fontname=Arial, fontsize=18];
      center=true;
      fontname="Comic Sans MS";
      fontsize = "20";
      rankdir="LR";
      size = "8,8";
      incompressible_flow [label="Incompressible \nFlow", href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/incompressible-flow.html"];
      incompressible_1 [label="Lid-driven cavity flow",shape=polygon,sides=4, color=blue, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-lid%E2%80%90driven-cavity-flow/lid%E2%80%90driven-cavity-flow.html"]; 
      incompressible_2 [label="Flow around a cylinder", shape=polygon,sides=4, color=blue, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-flow-around-cylinder/2d-flow-around-cylinder.html"];
      incompressible_3 [label="Taylor-Couette flow", shape=polygon,sides=4, color=blue, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-taylor-couette-flow/2d-taylor-couette-flow.html"];
      incompressible_4 [label="Flow past a backward-facing step", shape=polygon,sides=4, color=blue, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-backward-facing-step/2d-backward-facing-step.html"];
      incompressible_5 [label="Transient flow around an Ahmed body", shape=polygon,sides=4, color=blue, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-transient-around-ahmed-body/2d-transient-around-ahmed-body.html"];
      incompressible_6 [label="Transient flow around a cylinder", shape=polygon,sides=4, color=blue, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-transient-flow-around-cylinder/2d-transient-flow-around-cylinder.html"];
      incompressible_7 [label="Flow around NACA0012 at low Reynolds number", shape=polygon,sides=4, color=blue, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-naca0012-low-reynolds/2d-naca0012-low-reynolds.html"];
      incompressible_8 [label="Flow around a sphere", shape=polygon,sides=4, color=red, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-flow-around-sphere/flow-around-sphere.html"];
      incompressible_9 [label="Flow over periodic hills", shape=polygon,sides=4, color=red, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-flow-over-periodic-hills/3d-flow-over-periodic-hills.html"];
      incompressible_10 [label="Taylor-Couette flow using Nitsche immersed boundary", shape=polygon,sides=4, color=blue, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-taylor-couette-flow-nitsche/2d-taylor-couette-flow-nitsche.html"];
      incompressible_11 [label="Ribbon mixer using a single rotating reference frame", shape=polygon,sides=4, color=red, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-mixer-using-single-rotating-frame/3d-mixer-using-single-rotating-frame.html"];
      incompressible_12 [label="Mixer with pitched-blade turbine impeller using Nitsche immersed boundary", shape=polygon,sides=4, color=red, href="https://lethe-cfd.github.io/lethe/examples/incompressible-flow/3d-nitsche-mixer-with-pbt-impeller/nitsche-mixer-with-pbt-impeller.html"];

      incompressible_flow -> incompressible_1;
      incompressible_flow -> incompressible_2;
      incompressible_flow -> incompressible_3;
      incompressible_flow -> incompressible_4;
      incompressible_flow -> incompressible_5;
      incompressible_flow -> incompressible_6;
      incompressible_flow -> incompressible_7;
      incompressible_flow -> incompressible_8;
      incompressible_flow -> incompressible_9;
      incompressible_flow -> incompressible_10;
      incompressible_flow -> incompressible_11;
      incompressible_flow -> incompressible_12;

        subgraph cluster_key {
            two_dimensions [label="2D",shape=polygon,sides=4, color=blue]; 
            three_dimensions [label="3D",shape=polygon,sides=4, color=red]; 
        }
    }

.. graphviz:: 

    digraph multiphysics_diagram {
      graph [align=true];
      node [fontname=Arial, fontsize=18];
      center=true;
      fontname="Comic Sans MS";
      fontsize = "20";
      rankdir="LR";
      
      multiphysics [label="Multiphysics", href="https://lethe-cfd.github.io/lethe/examples/multiphysics/multiphysics.html"];
      multiphysics -> multiphysics_1;
      multiphysics -> multiphysics_2;
      multiphysics -> multiphysics_3;
    }

.. graphviz::

    digraph second_diagram {
      align = "right";
      rankdir="LR";
      multiphysics [label="Multiphysics", href="https://lethe-cfd.github.io/lethe/examples/multiphysics/multiphysics.html"];
      multiphysics -> multiphysics_1;
      multiphysics -> multiphysics_2;
      multiphysics -> multiphysics_3;

      dem [label="Discrete Element \nMethod", href="https://lethe-cfd.github.io/lethe/examples/dem/dem.html"];
      dem -> dem_1;
      dem -> dem_2;
      dem -> dem_3;
      dem -> dem_4;
      dem -> dem_5;
      dem -> dem_6;
      dem -> dem_7;
      dem -> dem_8;
      dem -> dem_9;
      dem -> dem_10;
      dem -> dem_11;
      dem -> dem_12;

      unresolved_cfd_dem [label="Unresolved \nCFD-DEM", href="https://lethe-cfd.github.io/lethe/examples/unresolved-cfd-dem/unresolved-cfd-dem.html"];
      unresolved_cfd_dem -> cfd_dem_1;
      unresolved_cfd_dem -> cfd_dem_2;
      unresolved_cfd_dem -> cfd_dem_3;
      unresolved_cfd_dem -> cfd_dem_4;

      sharp_immersed_boundary_solver [label="Sharp Immersed \nBoundary Solver", href="https://lethe-cfd.github.io/lethe/examples/sharp-immersed-boundary-solver/sharp-immersed-boundary-solver.html"];
      sharp_immersed_boundary_solver -> sharp_1;
      sharp_immersed_boundary_solver -> sharp_2;

      rpt [label="Radioactive Particle Tracking (RPT)", href="https://lethe-cfd.github.io/lethe/examples/rpt/rpt.html"];
      rpt -> rpt_1;
      rpt -> rpt_2;      
    }