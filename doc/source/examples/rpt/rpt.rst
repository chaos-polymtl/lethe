**************************************
Radioactive Particle Tracking (RPT)
**************************************

.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:

    photon-count-calculation-in-a-cylindrical-vessel/photon-count-calculation-in-a-cylindrical-vessel
    tuning-parameters-with-nomad/tuning-parameters-with-nomad
    particle-fem-reconstruction/particle-fem-reconstruction

.. graphviz:: 

    digraph rpt_diagram {
      graph [bgcolor="transparent", align=true, ranksep=1.5];
      node [fontname=Arial, fontsize=20, shape=box, fontcolor=royalblue, color=royalblue, height=1];
      edge [color=royalblue];
      rankdir="LR";
      size = "9,9";

      rpt [label="Radioactive Particle \nTracking (RPT)", href="https://lethe-cfd.github.io/lethe/documentation/examples/rpt/rpt.html"];

      rpt_1 [label="Photon Count Calculation \nin a Cylindrical Vessel", href="https://lethe-cfd.github.io/lethe/documentation/examples/rpt/photon-count-calculation-in-a-cylindrical-vessel/photon-count-calculation-in-a-cylindrical-vessel.html", tooltip="Photon count calculation in a cylindrical vessel"];

      rpt_2 [label="Tuning Count Calculation Model \nParameters with NOMAD", href="https://lethe-cfd.github.io/lethe/documentation/examples/rpt/tuning-parameters-with-nomad/tuning-parameters-with-nomad.html", tooltip="Tuning count calculation model parameters with NOMAD"];

      rpt_3 [label="Particle FEM reconstruction", href="https://lethe-cfd.github.io/lethe/documentation/examples/rpt/particle-fem-reconstruction/particle-fem-reconstruction.html", tooltip="Particle FEM reconstruction"];

      rpt -> rpt_1:w;
      rpt -> rpt_2:w;  
      rpt -> rpt_3:w;
    }