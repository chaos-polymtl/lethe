**************************************
Radioactive Particle Tracking (RPT)
**************************************

.. toctree::
    :maxdepth: 1
    :glob:
    :numbered:

    photon-count-calculation-in-a-cylindrical-vessel/photon-count-calculation-in-a-cylindrical-vessel
    tuning-parameters-with-nomad/tuning-parameters-with-nomad

.. graphviz:: 

    digraph rpt_diagram {
      graph [align=true];
      node [fontname=Arial, fontsize=18];
      center=true;
      fontname="Comic Sans MS";
      fontsize = "20";
      rankdir="LR";
      size = "9,9";
      
      rpt [label="Radioactive Particle \nTracking (RPT)", href="https://lethe-cfd.github.io/lethe/examples/rpt/rpt.html"];

      rpt_1 [label="1. Photon count calculation \nin a cylindrical vessel",shape=polygon,sides=4, color=black, href="https://lethe-cfd.github.io/lethe/examples/rpt/photon-count-calculation-in-a-cylindrical-vessel/photon-count-calculation-in-a-cylindrical-vessel.html"]; 

      rpt_2 [label="2. Tuning Count Calculation Model \nParameters with NOMAD",shape=polygon,sides=4, color=black, href="https://lethe-cfd.github.io/lethe/examples/rpt/tuning-parameters-with-nomad/tuning-parameters-with-nomad.html"]; 

      rpt -> rpt_1;
      rpt -> rpt_2;  
    }