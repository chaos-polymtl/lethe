****************************
Discrete Element Method
****************************

We organize the DEM examples from the simplest to the most complicated example: 


.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:
    :numbered:

    bouncing-particle/bouncing-particle
    packing-in-circle/packing-in-circle
    packing-in-ball/packing-in-ball
    small-scale-rotating-drum/small-scale-rotating-drum
    small-scale-rotating-drum-post-processing/small-scale-rotating-drum-post-processing
    rotating-drum/rotating-drum
    rotating-drum-with-post-processing/rotating-drum-with-post-processing
    rotation-of-box/rotation-of-box
    silo/silo
    rectangular-hopper/rectangular-hopper
    granular-dam-break/granular-dam-break
    granular-mixer/granular-mixer

.. graphviz:: 

    digraph dem_diagram {
      graph [align=true];
      node [fontname=Arial, fontsize=18, shape=box];
      center=true;
      fontname="Comic Sans MS";
      fontsize = "20";
      rankdir="LR";
      size = "9,9";
      
      dem [label="Discrete Element \nMethod", href="https://lethe-cfd.github.io/lethe/examples/dem/dem.html"];

      dem_1 [label="1. Bouncing particle", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/bouncing-particle/bouncing-particle.html"]; 

      dem_2 [label="2. Packing in circle", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/packing-in-circle/packing-in-circle.html"]; 

      dem_3 [label="3. Packing in ball", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/packing-in-ball/packing-in-ball.html"]; 

      dem_4 [label="4. Small scale rotating drum", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/small-scale-rotating-drum/small-scale-rotating-drum.html"]; 
    
      dem_5 [label="5. Small scale rotating drum post-processing", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/small-scale-rotating-drum-post-processing/small-scale-rotating-drum-post-processing.html"]; 

      dem_6 [label="6. Rotating drum", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/rotating-drum/rotating-drum.html"]; 

      dem_7 [label="7. Rotating drum with post-processing", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/rotating-drum-with-post-processing/rotating-drum-with-post-processing.html"]; 

      dem_8 [label="8. Rotation of box", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/rotation-of-box/rotation-of-box.html"]; 

      dem_9 [label="9. Silo", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/silo/silo.html"]; 

      dem_10 [label="10. Rectangular hopper", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/rectangular-hopper/rectangular-hopper.html"]; 

      dem_11 [label="11. Granular dam-break", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/granular-dam-break/granular-dam-break.html"]; 

      dem_12 [label="12. Granular mixer", color=black, href="https://lethe-cfd.github.io/lethe/examples/dem/granular-mixer/granular-mixer.html"]; 

      dem -> dem_1:w;
      dem -> dem_2:w;
      dem -> dem_3:w;
      dem -> dem_4:w;
      dem -> dem_5:w;
      dem -> dem_6:w;
      dem -> dem_7:w;
      dem -> dem_8:w;
      dem -> dem_9:w;
      dem -> dem_10:w;
      dem -> dem_11:w;
      dem -> dem_12:w;
    }