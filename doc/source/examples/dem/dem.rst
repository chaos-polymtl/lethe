****************************
Discrete Element Method
****************************

We organize the DEM examples from the simplest (first example) to the most complicated example (last example). `Example 1 <./packing-in-circle/packing-in-circle.html>`_ (packing in circle), shows a two-dimensional simulation of particles packing in a circle. This example is extended in `Example 2 <./packing-in-ball/packing-in-ball.html>`_ (packing in ball) to a three-dimensional simulation. `Example 3 <./rotating-drum/rotating-drum.html>`_ (rotating drum), shows a simulation in which the particles are packed in a cylinder and the cylinder is rotating). We simulate this rotation using a ``rotational`` boundary condition in Lethe-DEM. In `Example 4 <./rotating-drum-with-post-processing/rotating-drum-with-post-processing.html>`_ (rotating drum with post-processing) we explain the post-processing features of Lethe-DEM using the same rotating drum simulation in `Example 3 <./rotating-drum/rotating-drum.html>`_. In `Example 5 <./rotation-of-box/rotation-of-box.html>`_ (rotation of box), we simulate the dynamics of particles in a rotating box. The key difference between `Example 3 <./rotating-drum/rotating-drum.html>`_ and `Example 5 <./rotation-of-box/rotation-of-box.html>`_ is that `Example 3 <./rotating-drum/rotating-drum.html>`_ uses a ``rotational`` boundary condition to simulate the rotating drum, while `Example 5 <./rotation-of-box/rotation-of-box.html>`_ uses a rotation of triangulation (``grid motion``). In `Example 6 <./silo/silo.html>`_ (silo), we simulate the filling (using a ``floating wall``) and discharge of particles in a wedge-shaped silo.

.. toctree::
    :maxdepth: 1
    :glob:

    **/*
