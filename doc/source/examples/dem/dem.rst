****************************
Discrete Element Method
****************************

We organize the DEM examples from the simplest (first example) to the most complicated example (last example). Example 1 (packing in circle), shows a two-dimensional simulation of particles packing in a circle. This example is extended in example 2 (packing in ball) to a three-dimensional simulation. Example 3 (rotating drum), shows a simulation in which the particles are packed in a cylinder and the cylinder is rotating). We simulate this rotation using a ``rotational`` boundary condition in Lethe-DEM. In example 4 (rotating drum with post-processing) we explain the post-processing features of Lethe-DEM using the same rotating drum simulation in example 3. In example 5 (rotation of box), we simulate the dynamics of particles in a rotating box. The key difference between examples 3 and 5 is that example 3 uses a ``rotational`` boundary condition to simulate the rotating drum, while example 5 uses a rotation of triangulation (``grid motion``). In example 6 (silo), we simulate the filling (using a ``floating wall``) and discharge of particles in a wedge-shaped silo.

.. toctree::
    :maxdepth: 1
    :glob:

    **/*
