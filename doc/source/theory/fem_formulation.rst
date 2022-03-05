FEM formulation
####################

This section describes the FEM formulation used within Lethe. First we describe the weak-form, then we introduced the two approaches that are available in Lethe to solve the resulting set of equations.


Starting from the incompressible Navier-Stokes equations:

.. math::
    \partial_j u_j = 0 

    \partial_i+ u_i \partial_i u_j = -\frac{1}{p} \partial_i p + \nu \partial_j \partial_j u_i


