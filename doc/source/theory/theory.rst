############
Theory
############

Lethe uses the Finite Element Method (FEM) to solve multiple partial differential equations related to fluid dynamics, the most important one being the incompressible Navier-Stokes equations. The aim of this theory guide is to provide information on the exact equations solved within Lethe and give basic indications on the finite element strategies used (e.g. stabilization) to solve them. This theory guide is not intended to be an introduction to FEM nor the deal.II library. We expect the reader to have a working knowledge of FEM. For more information on FEM we refer the reader to the exquisite `video capsules <https://www.math.colostate.edu/~bangerth/videos.html>`_ created by Wolfgang Bangerth or to books on the topic such as the book by Larson `[1] <https://books.google.ca/books?id=Vek_AAAAQBAJ&lpg=PR3&ots=Ck6m3Q7VxP&dq=The%20finite%20element%20method%3A%20theory%2C%20implementation%2C%20and%20applications&lr&hl=pt-BR&pg=PR3#v=onepage&q=The%20finite%20element%20method:%20theory,%20implementation,%20and%20applications&f=false>`_, which is very accessible. To learn the basics of the deal.II library, we refer the reader to the `deal.II tutorials <https://www.dealii.org/current/doxygen/deal.II/Tutorial.html>`_.

To solve the flow of granular matter, Lethe uses the Discrete Element Method (DEM). This method is clearly not as well established as the FEM. The DEM component of Lethe is, however, fully described in `[2] <arXiv:2106.09576>`_. There are also multiple review articles which give an introduction to this method. Our research group has written a review article on DEM `[3] <https://doi.org/10.1002/cjce.23501>`_ and its application to complex chemical engineering problems `[4] <https://doi.org/10.1016/j.ces.2020.115646>`_. There are multiple other review articles that give a good introduction to DEM such as the one by Bertrand `[5] <https://doi.org/10.1016/j.ces.2004.11.048>`_ or Zhu `[6] <https://doi.org/10.1016/j.ces.2006.12.089>`_. There is also a great book on the topic of DEM and its coupling with CFD by Norouzi et al. `[7] <https://books.google.ca/books?id=7DQWDQAAQBAJ&lpg=PA11&ots=9K7iyGPERX&dq=Coupled%20CFD-DEM%20modeling%3A%20formulation%2C%20implementation%20an&lr&hl=pt-BR&pg=PA11#v=onepage&q=Coupled%20CFD-DEM%20modeling:%20formulation,%20implementation%20an&f=false>`_.

Finally, Lethe also contains an unresolved and a resolved coupling between its CFD and DEM solvers. For an introduction on this topic, we refer the reader to another introductory review article by our group `[8] <https://doi.org/10.1002/cjce.23686>`_.


This theory guide is divided into four core sections:

* Computational Fluid Dynamics (CFD) and other continuum based approaches

* Discrete Element Method for granular flows

* Unresolved CFD-DEM

* Resolved CFD-DEM

Computational Fluid Dynamics
=============================

.. toctree::
    :maxdepth: 3
    :glob:
    :titlesonly:

    fluid_dynamics/fluid_dynamics

Heat transfer
----------------

**Under construction**


Advection-diffusion of a passive tracer
----------------------------------------
**Under construction**


Volume-of-Fluid method for two-phase flows
-------------------------------------------
**Under construction**


Discrete Element Method
=============================
.. toctree::
    :maxdepth: 3
    :glob:
    :titlesonly:

    dem/dem


Unresolved CFD-DEM
=============================
.. toctree::
    :maxdepth: 3
    :glob:
    :titlesonly:

    unresolved_cfd-dem/unresolved_cfd-dem

Resolved CFD-DEM
=============================
**Under construction**


References
-------------
`[1] <https://books.google.ca/books?id=Vek_AAAAQBAJ&lpg=PR3&ots=Ck6m3Q7VxP&dq=The%20finite%20element%20method%3A%20theory%2C%20implementation%2C%20and%20applications&lr&hl=pt-BR&pg=PR3#v=onepage&q=The%20finite%20element%20method:%20theory,%20implementation,%20and%20applications&f=false>`_ Larson, Mats G., and Fredrik Bengzon. The finite element method: theory, implementation, and applications. Vol. 10. Springer Science & Business Media, 2013.

`[2] <arXiv:2106.09576>`_ Golshan, Shahab, et al. Lethe-DEM: An open-source parallel discrete element solver with load balancing. arXiv preprint arXiv:2106.09576 (2021).

`[3] <https://doi.org/10.1002/cjce.23501>`_ Blais, Bruno, et al. Experimental methods in chemical engineering: Discrete element method—DEM. The Canadian Journal of Chemical Engineering 97.7 (2019): 1964-1973.

`[4] <https://doi.org/10.1016/j.ces.2020.115646>`_ Golshan, Shahab, et al. Review and implementation of CFD-DEM applied to chemical process systems. Chemical Engineering Science 221 (2020): 115646.

`[5] <https://doi.org/10.1016/j.ces.2004.11.048>`_ Bertrand, F., L-A. Leclaire, and G. Levecque. DEM-based models for the mixing of granular materials. Chemical Engineering Science 60.8-9 (2005): 2517-2531.

`[6] <https://doi.org/10.1016/j.ces.2006.12.089>`_ Zhu, H. P., et al. Discrete particle simulation of particulate systems: theoretical developments." Chemical Engineering Science 62.13 (2007): 3378-3396.

`[7] <https://books.google.ca/books?id=7DQWDQAAQBAJ&lpg=PA11&ots=9K7iyGPERX&dq=Coupled%20CFD-DEM%20modeling%3A%20formulation%2C%20implementation%20an&lr&hl=pt-BR&pg=PA11#v=onepage&q=Coupled%20CFD-DEM%20modeling:%20formulation,%20implementation%20an&f=false>`_ Norouzi, Hamid Reza, et al. Coupled CFD-DEM modeling: formulation, implementation and application to multiphase flows. John Wiley & Sons, 2016.

`[8] <https://doi.org/10.1002/cjce.23686>`_ Bérard, Ariane, Gregory S. Patience, and Bruno Blais. Experimental methods in chemical engineering: Unresolved CFD‐DEM. The Canadian Journal of Chemical Engineering 98.2 (2020): 424-440.

