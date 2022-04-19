############
Theory
############

Lethe uses the Finite Element Method (FEM) to solve multiple partial differential equations related to fluid dynamics, the most important one being the incompressible Navier-Stokes equations. The aim of this theory guide is to provide information on the exact equations solved within Lethe and give basic indications on the finite element strategies used (e.g. stabilization) to solve them. This theory guide is not intended to be an introduction to FEM nor the deal.II library. We expect the reader to have a working knowledge of FEM. For more information on FEM we refer the reader to the exquisite `video capsules <https://www.math.colostate.edu/~bangerth/videos.html>`_ created by Wolfgang Bangerth or to books on the topic such as the book by Larson [1], which is very accessible. To learn the basics of the deal.II library, we refer the reader to the `deal.II tutorials <https://www.dealii.org/current/doxygen/deal.II/Tutorial.html>`_.

To solve the flow of granular matter, Lethe uses the Discrete Element Method (DEM). This method is clearly not as well established as the FEM. The DEM component of Lethe is, however, fully described in [2]. There are also multiple review articles which give an introduction to this method. Our research group has written a review article on DEM [3] and its application to complex chemical engineering problems [4]. There are multiple other review articles that give a good introduction to DEM such as the one by Bertrand [5] or Zhu [6]. There is also a great book on the topic of DEM and its coupling with CFD by Norouzi et al. [7].

Finally, Lethe also contains an unresolved and a resolved coupling between its CFD and DEM solvers. For an introduction on this topic, we refer the reader to another introductory review article by our group [8].


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
**Under construction**


Unresolved CFD-DEM
=============================
**Under construction**


Resolved CFD-DEM
=============================
**Under construction**


References
-------------
[1] Larson, Mats G., and Fredrik Bengzon. The finite element method: theory, implementation, and applications. Vol. 10. Springer Science & Business Media, 2013.

[2] Golshan, Shahab, et al. "Lethe-DEM: An open-source parallel discrete element solver with load balancing." arXiv preprint arXiv:2106.09576 (2021).

[3] Blais, Bruno, et al. "Experimental methods in chemical engineering: Discrete element method—DEM." The Canadian Journal of Chemical Engineering 97.7 (2019): 1964-1973.

[4] Golshan, Shahab, et al. "Review and implementation of CFD-DEM applied to chemical process systems." Chemical Engineering Science 221 (2020): 115646.

[5] Bertrand, F., L-A. Leclaire, and G. Levecque. "DEM-based models for the mixing of granular materials." Chemical Engineering Science 60.8-9 (2005): 2517-2531.

[6] Zhu, H. P., et al. "Discrete particle simulation of particulate systems: theoretical developments." Chemical Engineering Science 62.13 (2007): 3378-3396.

[7] Norouzi, Hamid Reza, et al. Coupled CFD-DEM modeling: formulation, implementation and application to multiphase flows. John Wiley & Sons, 2016.

[8] Bérard, Ariane, Gregory S. Patience, and Bruno Blais. "Experimental methods in chemical engineering: Unresolved CFD‐DEM." The Canadian Journal of Chemical Engineering 98.2 (2020): 424-440.

