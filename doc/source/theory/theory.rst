############
Theory
############

Lethe is a Computational Fluid Dynamics (CFD) simulation software that uses the Finite Element Method (FEM) to solve the multiple partial differential equations related to fluid dynamics, the most important one being the incompressible Navier-Stokes equations. The aim of this theory guide is to provide information on the exact equations solved within Lethe and give basic indications on the finite element strategies used (e.g. stabilization) to solve them. This theory guide is not intended to be an introduction to FEM nor the deal.II library. We expect the reader to have a working knowledge of FEM. For more information on FEM we refer the reader to the exquisite `video capsules <https://www.math.colostate.edu/~bangerth/videos.html>`_ created by Wolfgang Bangerth or to books on the topic such as the book by Larson `[1] <https://books.google.ca/books?id=Vek_AAAAQBAJ&lpg=PR3&ots=Ck6m3Q7VxP&dq=The%20finite%20element%20method%3A%20theory%2C%20implementation%2C%20and%20applications&lr&hl=pt-BR&pg=PR3#v=onepage&q=The%20finite%20element%20method:%20theory,%20implementation,%20and%20applications&f=false>`_, which is very accessible. To learn the basics of the deal.II library, we refer the reader to the `deal.II tutorials <https://www.dealii.org/current/doxygen/deal.II/Tutorial.html>`_.

To solve the flow of granular matter, Lethe uses the Discrete Element Method (DEM). This method is clearly not as well established as the FEM. The DEM component of Lethe is, however, fully described in `[2] <https://doi.org/10.1007/s40571-022-00478-6>`_. There are also multiple review articles which give an introduction to this method. Our research group has written a review article on DEM `[3] <https://doi.org/10.1002/cjce.23501>`_ and its application to complex chemical engineering problems `[4] <https://doi.org/10.1016/j.ces.2020.115646>`_. There are multiple other review articles that give a good introduction to DEM such as the one by Bertrand `[5] <https://doi.org/10.1016/j.ces.2004.11.048>`_ or Zhu `[6] <https://doi.org/10.1016/j.ces.2006.12.089>`_. There is also a great book on the topic of DEM and its coupling with CFD by Norouzi et al. `[7] <https://books.google.ca/books?id=7DQWDQAAQBAJ&lpg=PA11&ots=9K7iyGPERX&dq=Coupled%20CFD-DEM%20modeling%3A%20formulation%2C%20implementation%20an&lr&hl=pt-BR&pg=PA11#v=onepage&q=Coupled%20CFD-DEM%20modeling:%20formulation,%20implementation%20an&f=false>`_.

Finally, Lethe also contains an unresolved and a resolved coupling between its CFD and DEM solvers. For an introduction on this topic, we refer the reader to another introductory review article by our group `[8] <https://doi.org/10.1002/cjce.23686>`_.


This theory guide is divided into two global sections:

* Multiphysics (physics supported in Lethe)

* Multiphase (methods employed for multiphase numerical representation)



Multiphysics
=============================

.. toctree::
    :maxdepth: 3
    :glob:
    :titlesonly:

    multiphysics/multiphysics

Multiphase
=============================

.. toctree::
    :maxdepth: 3
    :glob:
    :titlesonly:

    multiphase/multiphase



References
-------------
`[1] <https://books.google.ca/books?id=Vek_AAAAQBAJ&lpg=PR3&ots=Ck6m3Q7VxP&dq=The%20finite%20element%20method%3A%20theory%2C%20implementation%2C%20and%20applications&lr&hl=pt-BR&pg=PR3#v=onepage&q=The%20finite%20element%20method:%20theory,%20implementation,%20and%20applications&f=false>`_ M. G. Larson and F. Bengzon, *The Finite Element Method: Theory, Implementation, and Applications*. Springer Science & Business Media, 2013.

`[2] <https://doi.org/10.1007/s40571-022-00478-6>`_ S. Golshan, P. Munch, R. Gassmöller, M. Kronbichler, and B. Blais, “Lethe-DEM: an open-source parallel discrete element solver with load balancing,” *Comput. Part. Mech.*, vol. 10, no. 1, pp. 77–96, Feb. 2023, doi: 10.1007/s40571-022-00478-6.

`[3] <https://doi.org/10.1002/cjce.23501>`_ B. Blais, D. Vidal, F. Bertrand, G. S. Patience, and J. Chaouki, “Experimental Methods in Chemical Engineering: Discrete Element Method—DEM,” *Can. J. Chem. Eng.*, vol. 97, no. 7, pp. 1964–1973, 2019, doi: 10.1002/cjce.23501.

`[4] <https://doi.org/10.1016/j.ces.2020.115646>`_ S. Golshan, R. Sotudeh-Gharebagh, R. Zarghami, N. Mostoufi, B. Blais, and J. A. M. Kuipers, “Review and implementation of CFD-DEM applied to chemical process systems,” *Chem. Eng. Sci.*, vol. 221, p. 115646, Aug. 2020, doi: 10.1016/j.ces.2020.115646.

`[5] <https://doi.org/10.1016/j.ces.2004.11.048>`_ F. Bertrand, L.-A. Leclaire, and G. Levecque, “DEM-based models for the mixing of granular materials,” *Chem. Eng. Sci.*, vol. 60, no. 8, pp. 2517–2531, Apr. 2005, doi: 10.1016/j.ces.2004.11.048.

`[6] <https://doi.org/10.1016/j.ces.2006.12.089>`_ H. P. Zhu, Z. Y. Zhou, R. Y. Yang, and A. B. Yu, “Discrete particle simulation of particulate systems: Theoretical developments,” *Chem. Eng. Sci.*, vol. 62, no. 13, pp. 3378–3396, Jul. 2007, doi: 10.1016/j.ces.2006.12.089.

`[7] <https://books.google.ca/books?id=7DQWDQAAQBAJ&lpg=PA11&ots=9K7iyGPERX&dq=Coupled%20CFD-DEM%20modeling%3A%20formulation%2C%20implementation%20an&lr&hl=pt-BR&pg=PA11#v=onepage&q=Coupled%20CFD-DEM%20modeling:%20formulation,%20implementation%20an&f=false>`_ H. R. Norouzi, R. Zarghami, R. Sotudeh-Gharebagh, and N. Mostoufi, *Coupled CFD-DEM Modeling: Formulation, Implementation and Application to Multiphase Flows*. John Wiley & Sons, 2016.

`[8] <https://doi.org/10.1002/cjce.23686>`_ A. Bérard, G. S. Patience, and B. Blais, “Experimental methods in chemical engineering: Unresolved CFD-DEM,” *Can. J. Chem. Eng.*, vol. 98, no. 2, pp. 424–440, 2020, doi: 10.1002/cjce.23686.
