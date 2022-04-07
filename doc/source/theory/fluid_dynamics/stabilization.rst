On the need for stabilization
###############################

It is well known that, in the absence of stabilization, the choice of the velocity and pressure finite element spaces must be done to ensure that the `Ladyzenskaya-Babuska-Brezzi (LBB) inf-sup <https://en.wikipedia.org/wiki/Ladyzhenskaya%E2%80%93Babu%C5%A1ka%E2%80%93Brezzi_condition>`_ condition is met. It is also known that the Galerkin approximation of the Navier–Stokes equations may fail in convection dominated flows, for which there are boundary layers where the velocity solution and its gradient exhibit rapid variation over short length scales. In these regions, the classical Galerkin approach may lead to numerical oscillations which may greatly pollute the solution over the entire domain. Stabilization of the elements, which is used in Lethe, can circumvent the limitations of the classical Galerkin approach and alleviate the  need to use LBB stable elements. This notably allows for the use of equal order elements (such as P1-P1).


SUPG/SSPG monolithic formulation
-----------------------------------
**Under construction**


Grad-Div block formulation
------------------------------------
**Under construction**


Galerkin Least-Squares formulation
-----------------------------------
**Under construction**

