[![Build Status](https://travis-ci.com/lethe-cfd/lethe.svg?branch=master)](https://travis-ci.com/lethe-cfd/lethe)
![Lethe](logo/logo_black.png?raw=true)
# What is Lethe?


Lethe is an Open-source Computational Fluid Dynamics (CFD) Software which uses high-order Continuous Galerkin formulations to solve the incompressible Navier-Stokes equations (among others). Lethe contains a family of solvers that are based on the deal.II open source framework for its finite elements formulation (https://github.com/dealii/dealii). Through deal.II, Lethe uses Trilinos for its sparse linear algebra and p4est for distribute adaptative quad/oct trees.

Lethe is named after the river of forgetfulness, which is one of the five rivers of the Greek underworld, the other four being Styx, Acheron (the river of sorrow), Cocytus (the river of lamentation) and Phlegethon (the river of fire) (https://en.wikipedia.org/wiki/Lethe). The shades of the dead were required to drink the waters of the Lethe in order to forget their earthly life.

Lethe originally started that as a week-end project, but ended up slowly growing into a real research CFD code. It is still an immature project, but is under active development.

Lethe is described in the following articles:
- https://doi.org/10.1016/j.softx.2020.100579

Lethe is featured in the following articles:


Note: Lethe would not exist without the thorough dedicated work of the deal.II authors. The authors of Lethe would like to emphasize that without deal.II, dedicated solvers like Lethe could not exist.

Authors:
--------
Main developer :
- Bruno Blais (https://www.polymtl.ca/expertises/en/blais-bruno)

Contributors (in alphabetical order):
- Antoine Avrit
- Lucka Barbeau
- Valérie Bibeau
- Audrey Collard-Daigneault
- Carole-Anne Daunais
- Toni El Geitani
- Simon Gauvin
- Shahab Golshan
- Rajeshwari Kamble
- Ghazaleh Mirakhori

Installation:
------------
Follow the instructions in the [wiki](https://github.com/lethe-cfd/lethe/wiki/Installation).

License:
--------
Please see the file [./LICENSE](LICENSE) for details
