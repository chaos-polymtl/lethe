#####
Lethe
#####

.. image:: ../../logo/lethe-logo-with-bkgd.png
  :alt: Lethe
  :align: center

**Lethe** (pronounced /ˈliːθiː/) is an open-source computational fluid dynamics
(CFD), discrete element method (DEM) and coupled CFD-DEM
software which uses high-order continuous Galerkin formulations to
simulate single and multiphase flows.
Lethe contains a family of solvers that are based on
`deal.II <https://www.dealii.org/>`_, a finite element library.
Through deal.II, Lethe uses `Trilinos <https://trilinos.github.io/>`_, for
its sparse linear algebra routines and `p4est <https://www.p4est.org/>`_
for its distributed adaptative quadtrees and octrees.

Lethe is named after the river of forgetfulness, which, according to
`Wikipedia <https://en.wikipedia.org/wiki/Lethe>`_ :

   is one of the five rivers of the Greek underworld[,] the other four
   [being] Acheron (the river of sorrow), Cocytus (the river of
   lamentation), Phlegethon (the river of fire) and Styx (the river that
   separates Earth and the Underworld).
   …
   The shades of the dead were required to drink the waters of the Lethe
   in order to forget their earthly life.

Lethe is described `here <https://doi.org/10.1016/j.softx.2020.100579>`_.
It originally began as a weekend project, but slowly grew into a robust and efficient CFD, DEM and CFD-DEM
software. Lethe is under active development.

.. note::
   Lethe would not exist without the dedicated work of the deal.II
   authors.
   The authors of Lethe would like to emphasize that without deal.II,
   dedicated solvers like Lethe could not exist.

Contents
--------

.. toctree::
   :maxdepth: 2

   installation/installation
   applications
   structure
   first_simulation
   parameters/parameters
   examples/examples
   theory/theory
   tools/tools
   publications
   referencing
   contributing/contributing

Installation
------------

Follow the instructions in the :doc:`installation/installation` section.

Authors
-------

Main developer:

* `Bruno Blais <https://www.polymtl.ca/expertises/en/blais-bruno>`_

Contributors (in alphabetical order):

* Amishga Alphonius
* Antoine Avrit
* Lucka Barbeau
* Valérie Bibeau
* Bianca Bugeaud
* Audrey Collard-Daigneault
* Carole-Anne Daunais
* Toni El Geitani
* Olivier Gaboriault
* Simon Gauvin
* Shahab Golshan
* Olivier Guévremont
* Marion Hanne
* Jeanne Joachim
* Rajeshwari Kamble
* Martin Kronbichler
* Maxime Landry
* Pierre Laurentin
* Charles Le Pailleur
* Oreste Marquis
* Ghazaleh Mirakhori
* Thomas Monatte
* Peter Munch
* Victor Oliveira Ferreira
* Hélène Papillon-Laroche
* Paul Patience
* Laura Prieto Saavedra
* Catherine Radburn
* Philippe Rivest
* Mikael Vaillant
