![Lethe](logo/logo_black.png?raw=true)
What is Lethe?
================

Lethe is an Open-source Computational Fluid dynamics Software which uses high-order Continuous Galerkin formulations to solve the incompressible Navier-Stokes equations (among others). Lethe contains a family of solvers that are based on the deal.II open source framework for its finite elements formulation(https://github.com/dealii/dealii). Through deal.II, Lethe uses Trilinos for its sparse linear algebra and p4est for distribute adaptative quad/oct trees.

Lethe is named after the river of forgetfulness, which is one of the five rivers of the Greek underworld, the other four being Styx, Acheron (the river of sorrow), Cocytus (the river of lamentation) and Phlegethon (the river of fire) (https://en.wikipedia.org/wiki/Lethe). The shades of the dead were required to drink the waters of the Lethe in order to forget their earthly life.

Lethe originally started that as a week-end project, but ended up slowly growing into a real research CFD code. It is still an immature project, but is under active development.

Note: Lethe would not exist without the thorough dedicated work of the deal.II. The authors of Lethe would like to emphasize that without deal.II, dedicade solvers like Lethe could not exist.

Authors:
--------
Main developer :
- Bruno Blais (https://www.polymtl.ca/expertises/en/blais-bruno)

Contributors:
- Rajeshwari Kamble
- Antoine Avrit

Installation:
------------
Lethe can be installed on Linux either on a dedicated machine or on a virtual machine.
Lethe has been tested on numerous distributions, including Ubuntu 18.04 LTS and Centos 7

The installation of Lethe consists in three major steps:
* Installation of deal.II dependency (p4est and trilinos)
* Installation of deal.II
* Installation of lethe

Lethe cannot be installed without p4est and Trilinos. Although Lethe can be run in serial and distributed (parallel) mode, it uses p4est and Trilinos for mesh decomposition and linear solvers respectively. Lethe generally supports the lastest release of the dealii library.

### 1. Installation of dealii dependencies
#### P4est
To install p4est, the usual installation procedure of dealii can be followed:
https://www.dealii.org/current/external-libs/p4est.html

#### Trilinos:
The installation of Trilinos should be done using the installation procedure of dealii:
https://www.dealii.org/current/external-libs/trilinos.html

### 2. Installation of dealii
* Clone dealii from the official repository (https://github.com/dealii/dealii)
`$ git clone https://github.com/dealii/dealii `
* Configure dealii in a build folder at the same level as the source code
`$ mkdir build`
`$ cd build`

Depending on how you have installed p4est and Trilinos, you may need to specify the installation folder of the two libraries

`$ cmake ../dealii -DDEAL_II_WITH_TRILINOS=ON -DTRILINOS_DIR=path/to/your/trilinos/installation  
    -DDEAL_II_WITH_P4EST=ON -DP4EST_DIR=path/to/your/p4est/installation -DCMAKE_INSTALL_PREFIX=/path/to/desired/installation`

* Compile dealii
`$ make -j install`

* Create an environment variable for the DEALII directory.
`$ export DEAL_II_DIR=/path/to/dealii/installation`

* It is generally recommended to add a variable to your bashrc so it is always loaded:
`$ echo "export DEAL_II_DIR=/path/to/dealii/installation" >> ~/.bashrc`


### 3. Installation of lethe

* Clone lethe from the official repository (https://github.com/dealii/dealii)
* Create a build folder at the same level as the lethe folder
`mkdir build`

* Compile Lethe choosing the compilation option (Debug or Release)

`cmake ../lethe -DCMAKE_BUILD_TYPE=Debug`
or
`cmake ../lethe -DCMAKE_BUILD_TYPE=Release`


License:
--------
Please see the file [./LICENSE](LICENSE) for details
