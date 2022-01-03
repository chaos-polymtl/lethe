############
Installation
############

.. toctree::
    :maxdepth: 2
    :glob:
    :titlesonly:

    docker


Lethe can be installed on Linux, Mac OS and under the Windows Subsystem for Linux (WSL). Lethe has been tested on numerous distributions, including Ubuntu 18.04 LTS, Ubuntu 20.04 LTS, Centos 7 and Manjaro. Lethe requires a modern version of the deal.II library. At the time of this writing, deal.II 9.3 and deal.II 10.0pre (the master branch version) are supported. The compatibility with these two branches is ensured by Continuous Integration (CI) using Github Actions. A `dealii fork <https://github.com/lethe-cfd/dealii>`_ is maintained by the Lethe organization. This fork does not include any modification to the library, but is the latest deal.II version against which Lethe was tested. We work hard to ensure compatibility with the latest deal.II version and we do not modify the library except through pull request on the official deal.II repository.

The installation of Lethe consists in three major steps:
1. Installation of dealii dependency (mpi, numdiff, p4est, trilinos and METIS)
2. Installation of the dealii library
3. Installation of lethe

There are two methods to install the dependencies. They can be either installed manually or through ``candi``. The candi.sh shell script can download, configure, build, and install all the dependencies along with dealii. This is the easiest way, since it requires much less manual intervention. It is more machine time consuming since all common dependencies of deal.II will be installed (On a single core computer with 8Gb of ram, count up to 8 hours for a first installation, and 3 hours for a dealii update. On a machine with 16 cores and 32 Gb or ram, this process can take an hour or so). This is the recommended way to install the dependencies.

Lethe cannot be installed if deal.II has not been configured with p4est, Trilinos and METIS. Although Lethe can be run in serial and distributed (parallel) mode, it uses p4est and Trilinos for mesh decomposition and the linear solvers respectively. 

Installing deal.II using candi (Step #1)
-----------------------------------------

To install the dependencies (mpi, numdiff, p4est, trilinos and METIS) all together by Candi, the procedure on the candi repository can be followed:

https://github.com/dealii/candi.git

Clone the candi git reposiroty in a folder of your choice  (e.g. `/home/username/Software`). You can edit candi.cfg if you want to alter which dependencies are compiled. This can notably be used to force the installation of the deal.II master version instead of the current stable version by setting the `STABLE_BUILD=false`

From the candi folder, the installation of candi can be launched using:
`./candi.sh -j $num_proc --prefix=$path`
where `$num_proc` is the number of threads you want to use to compile deal.II and `$path$` the installation prefix that is desired (e.g. `/home/username/Software/`). For a computer with 8Gb of RAM, 1 thread should be used. For 16 Gb, 4 threads is reasonable. For 32 Gb, 8 threads and for 64 Gb 16 threads or more can be used.


After installation, you should have a `deal.ii-candi` folder the installation prefix directory, with the dealii folder of the desired version (see section [Updating Lethe - with candi](https://github.com/lethe-cfd/lethe/wiki/Installation#with-candi)), as well as the required dependencies (p4est, trilinos, etc.).

After installation, add an environment variable to your `.bashrc`:
`$ echo "export DEAL_II_DIR=/home/<usr>/deal.ii-candi/deal.II-<version>" >> ~/.bashrc`

Installing Deal.II manually (Step #1)
--------------------------------------

Installation of dealii dependencies (Step #1)

MPI
~~~~~

MPI, for Message Passing Interface, is required for the installation of all components of Lethe. Therefore it should be the first thing you install. On most Linux based OS, MPI can be installed directly from your package manager. As an example, in Ubuntu:

`sudo apt-get install libopenmpi-dev openmpi-bin openmpi-common`

numdiff
~~~~~~~~

numdiff is used within the automatic testing procedure of Lethe to compare files obtained through floating point arithmetic. Without numdiff, Lethe automatic tests may fail when they should not. Numdiff can be installed directly from your package manager.

`sudo apt-get install numdiff`

P4est
~~~~~~~

To install p4est, the usual installation procedure of dealii can be followed:

https://www.dealii.org/current/external-libs/p4est.html

Trilinos
~~~~~~~~~

The installation of Trilinos should be done using the installation procedure of dealii:

https://www.dealii.org/current/external-libs/trilinos.html

METIS
~~~~~~~

METIS is used for mesh partitioning for parallel computing purposes, specifically in cases with simplex grids. Some parallel simulations are impossible without it. It can be downloaded here or through candi.

http://glaros.dtc.umn.edu/gkhome/metis/metis/download

Installation of dealii
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Clone deal.II from the official repository (https://github.com/dealii/dealii)

`$ git clone https://github.com/dealii/dealii `

* Configure deal.II in a build folder at the same level as the source code

`$ mkdir build`

`$ cd build`

Depending on how you have installed p4est, Trilinos and METIS, you may need to specify the installation folder of the three libraries

`$ cmake ../dealii -DDEAL_II_WITH_MPI=ON -DDEAL_II_WITH_TRILINOS=ON -DTRILINOS_DIR=path/to/your/trilinos/installation  
    -DDEAL_II_WITH_P4EST=ON -DP4EST_DIR=path/to/your/p4est/installation  -DDEAL_II_WITH_METIS=ON -DMETIS_DIR=path/to/your/metis/installation -DCMAKE_INSTALL_PREFIX=/path/to/desired/installation`

* Compile dealii

`$ make -j<nprocessor> install`

* Create an environment variable for the DEALII directory. 

`$ export DEAL_II_DIR=/path/to/dealii/installation`

* It is generally recommended to add variable to your bashrc so it is always loaded:

`$ echo "export DEAL_II_DIR=/path/to/dealii/installation" >> ~/.bashrc`

Installation of lethe (Step #2)
-------------------------------

* Clone lethe from the official repository (https://github.com/lethe-cfd/lethe)

* Create a build folder at the same level as the lethe folder

`mkdir build`

* Compile Lethe choosing the compilation option (Debug or Release)

`cmake ../lethe -DCMAKE_BUILD_TYPE=Debug`

or

`cmake ../lethe -DCMAKE_BUILD_TYPE=Release`

Compiling

`make -j<numprocs>`

# Testing your installation

The validity of your installation can be tested by running Lethe test suite. Within the build folder, the test suite can be launched with the following command:

`ctest -j<numprocs>`


Updating deal.II
-------------------

Through the git repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The deal.II version supported by Lethe is updated and tested every week, see the repository here: [dealii fork](https://github.com/lethe-cfd/dealii). If Lethe was installed with this forked version of deal.II, updating your deal.II installation is as simple as pulling the repository and recompiling the deal.ii library. If your deal.II was installed manually using the deal.II master repository, the same process can be used.

With candi
~~~~~~~~~~~~~
Tested on Linux Ubuntu 18.04.5 LTS
1. In the candi folder (for instance, `/home/<usr>/Softwares/candi`), modify the **candi.cfg** to get the latest dealii version, by changing the `DEAL_II_VERSION` variable in the case of an official release, or by changing the `STABLE_BUILD` in the case of a development release. When this tutorial is written, the official dealii version is 9.2, but Lethe uses features of the upcoming 9.3 under development, so the candi.cfg should show:
<pre><code># Install the following deal.II version:
DEAL_II_VERSION=v9.2.0

# Would you like to build stable version of deal.II?
# If STABLE_BUILD=false, then the development version of deal.II will be
# installed.
#STABLE_BUILD=true
STABLE_BUILD=false</pre></code>

2. Run the command `./candi.sh` to install the new version of dealii (takes about 3h on a laptop workstation)
3. In your `/home/deal.ii-candi` folder, you should have a new folder with the dealii updated version (specified in `DEAL_II_VERSION`, or `deal.II-master` in the case of development version)
4. Edit your `.bashrc` to change the environment variable: 
`$ echo "export DEAL_II_DIR=/home/<usr>/deal.ii-candi/deal.II-master" >> ~/.bashrc`

You could have to restart your computer and to remove your lethe build folder, then recreate it:
<pre><code>mkdir build
cmake ../lethe/ -DCMAKE_BUILD_TYPE=Release <or =Debug>
make -j<numprocs>
</code></pre>

You can now run the ctest (see [Testing your installation](https://github.com/lethe-cfd/lethe/wiki/Installation#testing-your-installation)) to check your new installation.
>>>>>>> 285fe47 (First draft of installation instructions):doc/source/installation.rst
