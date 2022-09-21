##############################
Regular installation on Linux
##############################

.. figure:: ./images/linux.png
   :height: 100px

.. important::
	Distributions compatibility: Ubuntu 18.04 LTS, Ubuntu 20.04 LTS, Ubuntu 22.04 LTS, Centos 7 and Manjaro

* Dependencies: deal.II library (`deal.II website <https://www.dealii.org/>`_) and its dependencies (MPI, numdiff, p4est, trilinos and METIS)
	* Lethe requires a modern version of the deal.II library. At the time of this writing, ``deal.II 9.4`` and ``deal.II 10.0pre`` (the ``master`` branch version) are supported. 
	* The compatibility with these two branches is ensured by Continuous Integration (CI) using Github Actions. 
	* A `dealii fork <https://github.com/lethe-cfd/dealii>`_ is maintained by Lethe team. This fork does not include any modification to deal.II library, but it is the latest version with which Lethe was tested. We work hard to ensure compatibility with the latest deal.II version and we do not modify the library except through pull requests on the official deal.II repository.

.. warning:: 
	Lethe cannot be installed if deal.II has not been configured with p4est, trilinos and METIS. Although Lethe can be run in serial and parallel mode (through MPI), it uses p4est, METIS and trilinos for mesh decomposition and the linear solvers respectively. 

* Lethe installation steps:
	1. Installation of deal.II dependencies (MPI, numdiff, p4est, trilinos and METIS)
	2. Installation of the deal.II library
	3. Installation of Lethe

* Methods to install deal.II and its dependencies:
	1. through candi shell script (`candi github page <https://github.com/dealii/candi>`_): **this is by far the easiest way to proceed**, since it requires much less manual intervention. Even if you do not want to use candi to install deal.II, you can use it for the dependencies.	
	
	2. manually 

.. warning::
	On a single core computer with 8GB of RAM, count up to 8 hours for a first installation, and 3 hours for a deal.II update. On a machine with 16 cores and 32GB of RAM, this process will take less than an hour or so. The installation of deal.II and its dependencies (especially trilinos), can be extremely RAM consuming. Installation on a machine with less than 8GB of RAM is difficult at best, impossible at worst.


Installing deal.II using candi (Step #1)
-----------------------------------------

To install the dependencies (MPI, p4est, trilinos and METIS) all together using candi, the `procedure <https://github.com/dealii/candi.git>`_ on the candi repository can be followed.

Clone the candi git repository in a folder of your choice  (e.g. ``/home/username/software``). You can edit the ``candi.cfg`` file if you want to alter which dependencies are compiled. This can notably be used to force the installation of the deal.II master version instead of the current stable version by setting the ``STABLE_BUILD=false``

From the candi folder, the installation of candi can be launched using:

.. code-block:: text
  :class: copy-button

  ./candi.sh -j $num_proc --prefix=$path


where ``$num_proc`` is the number of threads you want to use to compile deal.II and ``$path`` the installation prefix that is desired (e.g. ``/home/username/software/candi``). 

.. warning:: 
  For a computer with 8Gb of RAM, 1 thread (``num_proc=1``) should be used. For 16 Gb, 4 threads is reasonable. For 32 Gb, 16 threads or more can be used.


After installation, you should have a ``deal.II-candi`` folder in the installation prefix directory, with the dealii folder of the desired version (see section :ref:`update-dealii`), as well as the required dependencies (p4est, trilinos, etc.).

After installation, add an environment variable to your ``.bashrc`` either manually or through the following command:

.. code-block:: text
  :class: copy-button

  echo "export DEAL_II_DIR=/home/username/deal.ii-candi/deal.II-<version>" >> ~/.bashrc

Installing deal.II manually (Step #1)
--------------------------------------
.. note:: 
  If you have installed deal.II through candi, you can skip right away to :ref:`install-lethe`

We first need to install the deal.II dependencies.


MPI
~~~~~

MPI, for Message Passing Interface, is required for the installation of all components of Lethe. Therefore it should be the first thing you install. On Debian Linux-based OS, MPI can be installed directly from your package manager. 

As an example, in Ubuntu:

.. code-block:: text
  :class: copy-button

  sudo apt-get install libopenmpi-dev openmpi-bin openmpi-common

In Manjaro or other arch based linux distribution:

.. code-block:: text
  :class: copy-button

  sudo pacman -Sy openmpi


numdiff
~~~~~~~~

numdiff is used within the automatic testing procedure of Lethe to compare files obtained through floating point arithmetic. Without numdiff, Lethe automatic tests may fail when they should not. numdiff can be installed directly from your package manager.

.. code-block:: text
  :class: copy-button

  sudo apt-get install numdiff

Regrettably, numdiff is not available in the pacman package manager. It can be downloaded from the following `website <http://www.nongnu.org/numdiff/>`_. If you are using an arch distribution, we assume that you will already know how to carry on with the installation of numdiff. 

P4est
~~~~~~~

To install p4est, the usual `installation <https://www.dealii.org/current/external-libs/p4est.html>`_ of deal.II can be followed.



Trilinos
~~~~~~~~~

The installation of trilinos should be done using the `installation procedure <https://www.dealii.org/current/external-libs/trilinos.html>`_ of deal.II.



METIS
~~~~~~~

METIS is used for mesh partitioning for parallel computing purposes, specifically in cases with simplex grids. It can be downloaded in this `link <http://glaros.dtc.umn.edu/gkhome/metis/metis/download>`_ or through candi.



Installation of deal.II
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clone deal.II from the `deal.ii official repository <https://github.com/dealii/dealii>`_

.. code-block:: text
  :class: copy-button

  git clone https://github.com/dealii/dealii 

Configure deal.II in a build folder at the same level as the source code

.. code-block:: text
  :class: copy-button

  mkdir build
  cd build

Depending on how you have installed p4est, Trilinos and METIS, you may need to specify the installation folder of the three libraries

.. code-block:: text
  :class: copy-button

  cmake ../dealii -DDEAL_II_WITH_MPI=ON -DDEAL_II_WITH_TRILINOS=ON -DTRILINOS_DIR=path/to/your/trilinos/installation -DDEAL_II_WITH_P4EST=ON -DP4EST_DIR=path/to/your/p4est/installation  -DDEAL_II_WITH_METIS=ON -DMETIS_DIR=path/to/your/metis/installation -DCMAKE_INSTALL_PREFIX=/path/to/desired/installation`

Compile deal.II

.. code-block:: text
  :class: copy-button

  make -j<nprocessor> install

Create an environment variable for the deal.II directory

.. code-block:: text
  :class: copy-button

  export DEAL_II_DIR=/path/to/dealii/installation

It is generally recommended to add the variable to your bashrc so it is always loaded:

.. code-block:: text
  :class: copy-button

  echo "export DEAL_II_DIR=/path/to/dealii/installation" >> ~/.bashrc

.. _install-lethe:

Installing lethe (Step #2)
-------------------------------

Clone lethe from the `Lethe official repository <https://github.com/lethe-cfd/lethe>`_.

.. code-block:: text
  :class: copy-button

  git clone https://github.com/lethe-cfd/lethe 

Create a build folder at the same level as the lethe folder

.. code-block:: text
  :class: copy-button

  mkdir build
  cd build

Build Lethe choosing the compilation option (Debug or Release). You can also optionally specify a path to an installation directory of your choice. We recommend that you do so, since this makes using Lethe much more comfortable.

.. code-block:: text
  :class: copy-button

  cmake ../lethe -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/home/username/path/to/installation

or

.. code-block:: text
  :class: copy-button

  cmake ../lethe -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/username/path/to/installation

Then you can compile:

.. code-block:: text
  :class: copy-button

  make -j<numprocs>

Testing your installation (Step #3)
-------------------------------------

Lethe comes pre-packaged with an extensive test suit for all of its modules. It can be used to test the validity of your installation. Within the build folder, the test suite can be launched with the following command:

.. code-block:: text
  :class: copy-button

  ctest -j $numprocs

where $numprocs can be the number of physical cores on your machine.

.. _update-dealii:

Updating deal.II
-------------------

Through the git repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The deal.II version supported by Lethe is updated and tested every week or so, see the repository `here <https://github.com/lethe-cfd/dealii>`_. If Lethe was installed with this forked version of deal.II, updating your deal.II installation is as simple as pulling the repository and recompiling the deal.II library. If your deal.II was installed manually using the deal.II master repository, the same process can be used.

With candi
~~~~~~~~~~~~~
In the candi folder (for instance, ``/home/username/software/candi``), modify the ``candi.cfg`` to get the latest dealii version, by changing the ``DEAL_II_VERSION`` variable in the case of an official release, or by changing the ``STABLE_BUILD`` in the case of a development release. The ``candi.cfg`` should contain:

.. code-block:: text
  :class: copy-button

  # Install the following deal.II version:
  DEAL_II_VERSION=v9.3.0

  # Would you like to build stable version of deal.II?
  # If STABLE_BUILD=false, then the development version of deal.II will be  
  # installed.
  STABLE_BUILD=true
  #STABLE_BUILD=false

Run the command ``./candi.sh`` to install the new version of dealii.

In your ``/home/deal.ii-candi`` folder, you should have a new folder with the dealii updated version (specified in ``DEAL_II_VERSION``, or ``deal.II-master`` in the case of development version)

You might need to delete the build folder of Lethe and redo the installation process from scratch, but this is rarely the case.
