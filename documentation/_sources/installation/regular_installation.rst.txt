==============================
Regular Installation on Linux
==============================

.. figure:: ./images/linux.png
   :height: 100px

.. important::
  Distributions on which compatibility was tested are: Ubuntu 20.04 LTS, Ubuntu 22.04 LTS, Centos 7 and Manjaro.

Lethe requires a modern version of the `deal.II library <https://www.dealii.org/>`_ and its dependencies (MPI, numdiff, p4est, trilinos and METIS). At the time of this writing, ``deal.II 9.5`` and ``deal.II 9.6pre`` (the ``master`` branch version) are supported. A `dealii fork <https://github.com/chaos-polymtl/dealii>`_ is maintained by the Lethe team. This fork does not include any modification to deal.II library, but it is the latest version with which Lethe was tested. 

**Lethe installation steps:**
  
1. Installing deal.II  
2. :ref:`install-lethe`

**Installing deal.II and its dependencies:**
  
1. :ref:`install-deal.II-apt` (recommended for users)
2. :ref:`install-deal.II-candi` (recommended for developers) 
3. :ref:`install-deal.II-manually` (recommended for experienced developers)


.. _install-deal.II-apt:

Installing deal.II using apt 
-----------------------------------------

This is done following `this procedure <https://www.dealii.org/download.html#:~:text=page%20for%20details.-,Linux%20distributions,-Arch%20Linux>`_.

In case you are using Ubuntu, you will need to `update the backports <https://launchpad.net/~ginggs/+archive/ubuntu/deal.ii-9.5.1-backports>`_:

.. code-block:: text
  :class: copy-button

  sudo add-apt-repository ppa:ginggs/deal.ii-9.5.1-backports
  sudo apt update

A dependency required by Lethe, and that deal.II needs to be compiled with, is muParser:

.. code-block:: text
  :class: copy-button

  sudo apt-get install libmuparser-dev

To install deal.II, run:

.. code-block:: text
  :class: copy-button

  sudo apt-get install libdeal.ii-dev

To verify if the correct version of deal.II is installed, run:

.. code-block:: text
  :class: copy-button

  apt show libdeal.ii-dev

This should output several information about the installed version. Everything worked as expected if ``deal.ii-9.5.1`` is output

.. note::

  If the installed version is other than ``deal.ii-9.5.1``, follow `this link <https://github.com/dealii/dealii/wiki/Getting-deal.II>`_.


.. _install-deal.II-candi:

Installing deal.II using Candi 
-----------------------------------------

To install the dependencies (MPI, p4est, trilinos and METIS) all together using candi, the `procedure <https://github.com/dealii/candi.git>`_ on the candi repository can be followed.

Clone the candi git repository in a folder of your choice  (e.g. ``/home/username/software``). Edit the ``candi.cfg`` file to alter which dependencies are compiled. This file should notably be used to force the installation of the deal.II master version directly instead of the current stable version by setting ``DEAL_II_VERSION=master`` on line 97.

The following packages (which are specified after line 57) should be installed:
  
  .. code-block:: text
    
    PACKAGES="load:dealii-prepare"
    PACKAGES="${PACKAGES} once:numdiff"
    PACKAGES="${PACKAGES} once:opencascade"
    PACKAGES="${PACKAGES} once:parmetis"
    PACKAGES="${PACKAGES} once:p4est"
    PACKAGES="${PACKAGES} once:trilinos"
    PACKAGES="${PACKAGES} dealii"

A dependency required by Lethe, and that deal.II needs to be compiled with, is muParser:

.. code-block:: text
  :class: copy-button

  sudo apt-get install libmuparser-dev

Other packages can be disabled by simply commenting out the lines (adding a ``#`` at the beginning of the lines)

To ensure that the Lethe test suite works, deal.II must be configured with p4est version 2.3.6, the current default candi version of p4est. Otherwise, application tests that include restart files will fail.

From the candi folder, the installation of candi can be launched using:

.. code-block:: text
  :class: copy-button

  ./candi.sh -j$numproc --prefix=$path


where ``$numproc`` is the number of threads you want to use to compile deal.II and ``$path`` the installation prefix that is desired (e.g. ``/home/username/software/candi``).

.. tip:: 
  For a computer with 8Gb of RAM, 1 thread (``numproc=1``) should be used. For 16 Gb, 4 threads is reasonable. For 32 Gb, 16 threads or more can be used.


After installation, you should have a ``deal.II-candi`` folder in the installation prefix directory, with the dealii folder of the desired version (see section :ref:`update-dealii`), as well as the required dependencies (p4est, trilinos, etc.).

After installation, add the following lines variable to your ``.bashrc`` :

.. code-block:: text
  :class: copy-button
    
    source cand/install/prefix/configuration/enable.sh
    export DEAL_II_DIR=cand/install/prefix/deal.II-v.<version>


.. _install-lethe:

Installing Lethe 
-------------------------------

Clone Lethe from the `Lethe official repository <https://github.com/chaos-polymtl/lethe>`_.

.. code-block:: text
  :class: copy-button

  git clone https://github.com/chaos-polymtl/lethe --single-branch

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

  make -j$numprocs

.. warning:: 
  Hyper-threading will result in poor performance when running Lethe. Make sure that it is not turned on for your machine if you have more logical CPU(s) than the number of Core(s). To obtain this information, you can run the following command in a linux terminal:

  .. code-block:: text
    :class: copy-button

    lscpu | grep -E '^Thread|^Core|^Socket|^CPU\('

  To limit the number of thread you can set the environment variable ``OMP_NUM_THREADS`` to 1 in the terminal you are running the application or in your ``.bashrc`` file to make it the default:

  .. code-block:: text
    :class: copy-button

    export OMP_NUM_THREADS=1

Testing Your Installation 
-------------------------------------

Lethe comes pre-packaged with an extensive test suit for all of its modules. It can be used to test the validity of your installation. Within the build folder, the test suite can be launched with the following command:

.. code-block:: text
  :class: copy-button

  ctest -j$numprocs

where $numprocs can be the number of physical cores on your machine.

.. warning:: 
  The lethe test suites requires that deal.II be configured with p4est 2.3.6, otherwise the test that include restart files will fail.


.. _install-deal.II-manually:

Installing deal.II manually 
-----------------------------------------

Clone deal.II from the `deal.ii official repository <https://github.com/dealii/dealii>`_

.. code-block:: text
  :class: copy-button

  git clone https://github.com/dealii/dealii 

Configure deal.II in a build folder at the same level as the source code

.. code-block:: text
  :class: copy-button

  mkdir build
  cd build

Depending on how you have installed p4est, Trilinos and METIS, you may need to specify the installation folder of the three libraries. You also need to ensure that muParser is installed.

.. code-block:: text
  :class: copy-button

  cmake ../dealii -DDEAL_II_WITH_MPI=ON -DDEAL_II_WITH_TRILINOS=ON -DTRILINOS_DIR=path/to/your/trilinos/installation -DDEAL_II_WITH_P4EST=ON -DP4EST_DIR=path/to/your/p4est/installation  -DDEAL_II_WITH_METIS=ON -DMETIS_DIR=path/to/your/metis/installation -DCMAKE_INSTALL_PREFIX=/path/to/desired/installation`

Compile deal.II

.. code-block:: text
  :class: copy-button

  make -j$numprocs install

Create an environment variable for the deal.II directory

.. code-block:: text
  :class: copy-button

  export DEAL_II_DIR=/path/to/dealii/installation

It is generally recommended to add the variable to your bashrc so it is always loaded:

.. code-block:: text
  :class: copy-button

  echo "export DEAL_II_DIR=/path/to/dealii/installation" >> ~/.bashrc


.. _update-dealii:

Updating deal.II
-------------------

Through apt
~~~~~~~~~~~~~~~~~~~~~~~~~~~

As all other ``apt`` packages, run:

.. code-block:: text
  :class: copy-button

  sudo apt update
  sudo apt upgrade -y

Through the Git Repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The deal.II version supported by Lethe is updated and tested every week or so, see the repository `here <https://github.com/chaos-polymtl/dealii>`_. If Lethe was installed with this forked version of deal.II, updating your deal.II installation is as simple as pulling the repository and recompiling the deal.II library. If your deal.II was installed manually using the deal.II master repository, the same process can be used.

With Candi
~~~~~~~~~~~~~
In the candi folder (for instance, ``/home/username/software/candi``), modify the ``candi.cfg`` to get the latest dealii version, by changing the ``DEAL_II_VERSION`` variable in the case of an official release with its number (e.g. ``v9.5.2``), or by changing it to ``master`` in the case of a development release. The ``candi.cfg`` file should contain on lines 96-97:

.. code-block:: text
  :class: copy-button

  # Install the following deal.II version (choose master, v9.3.0, v9.2.0, ...)
  DEAL_II_VERSION=master

Run the command ``./candi.sh`` to install the new version of dealii.

In your ``/home/deal.ii-candi`` folder, you should have a new folder with the dealii updated version (specified in ``DEAL_II_VERSION``, or ``deal.II-master`` in the case of development version)

You might need to delete the build folder of Lethe and redo the installation process from scratch, but this is rarely the case.
