====================================
Installation on Apple ARM (M1 or M2)
====================================

.. figure:: ./images/apple.png
   :height: 100px

Lethe can be now be deployed on Apple ARM, namely M1 and M2 chips. The support for these chips is experimental, but all Lethe solvers can be deployed for this type of architecture. So far, we have found that Lethe performs very efficiently on Apple ARM architecture. 

The installation of Lethe consists in two steps:
1. Installation of deal.II using the Candi toolset
2. Installation of Lethe

The easiest way to install deal.II and its dependencies under Mac OS is through the Candi toolset. Consequently, this is the only procedure which is explained here. Even under these condition, we will see that the procedure is shaky at best.


Installing deal.II Using Candi (Step #1)
-----------------------------------------

To install the dependencies (mpi, p4est, trilinos and METIS) all together using Candi, the following `procedure <https://github.com/dealii/dealii/wiki/Apple-ARM-M1-OSX>`_ on the deal.II wiki can be followed.

Clone the candi git repository in a folder of your choice  (e.g. ``$HOME/sofware/``). You can edit the ``candi.cfg`` file if you want to force the installation of the deal.II master version instead of the current stable version by setting the ``STABLE_BUILD=false``. Under Apple ARM, we only recommend the installation of the required libraries, namely parmetis, trilinos and p4est.

From the candi folder, the installation of candi can be launched using:

.. code-block:: text
  :class: copy-button

  ./candi.sh -j $num_proc --prefix=$path --packages="parmetis p4est trilinos dealii"

After installation, you should have a ``deal.ii-candi`` folder in the installation prefix directory, with the dealii folder of the desired version (see section :ref:`update-dealii`), as well as the required dependencies (p4est, trilinos, etc.).

After installation, add an environment variable to your ``~/.zshrc`` either manually or through the following command:

.. code-block:: text
  :class: copy-button

  echo "export DEAL_II_DIR=/home/username/deal.ii-candi/deal.II-<version>" >> ~/.zshrc

Setting the Library and Include Paths
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The deal.II installation procedure might not set the correct path for the libraries agaisnt which it needs to link. These include parmetis, p4est and trilinos. To fix this issue, the include and library path must be manually added by appending the following lines to the ``~/.zshrc`` file.

.. code-block::
  :class: copy-button

  export PATH=/opt/homebrew/bin:$PATH
  export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:path/to/candi/install/trilinos-release-12-18-1/lib/
  export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:path/to/candi/install/p4est-2.3.2/FAST/lib
  export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:path/to/candi/install/p4est-2.3.2/DEBUG/lib
  export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:path/to/candi/install//parmetis-4.0.3/lib
  export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:path/to/candi/install/trilinos-release-12-18-1/include/
  export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:path/to/candi/install/p4est-2.3.2/FAST/include/
  export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:path/to/candi/install/parmetis-4.0.3/lib/

.. note::
  This is not a clean workaround, but so far this is the only solution we have found. We welcome any suggestions on that front!


Fixing Trilinos Includes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lethe compiles using the C++ 17 standard. There is a header file in the Trilinos library which does not adhere to this standard and uses a deprecated function. When compiling with GCC, this does not pose any problem, but it does when using the built-in Clang compiler of Mac OS.

The following file needs to be modified in the include dir of your Trilinos installation:

.. code-block::
  :class: copy-button

  path/to/candi/install/trilinos-release-12-18-1/include/Tpetra_Import_def.hpp

At line 1169, the following line:

.. code-block:: cpp

  const size_type numInvalidRemote =
  std::count_if (remoteProcIDs.begin (), remoteProcIDs.end (),
                 std::bind1st (std::equal_to<int> (), -1));

Should be replaced by:

.. code-block:: cpp
  :class: copy-button

  const size_type numInvalidRemote = 0;

.. warning::
  It is still unclear to us what are the consequences of doing this, but we have found so far that it does not affect the outcome of the Lethe solvers. We are aware that this type of manual manipulation is far from being clean...

Numdiff
~~~~~~~~

numdiff is used within the automatic testing procedure of Lethe to compare files obtained through floating point arithmetic. Without numdiff, Lethe automatic tests may fail when they should not. numdiff can be installed directly from your package manager.

.. code-block:: text
  :class: copy-button

  brew install numdiff




Installation of Lethe (Step #2)
-------------------------------

Clone lethe from the `official repository <https://github.com/lethe-cfd/lethe>`

.. code-block:: text
  :class: copy-button

  git clone https://github.com/lethe-cfd/lethe 

Create a build folder at the same level as the lethe folder

.. code-block:: text
  :class: copy-button

  mkdir build
  cd build

Compile Lethe choosing the compilation option (Debug or Release). You can also optionally specify a path to an installation directory of your choice. We recommend that you do so, since this makes using Lethe much more comfortable.

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
