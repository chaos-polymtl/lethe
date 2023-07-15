=============================================================
Running Lethe on Digital Research Alliance of Canada Clusters
=============================================================


Setting-up the Folder Structure
-------------------------------

In your ``$HOME``, create a "dealii" folder and a "lethe" folder, each containing "build" and "inst" folders:

.. code-block:: text
  :class: copy-button

  mkdir -p {dealii,lethe}/{build,inst}

The deal.II and Lethe projects can then be cloned in their corresponding folders, as indicated later in this tutorial.

After installation is complete, the folder structure will be, for deal.II (and likewise for Lethe):

* ``$HOME/dealii/dealii`` for deal.ii git,
* ``$HOME/dealii/build`` for compilation (``cmake`` command),
* ``$HOME/dealii/inst`` for installation (``make install`` command)

Folders can be open with the ``cd`` command (``cd $folder_path``).

For the sake of clarity, this is the folder structure considered for the rest of this tutorial.

Installing deal.II
------------------

On Niagara, Beluga, Narval, Graham or Cedar
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All operations must be performed on login nodes.

.. warning::
 If on **Niagara**, this additional module is needed beforehand: ``module load CCEnv``.

Load ``Trilinos``, ``Parmetis`` and ``P4est``, and their prerequisite modules:

.. code-block:: text
  :class: copy-button

  module load CCEnv #if on Niagara
  module load StdEnv/2020
  module load cmake/3.23.1
  module load gcc/9.3.0
  module load openmpi/4.0.3
  module load trilinos/13.0.1
  module load parmetis/4.0.3
  module load p4est/2.2
  module load opencascade/7.5.2

Then, we can clone and compile ``dealii``. Although Lethe always supports the master branch of deal.II, we maintain an identical deal.II fork on the lethe repository. This fork is always tested to make sure it works with lethe. To clone the deal.II fork github repository, execute in ``$HOME/dealii`` directory:

.. code-block:: text
  :class: copy-button

  git clone https://github.com/lethe-cfd/dealii.git

We can compile ``dealii`` in the ``$HOME/dealii/build`` folder, by defining the paths to installation folders of ``Trilinos``, ``Parmetis`` and ``P4est``. To increase the speed of this step, we skip ``dealii`` tests and compile in release mode only.

.. code-block:: text
  :class: copy-button

  cmake ../dealii -DDEAL_II_COMPONENT_EXAMPLES=OFF -DCMAKE_BUILD_TYPE=Release -DDEAL_II_WITH_OPENCASCADE=ON -DOPENCASCADE_DIR=$EBROOTOPENCASCADE -DDEAL_II_WITH_MPI=ON -DDEAL_II_WITH_TRILINOS=ON -DTRILINOS_DIR=$EBROOTTRILINOS -DDEAL_II_WITH_P4EST=ON -DP4EST_DIR=$EBROOTP4EST -DDEAL_II_WITH_METIS=ON -DMETIS_DIR=$EBROOTPARMETIS -DCMAKE_INSTALL_PREFIX=../inst -DTrilinos_FIND_COMPONENTS="Pike;PikeImplicit;PikeBlackBox;TrilinosCouplings;Panzer;PanzerMiniEM;PanzerAdaptersSTK;PanzerDiscFE;PanzerDofMgr;PanzerCore;Piro;ROL;Stokhos;Tempus;Rythmos;ShyLU;ShyLU_DD;ShyLU_DDCommon;ShyLU_DDFROSch;ShyLU_DDBDDC;Zoltan2;Zoltan2Sphynx;MueLu;Moertel;NOX;Phalanx;Percept;STK;STKExprEval;STKDoc_tests;STKUnit_tests;STKBalance;STKTools;STKTransfer;STKSearchUtil;STKSearch;STKUnit_test_utils;STKNGP_TEST;STKIO;STKMesh;STKTopology;STKSimd;STKUtil;STKMath;Compadre;Intrepid2;Intrepid;Teko;FEI;Stratimikos;Ifpack2;Anasazi;Komplex;SEACAS;SEACASEx2ex1v2;SEACASTxtexo;SEACASNumbers;SEACASNemspread;SEACASNemslice;SEACASMat2exo;SEACASMapvar-kd;SEACASMapvar;SEACASMapvarlib;SEACASExplore;SEACASGrepos;SEACASGenshell;SEACASGen3D;SEACASGjoin;SEACASFastq;SEACASEx1ex2v2;SEACASExo_format;SEACASExotxt;SEACASExomatlab;SEACASExodiff;SEACASExo2mat;SEACASEpu;SEACASEjoin;SEACASConjoin;SEACASBlot;SEACASAprepro;SEACASAlgebra;SEACASPLT;SEACASSVDI;SEACASSuplibCpp;SEACASSuplibC;SEACASSuplib;SEACASSupes;SEACASAprepro_lib;SEACASChaco;SEACASIoss;SEACASNemesis;SEACASExoIIv2for32;SEACASExodus_for;SEACASExodus;Amesos2;ShyLU_Node;ShyLU_NodeTacho;ShyLU_NodeHTS;Belos;ML;Ifpack;Zoltan2Core;Pamgen;Amesos;Galeri;AztecOO;Pliris;Isorropia;Xpetra;Thyra;ThyraTpetraAdapters;ThyraEpetraExtAdapters;ThyraEpetraAdapters;ThyraCore;Domi;TrilinosSS;Tpetra;TpetraCore;TpetraTSQR;TpetraClassic;EpetraExt;Triutils;Shards;Zoltan;Epetra;MiniTensor;Sacado;RTOp;KokkosKernels;Teuchos;TeuchosKokkosComm;TeuchosKokkosCompat;TeuchosRemainder;TeuchosNumerics;TeuchosComm;TeuchosParameterList;TeuchosParser;TeuchosCore;Kokkos;KokkosAlgorithms;KokkosContainers;KokkosCore;Gtest;TrilinosATDMConfigTests;TrilinosFrameworkTests"

This command is absolutely barbaric, but it is required to prevent the usage of the Zadelus library within Trilinos 13.0.1.

and:

.. code-block:: text
  :class: copy-button

  make -j8 install

The argument ``-jX`` specifies the number of processors used for the compilation. On login nodes, a maximum of 8 cores should be used in order to ensure that other users can continue using the cluster without slowdowns. In addition, the ``nice`` command can be used at the beginning of the call to give a lower priority to this task.

Installing Lethe
----------------

After installing deal.II, compiling Lethe is relatively straightforward, especially since all of these clusters share a very similar environment. To compile Lethe, the ``Trilinos``, ``Parmetis`` and ``P4est`` modules should be loaded.

In the ``$HOME/lethe`` directory, download Lethe:

.. code-block:: text
  :class: copy-button

  git clone https://github.com/lethe-cfd/lethe.git

To install Lethe in the ``$HOME/lethe/inst`` directory (applications will be in ``inst/bin``), run in the ``$HOME/lethe/build`` directory:

.. code-block:: text
  :class: copy-button

  cmake ../lethe  -DDEAL_II_DIR=../../dealii/inst -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../inst -DBUILD_TESTING=OFF
  make install -j8


Installing Numdiff to Enable Tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You will need to have numdiff installed to enable the test suite, otherwise you will have an error at the cmake step of Lethe's installation when using ``-DBUILD_TESTING=ON``, stating that this module is missing. To install the package manually use the following steps:

1. Download the `compressed folder <https://mirror.csclub.uwaterloo.ca/nongnu/numdiff/>`_ (ex/ numdiff-5.9.0.tar.gz)
2. Unzip it
3. Copy it with ``scp -r`` to your Compute Canada account on the chosen cluster (see :ref:`copying-local-files` section)
4. In the numdiff folder on the cluster, execute:

   .. code-block:: text
     :class: copy-button

     ./configure
     make

5. Add it to your path environment:

   .. code-block:: text
     :class: copy-button

     PATH=$PATH:$HOME/path/to/numdiff/folder


.. _copying-local-files:

Copying Local Files
-------------------

On Linux, use ``scp`` (for secure copy) to copy needed files for the simulation (``prm``, ``msh``):

.. code-block:: text
  :class: copy-button

  scp /home/path/in/your/computer/*.prm username@clustername.calculcanada.ca:/scratch/path/in/cluster

If you need to copy a folder, use ``scp -r``.

Simulation files must be in scratch. To get the address of your scratch folder, in your cluster account run:

.. code-block:: text
  :class: copy-button

  cd $SCRATCH
  pwd

On Windows, use third-party, such as ``PuTTY`` (see the `wiki page on Transferring data <https://docs.computecanada.ca/wiki/Transferring_data>`_))


Creating a .dealii
------------------

In order to call your deal.II local installation, it is convenient to create a ``.dealii`` file in your ``$HOME`` directory:

.. code-block:: text
  :class: copy-button

  nano .dealii

In the nano terminal, copy-paste (with ``Ctrl+Shift+V``):

.. code-block:: text
  :class: copy-button

  module load CCEnv #if on Niagara
  module load StdEnv/2020
  module load cmake/3.23.1
  module load gcc/9.3.0
  module load openmpi/4.0.3
  module load trilinos/13.0.1
  module load parmetis/4.0.3
  module load p4est/2.2
  module load opencascade/7.5.2

  export DEAL_II_DIR=$HOME/dealii/inst/
  export PATH=$PATH:$HOME/lethe/inst/bin/

Exit the nano mode with ``Ctrl+x`` and save the document by hitting ``y`` on the prompt "Save modify buffer?" (in the bottom). The prompt "File Name to Write: .dealii" should then appear, hit ``Enter``.

You can then source it on the terminal with:

.. code-block:: text
  :class: copy-button

  source $HOME/.dealii

and use it in your ``.sh`` script (see Launching simulations below).

Launching Simulations
---------------------

Simulations are sent to the scheduler via batch scripts. Visit the Digital Research Alliance of Canada (Alliance) wiki page for more information about the `scheduler <https://docs.alliancecan.ca/wiki/What_is_a_scheduler%3F>`_ and `running jobs <https://docs.alliancecan.ca/wiki/Running_jobs>`_. For your convenience, an example of ``job.sh`` used on Beluga is given below:

.. code-block:: text
  :class: copy-button

  #!/bin/bash
  #SBATCH --account=$yourgroupaccount
  #SBATCH --ntasks-per-node=$X #number of parallel tasks (as in mpirun -np X)
  #SBATCH --nodes=1 #number of whole nodes used (each with up to 40 tasks-per-node)
  #SBATCH --time=1:00:00 #maximum time for the simulation (hh:mm:ss)
  #SBATCH --mem=120G #memory usage per node. See cluster specification for maximal amount.
  #SBATCH --job-name=$yourjobname
  #SBATCH --mail-type=END #email preferences
  #SBATCH --mail-type=FAIL
  #SBATCH --mail-user=$your.email.adress@email.provider

  source $HOME/.dealii
  srun $HOME/lethe/inst/bin/$lethe_application_name_wanted $parameter_file_name.prm

The job is sent using:

.. code-block:: text
  :class: copy-button

  sbatch job.sh

Status can be followed with the ``sq`` command: under ``ST``, ``PD`` indicates a pending job, and ``R`` a running job.

Console outputs are written in ``slurm-$jobID.out``. For instance, to display the 20 last lines from this file, use:

.. code-block:: text
  :class: copy-button

  tail -n 20 slurm-$jobID.out

.. note::
 If you need to launch multiple simulations, such as with varying parameter, feel free to adapt one of the scripts provided on `lethe-utils <https://github.com/lethe-cfd/lethe-utils/tree/master/python/cluster>`_.


Saving a SSH Key (Linux)
------------------------

To save your key on the cluster, so that it is not asked for each log or ``scp``, generate your ssh-key with:

.. code-block:: text
  :class: copy-button

  ssh-keygen

which defaults to an RSA key. If you want to specify the key type you want to generate (i.e. ED25519 key), type

.. code-block:: text
  :class: copy-button

  ssh-keygen -t ed25519

.. note::
  ED25519 keys are preferred to RSA keys since they are more secure and performant. Seek more information in the `Gitlab Documentation<https://docs.gitlab.com/ee/user/ssh.html>`.

To upload this local key to your Compute Canada Database account (CCDB) use:

.. code-block:: text
  :class: copy-button

  ssh-copy-id username@clustername.computecanada.ca

.. warning::
 This command does not work on Niagara anymore. You may use the following:

 .. code-block:: text
  :class: copy-button

  cat ~/.ssh/$KEY_ID.pub

 where ``$KEY_ID.pub`` is the public key file located in ``~/.ssh/``. For more information, see `SSH documentation <https://docs.scinet.utoronto.ca/index.php/SSH#SSH_Keys>`_.
