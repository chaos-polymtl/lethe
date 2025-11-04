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
* ``$HOME/dealii/inst`` for installation (``ninja install`` command)

Folders can be open with the ``cd`` command (``cd $folder_path``).

For the sake of clarity, this is the folder structure considered for the rest of this tutorial.

Installing deal.II 
------------------

All operations must be performed on login nodes.

Load ``Trilinos``, ``Parmetis`` and ``P4est``, and their prerequisite modules and set the appropriate environment variables. It is convenient to create a ``.dealii`` file in your ``$HOME`` directory that contains the following lines to source the appropriate libraries:

.. code-block:: text
  :class: copy-button

  module load StdEnv/2023
  module load trilinos/15.1.1
  module load parmetis/4.0.3
  module load p4est/2.8.6
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$EBROOTTRILINOS/lib64/
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$EBROOTP4EST/lib64/
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$EBROOTFLEXIBLAS/lib64/
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$EBROOTIMKL//mkl/2023.2.0/lib/intel64/
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$EBROOTOPENMPI/lib/
  source /scinet/vast/etc/vastpreload-openmpi.bash # Only if on Trillium or Nibi

  export DEAL_II_DIR=$HOME/dealii/inst/
  export PATH=$PATH:$HOME/lethe/inst/bin/

This file needs to be sourced every time you launch a job or you compile deal.II and/or Lethe. Once the file has been created, you can then source it on the terminal with:

.. code-block:: text
  :class: copy-button

  source $HOME/.dealii

and use it in your ``.sh`` script when launching a job (see :ref:`Launching Simulations<Launching Simulations>` below). 

Although Lethe always supports the master branch of deal.II, we maintain an identical deal.II fork on the CHAOS laboratory organization. This fork is always tested to make sure it works with Lethe. To clone this deal.II fork, execute in ``$HOME/dealii`` directory:

.. code-block:: text
  :class: copy-button

  git clone https://github.com/chaos-polymtl/dealii.git

We can compile ``dealii`` in the ``$HOME/dealii/build`` folder, by defining the paths to installation folders of ``Trilinos``, ``Parmetis`` and ``P4est``. To increase the speed of this step, we skip ``dealii`` tests and compile in release mode only.

.. code-block:: text
  :class: copy-button

  cmake ../dealii -DDEAL_II_WITH_MPI=ON -DDEAL_II_WITH_TRILINOS=ON   -DTRILINOS_DIR=$EBROOTTRILINOS  -DDEAL_II_WITH_P4EST=ON  -DDEAL_II_WITH_METIS=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/dealii/inst/ -DDEAL_II_COMPONENT_EXAMPLES=OFF  -DCMAKE_CXX_FLAGS="-march=native" -G Ninja

.. tip::

  The -DCMAKE_CXX_FLAGS="-march=native" works on both Rorqual, Nibi, and Trillium. To ensure that the flag has worked correctly, the cmake output should contain the following information : ``Vectorization level:    512 bit (sse2 avx2 avx512*)``.

.. warning::

  If you wish to run simulations with over 4B (:math:`4\cdot 10^9`) degrees of freedom, you must compile with the ``DEAL_II_WITH_64BIT_INDICES = ON`` flag. Such large simulations should be carried out using the ``lethe-fluid-matrix-free`` application.

and:

.. code-block:: text
  :class: copy-button

  nice ninja -j6 install

The argument ``-jX`` specifies the number of processors used for the compilation. On login nodes, a maximum of 6 cores should be used in order to ensure that other users can continue using the cluster without slowdowns. If you use more than 6 cores, your compilation may be terminated automatically.

Installing Lethe
----------------

After installing deal.II, compiling Lethe is relatively straightforward. To compile Lethe, the ``Trilinos``, ``Parmetis`` and ``P4est`` modules should be loaded.

In the ``$HOME/lethe`` directory, download Lethe:

.. code-block:: text
  :class: copy-button

  git clone https://github.com/chaos-polymtl/lethe.git 

To install Lethe in the ``$HOME/lethe/inst`` directory (applications will be in ``inst/bin``), run in the ``$HOME/lethe/build`` directory:

.. code-block:: text
  :class: copy-button

  cmake ../lethe  -DDEAL_II_DIR=$HOME/dealii/inst -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../inst -DCMAKE_CXX_FLAGS="-march=native" -G Ninja
  nice ninja -j6 install


.. _copying-local-files:

Copying Local Files
-------------------

We use `Globus <https://alliancecan.ca/en/services/advanced-research-computing/national-services/data-movement-globus>`_ to transfer files between your local machine and the cluster. For more information, visit the `Globus documentation <https://docs.alliancecan.ca/wiki/Globus>`_.

.. _Launching Simulations:

Launching Simulations
---------------------

Simulations are sent to the scheduler via batch scripts. Visit the Digital Research Alliance of Canada (Alliance) wiki page for more information about the `scheduler <https://docs.alliancecan.ca/wiki/What_is_a_scheduler%3F>`_ and `running jobs <https://docs.alliancecan.ca/wiki/Running_jobs>`_. For your convenience, an example of ``job.sh`` is given below:

.. code-block:: text
  :class: copy-button

  #!/bin/bash
  #SBATCH --account=$yourgroupaccount
  #SBATCH --ntasks-per-node=$X #number of parallel tasks per node.
  #SBATCH --nodes=1 #number of whole nodes used 
  #SBATCH --time=1:00:00 #maximum time for the simulation (hh:mm:ss)
  #SBATCH --mem=120G #memory usage per node. See cluster specification for maximal amount.
  #SBATCH --job-name=$yourjobname
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=$your.email.adress@email.provider

  source $HOME/.dealii
  mpirun $HOME/lethe/inst/bin/$lethe_application_name_wanted $parameter_file_name.prm


.. tip::
  The ``--ntasks-per-node`` option is the number of parallel tasks per node. When using a full node, this should correspond to the number of cores available on the node. For example, on Narval, this should be set to 64.

.. tip::
    If you have jobs that need to be launched one after the other, you can add ``#SBATCH --dependency=$previous-slurm-job-id`` to your launching script. This will make sure that the job will only start once the previous job has finished.

The job is sent using:

.. code-block:: text
  :class: copy-button

  sbatch job.sh

Status can be followed with the ``sq`` command: under ``ST``, ``PD`` indicates a pending job, and ``R`` a running job.

Console outputs are written in ``slurm-$jobID.out``. For instance, to display the 20 last lines from this file, use:

.. code-block:: text
  :class: copy-button

  tail -n 20 slurm-$jobID.out

Clusters Specifications
------------------------

Please consult the documentation for the machine you are using for the specification of the nodes: 

+-----------------+---------------------+---------------------+----------------------------------------------+
| Cluster         | Tasks per Node      | Memory per Node     | URL                                          |
+=================+=====================+=====================+==============================================+
| Narval          | 64                  | 248 Go              | https://docs.alliancecan.ca/wiki/Narval/en   |
+-----------------+---------------------+---------------------+----------------------------------------------+
| Trillium        | 192                 | 755 Go              | https://docs.alliancecan.ca/wiki/Trillium/en |
+-----------------+---------------------+---------------------+----------------------------------------------+
| Rorqual         | 192                 | 760 Go              | https://docs.alliancecan.ca/wiki/Rorqual/en  |
+-----------------+---------------------+---------------------+----------------------------------------------+
| Nibi            | 192                 | 754 Go              | https://docs.alliancecan.ca/wiki/Nibi/en     |
+-----------------+---------------------+---------------------+----------------------------------------------+

Saving a SSH Key (Linux)
------------------------

To save your key on the cluster, so that it is not asked for each log, generate your ssh-key with:

.. code-block:: text
  :class: copy-button

  ssh-keygen

which defaults to an RSA key. If you want to specify the key type you want to generate (i.e. ED25519 key), type

.. code-block:: text
  :class: copy-button

  ssh-keygen -t ed25519

.. note::
  ED25519 keys are preferred to RSA keys since they are more secure and performant. Seek more information in the `GitLab Documentation <https://docs.gitlab.com/ee/user/ssh.html>`_.

To upload this local key to your Compute Canada Database account (CCDB) use:

.. code-block:: text
  :class: copy-button

    ssh-copy-id username@clustername.computecanada.ca