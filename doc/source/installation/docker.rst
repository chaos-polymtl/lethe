======
Docker
======

.. note::

    You will need to install `Docker <https://www.docker.com/get-started>`_ and have 5 GB of free disk space.



##############################################################
No Compilation Required: Using the Provided Lethe Docker Image
##############################################################

If you don't want to build Lethe and its dependencies, you can use the provided `Docker image <https://github.com/lethe-cfd/lethe/pkgs/container/lethe>`_.

For example, to launch the 2D Lid-Driven Cavity Flow simulation, run the following lines inside the root Lethe folder:

.. code-block:: shell
  :class: copy-button

  docker run --rm \
    -v $(pwd):/home/dealii \
    ghcr.io/lethe-cfd/lethe:master \
    lethe-fluid examples/incompressible_flow/2d_lid_driven_cavity/cavity.prm

In general, to run a lethe simulation you will just need to run the following command:

.. code-block:: shell
  :class: copy-button

  docker run --rm \
    -v $(pwd):/home/dealii \
    ghcr.io/lethe-cfd/lethe:master \
    lethe-executable parameter-file.prm

#################################
Using a deal.II  Docker Container
#################################

The deal.II library is a mature project containing `O(1,000,000)` lines of C++ code, itself depending on other similar libraries (e.g. Trilinos). Compiling deal.II with its exact dependencies can be a daunting task even for moderate C++ programmers - one with many knobs and dials, normally taking 4-6 hours.

Thankfully, the deal.II dev heroes also maintain some `Docker Hub <https://hub.docker.com/r/dealii/dealii/>`_ images containing pre-built deal.II libraries in a Ubuntu environment. For those unfamiliar with Docker, it allows the creation of completely isolated environments, most often pre-configured for specific apps; it is almost like a virtual machine, but running natively on your current OS (the technical term is *OS-level* virtualization, rather than *hardware-level*) - so scientific programs shouldn't have any degradation in performance. 

Below are some instructions for getting pre-configured deal.II environments, reducing the time needed to install lethe's dependencies from hours to minutes. There are 3 main parts in this tutorial:

- If you have experience with Docker, the :ref:`tl-dr` section should be enough.
- The next two sections have a brief introduction to Docker with enough information to get you started with deal.II containers.
- The final :ref:`saving-docker` section has a shell script you can run directly to start and manage deal.II containers.

Sometimes it's easier to effectively set up a completely new OS than compile millions of lines of C++...


.. _tl-dr:

0. TL;DR
--------

In some directory run:

.. code-block:: bash
  :class: copy-button

    local:~$ mkdir lethe-simulations && cd lethe-simulations
    local:~/lethe-simulations$ docker run -it --name lethe-container -v $(pwd):/home/docker-host dealii/dealii:master-focal
    To run a command as administrator (user "root"), use "sudo <command>".
    See "man sudo_root" for details.

    dealii@2e92a38a6223:~$ cd /home/docker-host/
    dealii@2e92a38a6223:/home/docker-host/$ # compile lethe and run simulations here

And do the final lethe compilation step.
This starts a new docker container with an interactive terminal (:code:`-it`) named :code:`lethe-container`, mounting the current directory on your machine to :code:`/home/docker-host`. Inside the container, run any simulations in the :code:`/home/docker-host` directory to save them locally too.

Exit the docker container with :code:`ctrl-d`; run :code:`docker ps -a` to see the current containers on your machine - you'll see the :code:`lethe-container`. Restart your container and attach to its terminal with:

.. code-block:: bash
  :class: copy-button

  docker start lethe-container && docker attach lethe-container

Have fun.


1. Get Docker
-------------

First, install Docker; complete instructions for Windows, Mac or Linux can be found `here <https://docs.docker.com/get-docker/>`_; this should be a relatively quick step.


2. Launch a deal.II Container
-----------------------------

Some Docker basics:

- A Docker **image** is a "frozen", shareable environment - e.g. a Ubuntu terminal with pre-installed deal.II and its dependencies.
- You can launch a **container** running that specific image, in which you can e.g. compile lethe and run simulations. Indeed, you can run multiple separate containers of the same image.

In your terminal, you can then run:

.. code-block:: bash
  :class: copy-button

  docker run -it --name lethe-container dealii/dealii:master-focal

This starts a new container with an interactive terminal (:code:`-it`) running the :code:`dealii/dealii:master-focal` image (`see here <https://hub.docker.com/r/dealii/dealii/tags>`_ for all the deal.II images) named :code:`lethe-container`. Press :code:`ctrl-d` to exit the container.

To see your current docker containers, run:

.. code-block:: bash
  :class: copy-button

  docker ps -a

  CONTAINER ID   IMAGE                        COMMAND   CREATED          STATUS                      PORTS     NAMES
  e3a7f71639f6   dealii/dealii:master-focal   "bash"    14 minutes ago   Exited (0) 35 seconds ago             lethe-container

This container saved your changes. You can restart and attach to the container's terminal by running:

.. code-block:: bash
  :class: copy-button

  docker start lethe-container && docker attach lethe-container

If you want to, you can remove the container with :code:`docker rm lethe-container`; you'll start a new fresh container by running the :code:`docker run...` command above.

However, any files saved in the container are only accessible inside it, and are lost when removing the container. For simulations on the other hand, we want their outputs to be saved and accessible on the local machine (e.g. to post-process them); for this, we will *mount* a directory from the local machine to the container with :code:`-v LOCAL_DIR:CONTAINER_DIR`. While in the container, anything you save to :code:`CONTAINER_DIR` will be accessible on your local machine in :code:`LOCAL_DIR`.

For example, on your local machine:

.. code-block:: bash
  :class: copy-button

  local:~$ mkdir ~/lethe-simulations
  local:~$ cd ~/lethe-simulations
  local:~/lethe-simulations$ ls

  local$ docker run -it --name lethe-container -v $(pwd):/home/docker-host dealii/dealii:master-focal
  To run a command as administrator (user "root"), use "sudo <command>".
  See "man sudo_root" for details.

  dealii@2e92a38a6223:~$ cd /home/docker-host/
  dealii@2e92a38a6223:/home/docker-host$ echo "Hello lethe!" > somefile.txt
  dealii@2e92a38a6223:/home/docker-host$ exit

  local:~/lethe-simulations$ ls
  somefile.txt

That's all the Docker-specific tutorial! Launch your container running a deal.II image, go to :code:`/home/docker-host` to save your changes locally too, download lethe, compile it, and run your simulations there.


.. _saving-docker:

3. Saving Docker Commands in a Bash Script
------------------------------------------

We can add all the commands above, plus some comments and helpful messages to a single shell script named `docker_lethe.sh`:

.. code-block:: bash
  :class: copy-button

  #!/bin/sh

  # Launch a persistent docker container from a given image, automatically re-attaching to it on
  # future runs.
  #
  # The current local directory is mounted in /home/docker-host within the container; run any
  # simulations there to save results on the local machine's current directory.

  DOCKER_IMAGE='dealii/dealii:master-focal'
  DOCKER_CONTAINER="$USER-${DOCKER_IMAGE##*/}"        # Remove repository prefix
  DOCKER_CONTAINER=${DOCKER_CONTAINER/:/-}            # Replace : with -

  HOST_DIR=$(pwd)
  REMOTE_DIR="/home/docker-host"


  printf "Image:     ${DOCKER_IMAGE}\nContainer: ${DOCKER_CONTAINER}\n\n"


  # If a container with this name already exists, re-attach to it
  if [ "$(docker ps -q -af name=${DOCKER_CONTAINER})" ]
  then
      printf "Found previous container with same name; starting and attaching...\n\n"
      docker start ${DOCKER_CONTAINER}
      docker attach ${DOCKER_CONTAINER}
  else
      printf "Launching new container...\n\n"
      docker run -it \
          --name $DOCKER_CONTAINER \
          -v $HOST_DIR:$REMOTE_DIR \
          $DOCKER_IMAGE
  fi

Then just execute the shell script:

.. code-block:: bash
  :class: copy-button

  local:~/lethe-simulations$ sh docker_lethe.sh


Final Notes
-----------

You can now download, run and manage Docker containers pre-configured with deal.II; it is a powerful tool that you can use for any other projects as well, without polluting your main programming environment or spending hours figuring out the specific libraries needed (`dependency hell <https://en.wikipedia.org/wiki/Dependency_hell>`_ is real).

You can now clone ``lethe``, compile it, and run large-scale, efficient multi-physics simulations!

