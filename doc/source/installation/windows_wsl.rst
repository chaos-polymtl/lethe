================
Windows with WSL
================

.. figure:: ./images/windows.png
   :height: 100px

.. important::
  Distributions compatibility: Windows 10 and Windows 11

.. |linux_shell| image:: ./images/linux.png
   :height: 15px

.. |win_shell| image:: ./images/windows.png
   :height: 15px

.. seealso::

  This tutorial is aimed at Windows users who have no prior knowledge of Linux. It installs Lethe in Ubuntu 24.04 LTS, running in the Windows Subsystem for Linux (WSL). If you are a developer or need more options, see :doc:`regular_installation`.

Throughout this tutorial:
  * |win_shell| indicates operations performed in the Windows session, and
  * |linux_shell| indicates operations performed in the Linux subsystem.

.. tip::
  To execute a command on a shell (Ubuntu or Windows command prompt), type or copy/paste the given command and hit ``Enter``. Multiple commands are given in multiple lines, or separated by ``;``: when copying/pasting, they will be executed one after the other.

Installing WSL and Ubuntu (Step #0)
------------------------------------

1. |win_shell| Install WSL (Windows Subsystem for Linux) and Ubuntu 24.04 LTS. Open PowerShell or the Windows Command Prompt in administrator mode by right-clicking on it and selecting "Run as administrator", enter the following command, then restart your computer:

.. code-block:: text
  :class: copy-button

  wsl --install -d Ubuntu-24.04

.. admonition:: Verify the installed version of WSL

  In the Windows command prompt (Start menu > ``cmd``):

  .. code-block:: text
      :class: copy-button

      wsl -l -v

  should indicate ``2`` in the ``VERSION`` column. If not, follow `these instructions <https://learn.microsoft.com/en-us/windows/wsl/install#upgrade-version-from-wsl-1-to-wsl-2>`_ to upgrade to WSL 2.

2. |win_shell| Launch Ubuntu 24.04 LTS from the Start menu. The first time, you will be asked to choose a username and a password for Ubuntu. Then, |linux_shell| update Ubuntu:

.. code-block:: text
  :class: copy-button

  sudo apt update
  sudo apt upgrade

.. tip::
  The ``sudo`` command will ask you to type your Ubuntu password. Note that Linux does not show any symbol while typing a password, contrary to Windows with ``*``: simply type your password and press ``Enter``.

When prompted "Do you want to continue?", proceed by typing ``y`` and hitting ``Enter``.

3. |win_shell| (optional) For better ease in the Linux terminal (better coloring, multiple tabs), use ``Windows Terminal`` to launch Ubuntu. It is installed by default on Windows 11 and can be downloaded from the Microsoft Store on Windows 10. In its settings, you can select the Ubuntu profile as the default profile. In the terminal, use ``Ctrl+Shift+C`` and ``Ctrl+Shift+V`` to copy and paste text or commands.

.. tip::
  A (very) few Linux commands useful for navigation:
    * ``mkdir $dir``: (make directory) create a directory with the name specified as ``$dir``
    * ``cd $dir``: (change directory) move to the directory ``$dir``
    * ``cd ..``: move up to the parent directory
    * ``pwd``: (print working directory) return the directory you are in
    * ``cd $HOME``: move to your home directory (``/home/<user_name>/``)
    * ``explorer.exe .``: open the directory you are in with the Windows File Explorer

  You can find `here <https://linuxconfig.org/linux-commands>`_ a thorough guide for the most basic Linux commands.


The following step is to install deal.II, the finite element library on which Lethe is built. This can be done through

  1. Advanced Packaging Tool (apt) (**this is by far the easiest way to proceed**) : :ref:`install-deal.II apt` (recommended for users)

  2. Candi shell script (`candi github page <https://github.com/dealii/candi>`_), which compiles deal.II and its dependencies from source: :ref:`install-deal.II candi` (recommended for developers)

.. important::
  Since September 2026, Lethe requires deal.II to be compiled with `magic_enum <https://github.com/Neargye/magic_enum>`_. The deal.II package installed with apt already includes it, and deal.II 9.8, installed by candi, ships a bundled copy of it: no action is required.

.. _install-deal.II apt:

Installing deal.II using apt (Step #1)
-----------------------------------------

This is done following `this procedure <https://www.dealii.org/download.html#:~:text=page%20for%20details.-,Linux%20distributions,-Arch%20Linux>`_.

1. |linux_shell| The deal.II version provided by Ubuntu 24.04 is too old for Lethe. Add the `deal.II backports <https://launchpad.net/~ginggs/+archive/ubuntu/deal.ii-9.7.1-backports>`_, which provide a more recent version:

.. code-block:: text
  :class: copy-button

  sudo add-apt-repository ppa:ginggs/deal.ii-9.7.1-backports
  sudo apt update

2. |linux_shell| Install deal.II, along with the tools required to compile and test Lethe:

.. code-block:: text
  :class: copy-button

  sudo apt install libdeal.ii-dev build-essential cmake git numdiff

To verify if the correct version of deal.II is installed, run:

.. code-block:: text
  :class: copy-button

  apt show libdeal.ii-dev

The ``Version`` field of the output should start with ``9.7.1``.

.. note::

  If the installed version is other than ``9.7.1``, follow `this link <https://github.com/dealii/dealii/wiki/Getting-deal.II>`_.

You can now proceed to :ref:`install-lethe-wsl`.

.. _install-deal.II candi:

Installing deal.II using Candi (Step #1)
-----------------------------------------

.. important::
  Candi compiles deal.II and all its dependencies from source, which takes from 1 to 3 hours. Read and follow each step carefully.

1. |linux_shell| Install the packages required by candi and Lethe:

.. code-block:: text
  :class: copy-button

  sudo apt install lsb-release git subversion wget \
  bc libgmp-dev build-essential autoconf automake cmake \
  libtool gfortran libboost-all-dev zlib1g-dev openmpi-bin \
  openmpi-common libopenmpi-dev libblas3 libblas-dev \
  liblapack3 liblapack-dev libsuitesparse-dev numdiff

.. tip::
  The symbols ``\`` indicate that this a single command written on multiple lines.

2. |linux_shell| Create a ``Software`` folder and download candi in it:

.. code-block:: text
  :class: copy-button

  mkdir -p $HOME/Software; cd $HOME/Software
  git clone https://github.com/dealii/candi.git
  cd candi

Note the use of ``;`` which enables to serialize operations on a single execution line.

3. |linux_shell| Open the candi folder with the Windows File Explorer:

.. code-block:: text
  :class: copy-button

  explorer.exe .

.. tip::
  If the Windows File Explorer does not open, you can modify the files below directly in the Ubuntu terminal: :ref:`modify candi installation parameters with nano`.

|win_shell| Then, modify the installation parameters with Notepad (or any other text editor):

  * open the ``deal.II-toolchain/packages/p4est.package`` file. To ensure that the Lethe test suite works, deal.II must be configured with p4est 2.3.6, whereas candi installs a more recent version by default. Uncomment (remove the leading ``#``) the ``VERSION=2.3.6`` line and the ``CHECKSUM`` line that follows it, and comment (add a leading ``#``) the ``VERSION=2.8.7`` line and the ``CHECKSUM`` line that follows it. These lines should then read:

    .. code-block:: text

      VERSION=2.3.6
      CHECKSUM=4b35d9cc374e3b05cd29c552070940124f04af8f8e5e01ff046e39833de5e153

      #VERSION=2.8.7
      #CHECKSUM=0a1e912f3529999ca6d62fee335d51f24b5650b586e95a03ef39ebf73936d7f4

  * save and close

.. note::
  The ``DEAL_II_VERSION`` variable of the ``candi.cfg`` file sets the version of deal.II that is installed. Keep its default value, which is the latest deal.II release (``v9.8.0`` at the time of writing).

4. |linux_shell| Still in the candi folder, run candi installation script:

.. code-block:: text
  :class: copy-button

  ./candi.sh -j$numprocs

Where ``$numprocs`` corresponds to the number of processors used for the compilation:
  * if you have less than 8Gb of RAM, use 2 procs: ``./candi.sh -j2``
  * if you have 16Gb of RAM and above, ``$numprocs`` can be the number of physical cores minus 1. For instance, for a computer with 6 physical cores: ``./candi.sh -j5``

.. tip::

  Candi will print messages asking you if you installed the dependency. Hit ``Enter`` two times to validate and the installation will launch. If new lines are written in the console, this means the installation is going on correctly.

  If the installation is stuck (no change on the console for a few minutes), hitting ``Enter`` can unstuck it.

  You can exit the installation at any time hitting ``Ctrl+C`` 2-3 times.

5. |win_shell| At the end of the installation, check that you have deal.II and its dependencies installed: in the ``/home/<user_name>/dealii-candi`` folder, you should have a ``deal.II-v9.8.0`` folder (named after the installed deal.II version), as well as folders for the dependencies, namely: p4est, parmetis, petsc and trilinos.

6. |linux_shell| Load the environment of candi, which notably defines the ``DEAL_II_DIR`` variable used by Lethe to find deal.II, every time you open a terminal:

.. code-block:: text
  :class: copy-button

  echo "source $HOME/dealii-candi/configuration/enable.sh" >> ~/.bashrc
  source ~/.bashrc

.. note::

  Even if we use a ``echo`` command, nothing will be outputted in the terminal: the text is written directly at the end the ``.bashrc`` file, which is executed every time a terminal is opened.

.. _install-lethe-wsl:

Installing Lethe (Step #2)
-------------------------------------

1. |linux_shell| Create the folder structure and download Lethe:

.. code-block:: text
  :class: copy-button

  mkdir -p $HOME/Software/lethe; cd $HOME/Software/lethe
  git clone https://github.com/chaos-polymtl/lethe --single-branch git
  mkdir build inst

The ``lethe`` folder then contains:

* ``git`` with the source code of Lethe,
* ``build`` for the compilation files (``cmake`` and ``make`` commands),
* ``inst`` for the installed executables (``make install`` command).

2. |linux_shell| Configure Lethe:

.. code-block:: text
  :class: copy-button

  cd build
  cmake ../git -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../inst

3. |linux_shell| Compile and install Lethe:

.. code-block:: text
  :class: copy-button

  make -j$numprocs install

Where ``$numprocs`` corresponds to the number of processors used for the compilation:
  * if you have less than 8Gb of RAM, use 1 to 2 procs: ``make -j1 install`` or ``make -j2 install``
  * if you have 16Gb of RAM and above, ``$numprocs`` can be the number of physical cores minus 1. For instance, for a computer with 6 physical cores: ``make -j5 install``

4. |linux_shell| Add the Lethe executables to your ``PATH``, so that they can be launched from any folder:

.. code-block:: text
  :class: copy-button

  echo 'export PATH=$PATH:$HOME/Software/lethe/inst/bin' >> ~/.bashrc
  source ~/.bashrc

5. |linux_shell| (optional) Finally, it is recommended to test your installation. In the ``build`` folder, run:

.. code-block:: text
  :class: copy-button

  ctest -j$numprocs

This will take from a few minutes to an hour, depending on your hardware. At the end, you should have this message on the console:

  .. code-block:: text

    100% tests passed

.. warning::
  The Lethe test suite requires that deal.II be configured with p4est 2.3.6, otherwise the tests that include restart files or that use the ``lethe-fluid-vans`` or ``lethe-fluid-particles`` executables will fail. This is the case if you installed deal.II with apt, or with candi following the instructions above. Even if the tests fail, Lethe should work as expected (including the restart capabilities).

Congratulations, you are now ready to use Lethe! For instance, proceed to :doc:`../first_simulation`.

Updating deal.II and Lethe
-------------------------------------

If you have already installed deal.II and Lethe, you can update them without doing the entire installation from scratch.

Updating deal.II through apt
+++++++++++++++++++++++++++++++++

|linux_shell| As all other ``apt`` packages, run:

.. code-block:: text
  :class: copy-button

  sudo apt update
  sudo apt upgrade -y

Updating deal.II with Candi
+++++++++++++++++++++++++++++++++

|linux_shell| Set ``DEAL_II_VERSION`` to the desired deal.II version in the ``candi.cfg`` file, then run candi again:

.. code-block:: text
  :class: copy-button

  cd $HOME/Software/candi
  ./candi.sh -j$numprocs

The new version is installed in a new ``deal.II-<version>`` folder. Remove the configuration file of the previous version (``$HOME/dealii-candi/configuration/deal.II-<previous version>``) so that ``DEAL_II_DIR`` points to the new version, and open a new terminal.

Updating Lethe
+++++++++++++++++++++++++++++++++

|linux_shell| Download the latest version of Lethe, then compile and install it:

.. code-block:: text
  :class: copy-button

  cd $HOME/Software/lethe/git
  git pull
  cd ../build
  cmake ../git -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../inst
  make -j$numprocs install

If deal.II was updated, you might need to empty the ``build`` folder (``rm -rf $HOME/Software/lethe/build/*``) and run these commands again, but this is rarely the case.


Troubleshooting
-------------------------------------

.. _modify candi installation parameters with nano:

Modify Candi Installation Parameters with Nano
+++++++++++++++++++++++++++++++++++++++++++++++

|linux_shell| If the Windows File Explorer does not open, you can modify the candi parameter files in the Ubuntu terminal directly.

.. note::
  You cannot click, so use the keyboard arrows to move inside the text.

1. Open the desired file in the terminal with ``nano`` (built-in text editor):

.. code-block:: text

  cd <folder_name>
  nano <file_name>

.. admonition:: Example for the p4est.package

  .. code-block:: text

    cd $HOME/Software/candi/deal.II-toolchain/packages
    nano p4est.package

2. Modify the text in the file, using only the keyboard.

3. Save the file:

  * hit ``Ctrl + X``
  * a prompt will appear at the bottom of the terminal asking ``Save modified buffer?``
  * confirm by hitting ``y``
  * a prompt will appear at the bottom of the terminal to recall the file name
  * hit ``Enter`` to confirm
  * the file will be closed automatically and you will be back on the Ubuntu terminal
