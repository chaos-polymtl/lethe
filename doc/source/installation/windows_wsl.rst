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

  This tutorial is aimed at Windows users who have no prior knowledge of Linux. To keep it simple, all dependencies are installed using candi. Installation options given in this tutorial are well suited for lethe users. If you are a developer or need more options, see :doc:`regular_installation`.

Throughout this tutorial:
  * |win_shell| indicates operations performed in the Windows session, and
  * |linux_shell| indicates operations performed in the Linux subsystem.

.. tip::
  To execute a command on a shell (Ubuntu or Windows command prompt), type or copy/paste the given command and hit ``Enter``. Multiple commands are given in multiple lines, or separated by ``;``: when copying/pasting, they will be executed one after the other.

Installing WSL and Ubuntu (Step #0)
------------------------------------

1. |win_shell| Install WSL (Windows Subsystem for Linux). Open PowerShell or Windows Command Prompt in administrator mode by right-clicking and selecting "Run as administrator", enter the wsl --install command, then restart your machine.

.. code-block:: text
  :class: copy-button

  wsl --install

2. |win_shell| Enable WSL and Ubuntu 22.04 LTS from the microsoft store, following the steps indicated `in this tutorial <https://linuxconfig.org/ubuntu-22-04-on-wsl-windows-subsystem-for-linux>`_

.. admonition:: Verify the installed version of WSL

  In the windows command prompt (Start menu > ``cmd``):

  .. code-block:: text
      :class: copy-button

      wsl -l -v

  should indicate ``version 2``. If not, follow this to update WSL: https://docs.microsoft.com/fr-fr/windows/wsl/install#upgrade-version-from-wsl-1-to-wsl-2

3. |win_shell| Launch Ubuntu (from the start menu) and |linux_shell| update Ubuntu: 

.. code-block:: text
  :class: copy-button

  sudo apt update
  sudo apt upgrade

.. tip::
  The ``sudo`` command will ask you to type your user password, as defined during Ubuntu installation. Note that Linux does not show any symbol while typing a password, contrary to Windows with ``*``: simply type your password and press ``Enter``.

When prompted "do you want to continue?", proceed by typing ``y`` and hitting ``Enter``.

4. |win_shell| (optional) To activate copy/paste in the Linux sub-terminal (`tutorial with screenshots <https://defragged.org/2020/10/29/how-to-copy-paste-in-windows-subsystem-for-linux-wsl/>`_):

  * right-click on the Ubuntu Window pane header
  * in ``Properties``, select ``Use Ctrl+Shift+C/V as Copy/Paste``
  * you can then use ``Ctrl+Shift+V`` to paste text or commands in the Linux sub-terminal

5. |win_shell| (optional) For better ease in the Linux terminal (better coloring, multiple tabs), change the default terminal:

  * in the microsoft store, download ``Windows Terminal``
  * in the ``parameters`` of ``Windows Terminal``, select on the left pannel "start": change default profile with ``Ubuntu-22.04``
  * from now on, you can use this application instead to launch Ubuntu terminal

.. tip::
  A (very) few Linux commands useful for navigation:
    * ``mkdir $dir``: (make directory) create a directory with the name specified as ``$dir``
    * ``cd $dir``: (change directory) move to the directory ``$dir``
    * ``cd ..``: move up to the parent directory
    * ``pwd``: (print working directory) return the directory you are in
    * ``cd $HOME``: move to your home directory (``/home/<user_name>/``)

  You can find `here <https://linuxconfig.org/linux-commands>`_ a thorough guide for the most basic Linux commands.


The following step is to install deal.II. This can be done through

  1. Advanced Packaging Tool (apt) (**this is by far the easiest way to proceed**) : :ref:`install-deal.II apt` (recommended for users)

  2. Candi shell script (`candi github page <https://github.com/dealii/candi>`_): :ref:`install-deal.II candi` (recommended for developers)

.. _install-deal.II apt:

Installing deal.II using apt (Step #1)
-----------------------------------------

This is done following `this procedure <https://www.dealii.org/download.html#:~:text=page%20for%20details.-,Linux%20distributions,-Arch%20Linux>`_.

1. |linux_shell| In case you are using Ubuntu, you will need to `update the backports <https://launchpad.net/~ginggs/+archive/ubuntu/deal.ii-9.5.1-backports>`_:

.. code-block:: text
  :class: copy-button

  sudo add-apt-repository ppa:ginggs/deal.ii-9.5.1-backports
  sudo apt update

2. |linux_shell| To install deal.II, run:

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

.. _install-deal.II candi:

Installing deal.II using Candi (Step #1)
-----------------------------------------

.. important::
  This step is by far the most troublesome in all Lethe installation. Read and follow each step carefully.

1. |linux_shell| Install candi required packages:

.. code-block:: text
  :class: copy-button

  sudo apt-get install lsb-release git subversion wget \
  bc libgmp-dev build-essential autoconf automake cmake \
  libtool gfortran libboost-all-dev zlib1g-dev openmpi-bin \
  openmpi-common libopenmpi-dev libblas3 libblas-dev \
  liblapack3 liblapack-dev libsuitesparse-dev

.. tip::
  The symbols ``\`` indicate that this a single command written on multiple lines.

2. |linux_shell| Install compilers:

.. code-block:: text
  :class: copy-button

  sudo apt-get install gcc-10 g++-10 gfortran-10

.. admonition:: Check the default version of the compilers

  In the Ubuntu terminal:

  .. code-block:: text

    gcc --version

  should return ``gcc (Ubuntu 10.X.X...) 10.X.X`` or higher. If not, go to :ref:`change compiler default version`.

3. |linux_shell| Create folders (suggested structure):

.. code-block:: text
  :class: copy-button

  mkdir Software; cd Software
  mkdir candi; cd candi

Note the use of ``;`` which enable to serialize operations on a single execution line.

4. |linux_shell| Download candi:

.. code-block:: text
  :class: copy-button

  git clone https://github.com/dealii/candi.git .

Do not forget the ``.`` at the end of the command, which means "here".

5. |win_shell| Modify installation parameters (deal.II version and trilinos version):

  * open Windows file manager, and on the left panel (along with ``Files``, ``Computer`` etc.) click on the ``Ubuntu`` mount.

  .. tip::
    If you do not see any ``Ubuntu`` mount, use this alternative method: :ref:`modify candi installation parameters with nano`.

  * navigate to reach the candi folder, in: ``/home/<user_name>/Software/candi``
  * open the ``candi.cfg`` file with notepad (or other text editor) and change the following lines:

  +--------+------------------------------------------+----------------------------------------+
  | line # | initial parameter                        | changed parameter                      |
  +========+==========================================+========================================+
  |      7 | ``CLEAN_BUILD=OFF``                      | ``CLEAN_BUILD=ON``                     |
  +--------+------------------------------------------+----------------------------------------+
  |     86 | ``# PACKAGES="${PACKAGES} once:netcdf"`` | ``PACKAGES="${PACKAGES} once:netcdf"`` |
  +--------+------------------------------------------+----------------------------------------+
  |     97 | ``DEAL_II_VERSION=v9.5.0``               | ``DEAL_II_VERSION=master``             |
  +--------+------------------------------------------+----------------------------------------+

  * save and close
  * navigate to reach the following subfolder: ``deal.II-toolchain/packages/``
  * open the ``trilinos.package`` file with notepad and change the following lines:

  .. tip::
    The prefix ``#`` is used to comment a line. Here we are simply commenting lines 44 and 45, and uncommenting lines 50 and 51, to change the trilinos version.

  +--------+------------------------------------------------+-----------------------------------------------+
  | line # | initial parameter                              | changed parameter                             |
  +========+================================================+===============================================+
  |     44 | ``VERSION=12-18-1``                            | ``#VERSION=12-18-1``                          |
  +--------+------------------------------------------------+-----------------------------------------------+
  |     45 | ``CHECKSUM=9c1d151169949bca6cf203831e4d6aee``  | ``#CHECKSUM=9c1d151169949bca6cf203831e4d6aee``|
  +--------+------------------------------------------------+-----------------------------------------------+
  |     50 | ``#VERSION=12-12-1``                           | ``VERSION=12-12-1``                           |
  +--------+------------------------------------------------+-----------------------------------------------+
  |     51 | ``#CHECKSUM=ecd4606fa332212433c98bf950a69cc7`` | ``CHECKSUM=ecd4606fa332212433c98bf950a69cc7`` |
  +--------+------------------------------------------------+-----------------------------------------------+

  * save and close
  * still in the subfolder ``deal.II-toolchain/packages/``, open the ``p4est.package`` file with notepad and change the following lines:

  .. tip::
    The prefix ``#`` is used to comment a line. Here we are simply uncommenting line 7, and commenting lines 9 to 12, to change the p4est version.

  +--------+------------------------------------------------+-----------------------------------------------+
  | line # | initial parameter                              | changed parameter                             |
  +========+================================================+===============================================+
  |     7  | ``#VERSION=2.2;CHECKSUM=6943949a...``          | ``VERSION=2.2;CHECKSUM=6943949a...``          |
  +--------+------------------------------------------------+-----------------------------------------------+
  |     9  | ``VERSION=2.3.2``                              | ``#VERSION=2.3.2``                            |
  +--------+------------------------------------------------+-----------------------------------------------+
  |     10 | ``CHECKSUM=076df9e...``                        | ``#CHECKSUM=076df9e...``                      |
  +--------+------------------------------------------------+-----------------------------------------------+
  |     11 | ``CHECKSUM="${CHECKSUM} b41c8ef29ca...``       | ``#CHECKSUM="${CHECKSUM} b41c8ef29ca...``     |
  +--------+------------------------------------------------+-----------------------------------------------+
  |     12 | ``CHECKSUM="${CHECKSUM} 0ea6e4806b6...``       | ``#CHECKSUM="${CHECKSUM} 0ea6e4806b6...``     |
  +--------+------------------------------------------------+-----------------------------------------------+

  * save and close

6. |linux_shell| Still in the candi subfolder, run candi installation script:

.. code-block:: text
  :class: copy-button

  ./candi.sh -j$numprocs

Where ``$numprocs`` corresponds to the number of processors used for the compilation:
  * if you have less than 8Gb of RAM, use 1 to 2 procs: ``./candi.sh -j1`` or ``./candi.sh -j2``
  * if you have 16Gb of RAM and above, ``$numprocs`` can be the number of physical cores minus 1. For instance, for a computer with 6 physical cores: ``./candi.sh -j5``

.. tip::

  Candi will print messages asking you if you installed the dependency. Hit ``Enter`` two times to validate and the installation will launch. If new lines are written in the console, this means the installation is going on correctly. The installation will take from 1 to 3 hours depending on your hardware.

  If the installation is stuck (no change on the console for a few minutes), hitting ``Enter`` can unstuck it.

  You can exit the installation at any time hitting ``Ctrl+C`` 2-3 times.

7. |win_shell| At the end of the installation, check that you have deal.II and its dependencies installed:

  * on Windows file manager, go to the Ubuntu mount
  * in ``home/<user_name>`` you should have a folder ``deal.ii-candi``, or ``dealii-candi``
  * inside this folder, you should have folders for the dependencies, namely: p4est, petsc, parmetis, trilinos
  * you should also see this folder: ``deal.II-master``

8. |linux_shell| Add a deal.II environment variable in Ubuntu through the following command:

.. code-block:: text
  :class: copy-button

  echo "export DEAL_II_DIR=$HOME/dealii-candi/deal.II-master" >> ~/.bashrc

.. note::

  Even if we use a ``echo`` command, nothing will be outputted in the terminal: the text is written directly at the end the ``.bashrc`` file.

.. warning::

  For this change to be effective, you may need to restart your Ubuntu terminal.


Installing Lethe (Step #2)
-------------------------------------

1. |linux_shell| Set-up the folder structure. Create the ``Software`` folder (if you are doing the candi installation, this folder should alredy exist from Step #1).

.. code-block:: text
  :class: copy-button

  mkdir Software; cd Software 

In the ``Software`` folder created (if you are in the candi folder, type ``cd ..``), type:

.. code-block:: text
  :class: copy-button

  mkdir -p lethe/{git,build,inst}

After installation is complete, the folder structure will be:

* ``lethe/git`` with lethe downloaded files (git),
* ``lethe/build`` for compilation files (``cmake`` command),
* ``lethe/inst`` for installation files (``make install`` command).

2. |linux_shell| Download lethe:

.. code-block:: text
  :class: copy-button

  cd lethe
  git clone https://github.com/chaos-polymtl/lethe git

3. |linux_shell| Build lethe:

.. code-block:: text
  :class: copy-button

  cd build
  cmake ../git -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../inst/

4. |linux_shell| Compile lethe:

.. code-block:: text
  :class: copy-button

  make -j$numprocs

Where ``$numprocs`` corresponds to the number of processors used for the compilation:
  * if you have less than 8Gb of RAM, use 1 to 2 procs: ``make -j1`` or ``make -j2``
  * if you have 16Gb of RAM and above, ``$numprocs`` can be the number of physical cores minus 1. For instance, for a computer with 6 physical cores: ``make -j5``

5. |linux_shell| (optional) Finally, it is recommended to test your installation:


Run the tests in the build folder:

.. code-block:: text
  :class: copy-button

  ctest -j$numprocs

This will take from a few minutes to an hour, depending on your hardware. At the end, you should have this message on the console:

  .. code-block:: text

    100% tests passed

.. note:: If you are running these tests for the first time, install ``numdiff`` (if you need superuser privilege, use sudo):

  .. code-block:: text
    :class: copy-button
    
    apt-get numdiff

  or

  .. code-block:: text
    :class: copy-button
    
    apt install numdiff

.. warning:: 
  The lethe test suites requires that deal.II be configured with p4est 2.2.1, otherwise the test which include restart files will fail.

Congratulations, you are now ready to use lethe! For instance, proceed to :doc:`../first_simulation`.

Updating deal.II and Lethe
-------------------------------------

If you have already installed deal.II and lethe, you can update them without doing the entire installation from scratch:

Through apt
+++++++++++++++++++++++++++++++++

1. |linux_shell| As all other ``apt`` packages, run:

.. code-block:: text
  :class: copy-button

  sudo apt update
  sudo apt upgrade -y

2. |linux_shell| Then, update lethe:

.. code-block:: text
  :class: copy-button

  cd ../lethe/git
  git pull
  cd ../build
  cmake ../git -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../inst/
  make -j$numprocs

With Candi
+++++++++++++++++++++++++++++++++
1. |linux_shell| Update deal.ii by typing, from your home directory:

.. code-block:: text
  :class: copy-button

  cd Software/candi
  ./candi.sh -j$numprocs

2. |linux_shell| Then, update lethe:

.. code-block:: text
  :class: copy-button

  cd ../lethe/git
  git pull
  cd ../build
  cmake ../git -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../inst/
  make -j$numprocs


Troubleshooting
-------------------------------------

.. _change compiler default version:

Change Compiler Default Version
+++++++++++++++++++++++++++++++++++++

|linux_shell| After you installed ``gcc-10``, ``g++-10`` and ``gfortran-10``, manually update default versions in the terminal:

.. code-block:: text

  sudo update-alternatives --remove-all gcc
  sudo update-alternatives --remove-all g++
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 10
  sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 10
  sudo update-alternatives --set cc /usr/bin/gcc
  sudo update-alternatives --set c++ /usr/bin/g++

Then, check again the version used:

.. code-block:: text

  gcc --version

Should return ``gcc (Ubuntu 10.X.X...) 10.X.X``.


.. _modify candi installation parameters with nano:

Modify Candi Installation Parameters with Nano
+++++++++++++++++++++++++++++++++++++++++++++++

|linux_shell| If you do not see the Ubuntu mount in the Windows file manager, you can modify the candi parameter files in the Ubuntu terminal directly. 

.. note::
  You cannot click, so use the keyboard arrows to move inside the text.

1. Open the desired file in the terminal with ``nano`` (built-in text editor):

.. code-block:: text

  cd <folder_name>
  nano <file_name>

.. admonition:: Example for the candi.cfg

  .. code-block:: text

    cd /home/<user_name>/Software/candi
    nano candi.cfg

2. Modify the text in the file, using only the keyboard. 

3. Save the file: 

  * hit ``Ctrl + X``
  * a prompt will appear at the bottom of the terminal asking ``Save modified buffer?``
  * confirm by hitting ``y``
  * a prompt will appear at the bottom of the terminal to recall the file name
  * hit ``Enter`` to confirm
  * the file will be closed automatically and you will be back on the Ubuntu terminal


