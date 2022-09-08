####################
Windows with WSL
####################

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

1. |win_shell| Install WSL (Windows Subsystem for Linux), and Ubuntu 22.04 LTS from the microsoft store, following the steps indicated `in this tutorial <https://linuxconfig.org/ubuntu-22-04-on-wsl-windows-subsystem-for-linux>`_

.. admonition:: Verify the installed version of WSL

	In the windows command prompt (Start menu > ``cmd``):

	.. code-block:: text

		wsl -l -v

	should indicate ``version 2``. If not, follow this to update WSL: https://docs.microsoft.com/fr-fr/windows/wsl/install#upgrade-version-from-wsl-1-to-wsl-2

2. |win_shell| Launch Ubuntu (from the start menu) and |linux_shell| update Ubuntu: 

.. code-block:: text

	sudo apt update
	sudo apt upgrade

.. tip::
	The ``sudo`` command will ask you to type your user password, as defined during Ubuntu installation. Note that Linux does not show any symbol while typing a password, contrary to Windows with ``*``: simply type your password and press ``Enter``.

When prompted "do you want to continue?", proceed by typing ``y`` and hitting ``Enter``.

3. |win_shell| (optional) To activate copy/paste in the Linux sub-terminal (`tutorial with screenshots <https://defragged.org/2020/10/29/how-to-copy-paste-in-windows-subsystem-for-linux-wsl/>`_):
	* right-click on the Ubuntu Window pane header 
	* in ``Properties``, select ``Use Ctrl+Shift+C/V as Copy/Paste``
	* you can then use ``Ctrl+Shift+V`` to paste text or commands in the Linux sub-terminal

4. |win_shell| (optional) For better ease in the Linux terminal (better coloring, multiple tabs), change the default terminal:
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


Installing deal.II using candi (Step #1)
-----------------------------------------

.. important::
	This step is by far the most troublesome in all Lethe installation. Read and follow each step carefully.

1. |linux_shell| Install candi required packages:

.. code-block:: text

	sudo apt-get install lsb-release git subversion wget \
	bc libgmp-dev build-essential autoconf automake cmake \
	libtool gfortran libboost-all-dev zlib1g-dev openmpi-bin \
	openmpi-common libopenmpi-dev libblas3 libblas-dev \
	liblapack3 liblapack-dev libsuitesparse-dev

.. tip::
	The symbols ``\`` indicate that this a single command written on multiple lines.

2. |linux_shell| Install compilers:

.. code-block:: text

	sudo apt-get install gcc-10 g++-10 gfortran-10

3. |linux_shell| Create folders (suggested structure):

.. code-block:: text

	mkdir Softwares; cd Softwares
	mkdir candi; cd candi

Note the use of ``;`` which enable to serialize operations on a single execution line.

4. |linux_shell| Download candi:

.. code-block:: text

	git clone https://github.com/dealii/candi.git .

Do not forget the ``.`` at the end of the command, which means "here".

5. |win_shell| Modify installation parameters (deal.II version and trilinos version):
	* open Windows file manager, and on the left panel (along with ``Files``, ``Computer`` etc.) click on the ``Ubuntu`` mount.
	* navigate to reach the candi folder, in: ``/home/<user_name>/Softwares/candi``
	* open the ``candi.cfg`` file with notepad (or other text editor) and change the following lines:

	+--------+------------------------------------------+----------------------------------------+
	| line # | initial parameter                        | changed parameter                      |
	+========+==========================================+========================================+
	|      7 | ``CLEAN_BUILD=OFF``                      | ``CLEAN_BUILD=ON``                     |
	+--------+------------------------------------------+----------------------------------------+
	|     86 | ``# PACKAGES="${PACKAGES} once:netcdf"`` | ``PACKAGES="${PACKAGES} once:netcdf"`` |
	+--------+------------------------------------------+----------------------------------------+
	|     97 | ``DEAL_II_VERSION=v9.4.0``               | ``DEAL_II_VERSION=master``             |
	+--------+------------------------------------------+----------------------------------------+

	* save and close 
	* navigate to reach the following subfolder: ``deal.II-toolchain/packages/``
	* open the ``trilinos.package`` file with notepad and change the following lines:

	.. tip::
		The prefix ``#`` is used to comment a line. Here we are simply commenting lines 19 and 20, and uncommenting lines 25 and 26, to change the trilinos version.

	+--------+------------------------------------------------+-----------------------------------------------+
	| line # | initial parameter                              | changed parameter                             |
	+========+================================================+===============================================+
	|     19 | ``VERSION=12-18-1``                            | ``#VERSION=12-18-1``                          |
	+--------+------------------------------------------------+-----------------------------------------------+
	|     20 | ``CHECKSUM=9c1d151169949bca6cf203831e4d6aee``  | ``#CHECKSUM=9c1d151169949bca6cf203831e4d6aee``|
	+--------+------------------------------------------------+-----------------------------------------------+
	|     25 | ``#VERSION=12-12-1``                           | ``VERSION=12-12-1``                           |
	+--------+------------------------------------------------+-----------------------------------------------+
	|     26 | ``#CHECKSUM=ecd4606fa332212433c98bf950a69cc7`` | ``CHECKSUM=ecd4606fa332212433c98bf950a69cc7`` |
	+--------+------------------------------------------------+-----------------------------------------------+

	* save and close 

6. |linux_shell| Still in the candi subfolder, run candi installation script:

.. code-block:: text

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
	* in ``home/<user_name>`` you should have a folder ``deal.ii-candi``
	* inside this folder, you should have folders for the dependencies, namely: p4est, petsc, parmetis, trilinos
	* you should also see this folder: ``deal.II-master``

8. |linux_shell| Add a deal.II environment variable in Ubuntu through the following command:

.. code-block:: text

	echo "export DEAL_II_DIR=$HOME/deal.ii-candi/deal.II-master" >> ~/.bashrc


Installing Lethe (Step #2)
-------------------------------------

1. |linux_shell| Set-up the folder structure: in the ``Softwares`` folder created at the beginning of Step #1 (if you are in the candi folder, type ``cd ..``), type:

.. code-block:: text

	mkdir -p lethe/{git,build,inst}

After installation is complete, the folder structure will be:

* ``lethe/git`` with lethe downloaded files (git),
* ``lethe/build`` for compilation files (``cmake`` command),
* ``lethe/inst`` for installation files (``make install`` command).

2. |linux_shell| Download lethe:

.. code-block:: text

	cd lethe
	git clone https://github.com/lethe-cfd/lethe git

3. |linux_shell| Build lethe:

.. code-block:: text

	cd build
	cmake ../git -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=../inst/

4. |linux_shell| Compile lethe:

.. code-block:: text

	make -j$numprocs

Where ``$numprocs`` corresponds to the number of processors used for the compilation:
	* if you have less than 8Gb of RAM, use 1 to 2 procs: ``make -j1`` or ``make -j2``
	* if you have 16Gb of RAM and above, ``$numprocs`` can be the number of physical cores minus 1. For instance, for a computer with 6 physical cores: ``make -j5``

5. |linux_shell| (optional) Test your installation, still in the build folder:

.. code-block:: text

	ctest -j$numprocs

This will take from a few minutes to an hour, depending on your hardware. At the end, you should have this message on the console:

.. code-block:: text

	100% tests passed

Congratulations, you are ready to use lethe! You are now ready for :doc:`../first_simulation`.

Updating deal.II and lethe
-------------------------------------

If you have already installed deal.II and lethe, you can update them without doing the entire installation from scratch:

1. |linux_shell| Update deal.ii by typing, from your home directory:

.. code-block:: text

	cd Softwares/candi
	./candi.sh -j$numprocs

2. |linux_shell| Then, update lethe:

.. code-block:: text

	cd ../lethe/git
	git pull
	cd ../build
	cmake ../git -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=../inst/
	make -j$numprocs
	

