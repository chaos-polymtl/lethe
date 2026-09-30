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

  This tutorial is aimed at Windows users who have no prior knowledge of Linux. It explains how to install Ubuntu 24.04 LTS in the Windows Subsystem for Linux (WSL). Lethe is then installed following the :doc:`regular_installation` instructions.

Throughout this tutorial:
  * |win_shell| indicates operations performed in the Windows session, and
  * |linux_shell| indicates operations performed in the Linux subsystem.

.. tip::
  To execute a command on a shell (Ubuntu or Windows command prompt), type or copy/paste the given command and hit ``Enter``. Multiple commands are given in multiple lines, or separated by ``;``: when copying/pasting, they will be executed one after the other.

Installing WSL and Ubuntu
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

Installing deal.II and Lethe
------------------------------------

From this point on, Ubuntu running in WSL behaves like any Linux installation: follow the :doc:`regular_installation` instructions, running every command in the Ubuntu terminal (|linux_shell|). Installing deal.II with apt (:ref:`install-deal.II-apt`) is by far the easiest way to proceed.

.. tip::
  A fresh Ubuntu installation does not include the tools needed to compile and test Lethe. Install them first:

  .. code-block:: text
    :class: copy-button

    sudo apt install build-essential cmake git numdiff

  To edit a file of the Ubuntu file system (e.g. the ``candi.cfg`` file of candi), open its folder with ``explorer.exe .`` and use any Windows text editor, or edit it directly in the terminal with ``nano <file_name>``.
