..
  SPDX-FileCopyrightText: Copyright (c) 2022-2023, 2025-2026 The Lethe Authors
  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

##########################
Parameters Guide
##########################

Launching an application requires an executable of the required solver, and a parameters file (with extension ``.prm``). The parameter file controls all aspects of Lethe and drives the simulation. The parameter file is structured with different subsections, as explained in the following: 

.. toctree::
    :maxdepth: 1
    :glob:
    :titlesonly:

    cfd/cfd
    dem/dem
    unresolved-cfd-dem/unresolved-cfd-dem
    sharp-immersed-boundary/sharp-immersed-boundary

Lethe also has a few head level parameters that are compatible with all applications:

.. code-block:: text

  set dimension        = 0
  set print parameters = none
  set comment message  =

- ``dimension``: Dimension of the simulated problem. Lethe only supports 2D (``dimension = 2``) and 3D (``dimension = 3``) simulations at the moment.

  .. attention::
    This parameter must be defined in the parameters file. Otherwise, an error message will be printed:

    .. code-block:: text

      While reading the dimension from the input file, Lethe found a value that is
      neither 2 or 3. Since August 2023, Lethe requires that the user explicitly
      specify the dimension of the problem within the parameter file. This can be
      achieved by adding set dimension = 2 or set dimension = 3 within the parameter
      file.

- ``print parameters``: Prints the parameters of the simulation on the console before starting the simulation when set to ``changed only`` or ``all``. The latter (``all``) prints the values of all the parameters declared and parsed by the application, while the former (``changed only``) prints only the values of the parameters that have been modified (non-default values). By default, ``print parameters`` is set to ``none``, therefore no parameters are printed.

- ``comment message``: Prints a message on the console before starting the simulation when specified. On the console, the message is preceded by a line with ``User comment:`` .

  .. admonition:: Example

    In the parameters file, when

    .. code-block:: text

      set comment message = This test checks that the adaptive capillary time-step constraint is well \napplied for cases with a maximum capillary time-step ratio set to 10 (Δt/Δt_σ). \nThe mesh, the densities and the surface tension considered here remain constant.

    On the console, the following is displayed:

    .. code-block:: text

      User comment:
      This test checks that the adaptive capillary time-step constraint is well
      applied for cases with a maximum capillary time-step ratio set to 10 (Δt/Δt_σ).
      The mesh, the densities and the surface tension considered here remain constant.

    .. tip::
      A few escape sequences are supported by the ``comment message``:

      * ``\n``: Newline
      * ``\t``: Horizontal tab
      * ``\r``: Carriage return
      * ``\b``: Backspace
      * ``\v``: Vertical tab

.. tip:: 
	The complete options for the parameters file, with information regarding the parameters that can be used, can be obtained by doing the following:

	* Go to ``/lethe/build/applications/navier_stokes_parameter_template/``
	* Run the ``navier_stokes_parameter_template`` using the command ``./navier_stokes_parameter_template``

	Two parameter files, with default values for all the parameters, are obtained for 2D and 3D respectively.

.. note::
	* The parameters are established one by one using the following syntax: ``set parameter name = value``
	* Comments are preceded by the sharp symbol (e.g. ``#comment``)
	* The parameter file format has a sanity checking mechanism in place and will throw an error if an unknown parameter or an invalid value is entered.

.. warning::
  The radioactive particle-tracking (RPT) applications of Lethe has been migrated to a separate repository which is available `here <https://github.com/chaos-polymtl/lethe-rpt>`_.