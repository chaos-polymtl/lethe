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
    rpt/rpt

.. tip:: 
	The complete options for the parameter file, with information regarding the parameters that can be used, can be obtained by doing the following:

	* Go to ``/lethe/build/applications/navier_stokes_parameter_template/``
	* Run the ``navier_stokes_parameter_template`` using the command ``./navier_stokes_parameter_template``

	Two parameter files, with default values for all the parameters, are obtained for 2D and 3D respectively.

.. note::
	* The parameters are established one by one using the following syntax: ``set parameter name = value``
	* Comments are preceded by the sharp symbol (e.g. ``#comment``)
	* The parameter file format has a sanity checking mechanism in place and will throw an error if an unknown parameter or an invalid value is entered.
