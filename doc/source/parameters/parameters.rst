##########################
Parameters guide
##########################

Launching an application requires an executable of the required solver, and a parameters file (with extension .prm). The parameter file controls all aspects of Lethe and drives the simulation. This section aims at describing the various parameters available within Lethe.

The parameter file is composed of different subsections that have meaningful titles. In the following, the principal subsections of a parameter template file are explained. The complete options for the parameter file can be obtained by doing the following:
* Go to /lethe/build/applications/navier_stokes_parameter_template/
* Run the navier_stokes_parameter_template using the command ./navier_stokes_parameter_template.
Two parameter files are obtained for 2D and 3D respectively. These files contain all options and information regarding the parameters that can be used to run a simulation.

In the parameter file format, the parameters are established one by one using the following syntax: ``set parameter name = value``. For example, ``set time end = 10`` would fix the end time to 10. The arguments can be either doubles, integers, strings, booleans or a choice between a predefined number of variables. In the parameter files, comments are preceded by the sharp symbol (e.g. '#comment'). The parameter file format has a sanity checking mechanism in place and will throw an error if an unknown parameter or an invalid value are entered.


.. toctree::
    :maxdepth: 1
    :glob:
    :titlesonly:

    cfd/cfd
    dem/dem
    unresolved_cfd-dem/unresolved_cfd-dem
    resolved_cfd-dem/resolved_cfd-dem
