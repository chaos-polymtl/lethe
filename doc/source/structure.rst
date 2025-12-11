##########################
Folder Structure of Lethe
##########################

The folder structure of Lethe is relatively complex. Once you have cloned and compiled Lethe, it is normal that you might be confused as to the meaning of the various folders in the cloned directory (which we will refer to as the lethe folder) or in the build folder. This section of the documentation aims at demystifying these elements.

================
The Lethe Folder
================
In this folder, the following structure is present:

* ``/applications``
* ``/applications_tests``
* ``/contrib``
* ``/doc``
* ``/examples``
* ``/include``
* ``/logo``
* ``/prototypes``
* ``/source`` 
* ``/tests``

In what follows, we explain the content of each folder and the logic behind this complex structure. We explain the folder not in alphabetical order but in a more "logical" order.

Applications and Their Tests
----------------------------

Lethe is designed to contain a number of solvers for single or multiphysics problems. These various solvers take the form of multiple executables which are named applications. They are stored in the ``/applications`` folder. Each application is housed in its own separate folder. For example, the folder ``/applications/lethe-fluid`` contains the source which instantiates the GLS Navier-Stokes solver.

In Lethe, solvers are minimal applications with only a main file. They contain little source code (generally < 100 lines). Their only goal is to instantiate the solver classes with the appropriate dimension of the problem. Since all of Lethe solvers are templated for the dimension of the problem (int dim), the same source code is re-used in both 2D and 3D. Applications are made to be comprehensive. They are fully controlled from text files (.prm and eventually .json).

The ``/applications_tests`` folder contains the functional tests for the applications. These tests are grouped in folders that correspond to the name of the application. For example, the folder ``/applications_tests/lethe-fluid`` folder contains multiple tests for the GLS Navier-Stokes solver. Each test is two or more files and is identified with a single prefix. It must necessarily contain a ".prm" file, which is the parameter file that is used to launch the test. It must also contain a ".output" file, which is the expected output of the solver. Additional files may be present if the test is to be run with multiple processors combination or if the test requires a mesh file.

For example, the *cylinder_gls* test is made of three files: ``cylinder_gls.output``, which is the expected output, ``cylinder_gls.prm``, which is the parameter file used to control the simulation (with the appropriate boundary conditions, etc.) and the ``cylinder_structured.msh`` file which is a GMSH mesh that is required for the simulation to run.

You may notice that some outputs are preceded by a mpirun=2 prefix. For example, the ``mms-2d_gls.mpirun=2.output`` file. This additional prefix indicates that the test will be run with two processors. If this prefix is not present, the test will always be run with a single processor. To ensure enhanced compatibility with the CI, tests should not use over 3 processors.


The ``include`` and the ``source`` Folders
------------------------------------------

The ``/include`` folder houses the prototypes for all the functions and the classes within Lethe whereas the ``/source`` folder houses the source code. All the files within the include folder are .h files whereas all files within the source folder are .cc, as per the deal.II standards.

It can readily be noticed that the include and source folder are subdivided into additional sub-directories. At the time of this writing the sub-directories are :

* ``/core``
* ``/dem``
* ``/fem-dem``
* ``/solvers``

The content of Lethe is subdivided into modules that are generally independent of one another. This encapsulates furthermore the content of the software. The core module contains the core functionality of Lethe, whereas the solvers module contains the applications that are related to the finite element solution of multiphysics problems. The dependency between the module can only be unidirectional. For example, the solvers module requires the core module, but the core module could not require the solvers module. 


The ``tests`` Folder
--------------------

The ``/tests`` folder houses the unit tests for all of the modules of Lethe. The unit tests are different from the application's tests. Their scope is only to verify that a single class or functions behave in the expected manner. They are used to test individual components of the Lethe and not complete applications. Hence, these tests are grouped by modules.

The ``prototypes`` Folder
-------------------------

The ``/prototypes`` folder houses temporary applications in Lethe. These applications are prototypes that are very close to the steps of deal.II. They are used to test novel development before they are fully integrated within Lethe. The prototypes are not expected to be stable.

The prototypes folder contains a template folder which contains a virgin prototype. Each prototype should be implemented within their own folder. When implementing a new prototype, the following should be done:

* Copy the template folder to a new folder within the prototypes folder that bears the name of what you want to develop
* Within the new folder, modify the ``CMakeLists.txt`` file to ensure that the name of the application has the name you want it to have. This is achieved by modifying the ``SET(TARGET "protype_template")`` line to ``SET(TARGET "my_project_name")``
* Modify the name of the .cc files so that it corresponds to the name of your prototype. All .cc and .h files will be detected automatically by the ``CMakeLists.txt`` GLOB instructions
* Add the new folder to the list of prototypes to be taken into account by modifying the ``prototypes/CMakeLists.txt`` and adding a ``ADD_SUBDIRECTORY(my_prototype)``


The ``contrib`` Folder
----------------------

The ``/contrib`` folder houses various scripts that help for the maintenance of Lethe. It notably contains the ``contrib/utilities/indent-all`` script which can be used to automatically indent the entire Lethe source code to ensure that it adheres to the clang-format rules contained in the ``/contrib`` folder.

The ``doc`` Folder
------------------

The ``/doc`` folder contains the source files of this documentation page of Lethe. To contribute or compile the documentation on your own machine follow the instructions of the :doc:`contributing/documentation` tab.

The ``examples`` Folder
-----------------------

The ``/examples`` folder includes the parameter file and the post-processing scripts of examples using different applications in Lethe. It is subdivided into additional sub-directories, namely:

* ``/cfd_dem``
* ``/dem``
* ``/incompressible_flow``
* ``/multiphysics``

Detailed descriptions of most of these examples can be found on the :doc:`examples/examples` tab of this page.

The ``logo`` Folder
-------------------

The ``/logo`` folder is the simplest one, it houses the logo of Lethe in various formats :)!

====================
The ``build`` Folder
====================

This folder is obtained after compiling Lethe and it contains all the relevant executables for all applications and tests available. 



