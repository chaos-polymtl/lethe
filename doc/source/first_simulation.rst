################################
Launching Your First Simulation
################################

The objective of this section is to highlight the steps that are necessary to launch your first simulation after you have cloned and compiled Lethe.

Launching an application requires an executable of the required solver, and a parameters file (with extension .prm). After building Lethe, the solver executable files can be found in : ``$BUILD_FOLDER/applications`` directory.

For instance, ``lethe-fluid-block`` executable is located in the ``$BUILD_FOLDER/applications/lethe-fluid-block`` directory.

The executable for the solvers can be used directly from the folder it is compiled to. This can be achieved by:

* Writing the absolute path of the solver (e.g. ``$BUILD_FOLDER/applications/lethe-fluid-block/lethe-fluid-block``);
* Adding the lethe folder paths to your ``PATH`` environment variable;
* Specifying an installation folder when you compile Lethe to ensure that all the applications are grouped within a single folder;
* Locally copying the executable to the folder you are running your simulation from.

All these workflows can achieve the same result.

To launch a simulation, you must specify the solver executable and the parameter file in the following format: ``solver parameter_file``. For example, ``lethe-fluid poiseuille2d.prm``

In what follows, we describe a simple procedure to launch your first simulation using Lethe.

Step 1: Copying an Example
---------------------------

The source folder of lethe contains an examples folder. This folder contains ready to run examples. Some examples use the mesh generation capacity of Lethe and only require a parameter file, whereas others contain an additional .msh file to describe the mesh. In the present case, we copy the lid driven cavity example to a new destination of your choice using the terminal:

.. code-block:: text

 cp -r $SOURCE_FOLDER/examples/incompressible_flow/2d_lid_driven_cavity destination/first_simulation

Step 2: Launching the Example
-----------------------------

The cavity example we are launching uses the *lethe-fluid* solver. All of the solvers of Lethe can be found in the build folder where you have compiled Lethe or within the installation folder. Inside of your build folder, six sub folders should be found:

* ``/applications``
* ``/applications_tests``
* ``/CMakeFiles``
* ``/prototypes``
* ``/source``
* ``/tests``

Inside the ``/applications`` folder, there is one folder for each solver of Lethe. In the ``/lethe-fluid`` folder, we find the executable file with the same name as the folder: ``lethe-fluid``. This solver solves the 2D incompressible Navier-Stokes equations using a Galerkin Least-Square formulation.

From the ``/first_simulation`` folder we have created, we can launch the simulation directly. If you have decided to copy the executable to the ``first_simulation`` folder, you can launch using the following command: ``./lethe-fluid cavity.prm``. You can also launch the simulation using the absolute path of the executable: ``$BUILD_FOLDER/applications/lethe-fluid/lethe-fluid cavity.prm``.


Step 3: Post-processing the Results
------------------------------------

Once the application has run, the simulation results can be looked at by opening the .pvd file using Paraview.

Understanding the Examples
---------------------------

Lethe comes pre-packaged with some examples which are documented on the present documentation in the :doc:`examples/examples` tab. We greatly encourage you to look at these examples to understand how Lethe can be used to solve different problems. For a more in-depth understanding of the parameter file, the reader can take a look at the :doc:`parameters/parameters` section for general and application specific parameters. Tutorial examples offer basic cases, focusing on a specific aspect, whereas Engineering Applications provide complex cases to showcase some of the features of Lethe.
