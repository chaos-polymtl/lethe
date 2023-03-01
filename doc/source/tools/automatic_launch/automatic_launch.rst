======================================================================
How to automatically create and launch Lethe simulations
======================================================================

.. seealso::
	All files used in this example are available in the `lethe-utils github <https://github.com/lethe-cfd/lethe-utils/>`_ in the ``cases`` folder under ``automatic_launch``.

-------------------------------------
Generate automatically multiple cases
-------------------------------------
Lets say that you are simulating a flow around a cylinder and you want to see how the inlet velocity impacts the force around the sphere.
Lasy as we are, we want to automatically generate multiple copies of the cylinder case, but change the parameter file so the boundary condition 
imposing the inlet velocity is different for each case.

You will need all these files from the `2D flow around a cylinder example <https://lethe-cfd.github.io/lethe/examples/incompressible-flow/2d-flow-around-cylinder/2d-flow-around-cylinder.html>`_:
- ``cylinder.prm``
- ``cylinder-structured.geo``

Here are the boundary conditions of the flow around the cylinder.

.. image:: images/geometry-bc.png
    :alt: The boundary conditions
    :align: center
    :name: geometry_bc


In the ``.prm`` file, we need to change the function expression of ``bc = 1``, which represents the velocity :math:`u` in the :math:`x` direction.

To do so, we use the `Jinja2 <https://jinja.palletsprojects.com/en/3.1.x/>`_ python package.
It allows us to create a parameter file template, where some parameter variables, in double brackets ``{{}}``, can be replace to any value we want.

The boundary conditions section in the ``.prm`` file becomes as follow.

.. code-block:: text
    
    subsection boundary conditions
      set number = 3
      subsection bc 0
        set type = noslip
      end
      subsection bc 1
        set type = function
        subsection u
          set Function expression = {{velocity_x}}
        end
        subsection v
          set Function expression = 0
        end
        subsection w
          set Function expression = 0
        end
      end
      subsection bc 2
        set type = slip
      end
    end

.. note::
	Note the ``{{velocity_x}}`` parameter variable in the Jinja2 format.
	This will allow us to insert a specified value for a case of the 2D around a cylinder.

The 

We will now present how to generate multiple folders, containing different parameters files, to ultimately launch them as separated cases.
We can generate these folders locally and even on a cluster.

""""""""""""""""""""""""""""""""""
Locally
""""""""""""""""""""""""""""""""""
Before importing the packages, we need to install them using ``pip`` and ``requirements.txt`` (available in the lethe-utils folder).

.. code-block:: text
    
	pip install -r requirements.txt

Then, we can import the right packages to launch the script.

.. code-block:: python
    
	from jinja2 import Template
	import os
	import numpy as np
	import shutil

The first thing to do, is set up the constants of the script.

.. code-block:: python
    
	PATH = os.getcwd()
	PATH_PREFIX = 'cylinder_u_'
	PRM_FILE = 'cylinder.prm'
	MESH_FILE = 'cylinder-structured.msh'

The ``PATH`` is the current path of the user directory.
The ``PATH_PREFIX`` will specify how we want to name each folder.
The ``PRM_FILE`` is the name of the parameter file of the Lethe simulation.
The ``MESH_FILE`` is the name of the mesh used for the simulations.

.. warning::
	The ``.msh`` file is not available as it is. You will need to run ``gmsh`` in order to generate the mesh around the cylinder from the ``.geo`` file.

Then we specify the range of velocity we want to explore.
In this example, we will generate 20 cases of the flow around a cylinder, where the inlet velocity varies from 1 to 10 :math:`m/s`.

.. code-block:: python
    
	number_of_cases = 20
	first_velocity = 1
	last_velocity = 10
	velocity = np.linspace(1, 10, number_of_cases)

Now, the fun begins.

For each velocity in the range specified above,

.. code-block:: python
    
	for u in velocity:

we will:

1. Open the parameter file (with the ``{{}}`` variables) and read it.

.. code-block:: python

	fic_prm = open(PRM_FILE, 'r')
	read_prm = fic_prm.read()

2. Create a Jinja2 ``Template``.
   
.. code-block:: python

	template = Template(read_prm)

4. Render the template with the right value and close the reading state of the parameter file.

.. code-block:: python

	parameters = template.render(velocity_x = u)

	fic_prm.close()

.. warning::
	In the render step, it is really important to use the same variable name as the template file.

Then, we will need to copy in the ``case_path`` (containing the folder of one case) all the files we need for the simulation.

1. Name the ``case_path``.
   
.. code-block:: python

	case_folder_name = f'{PATH_PREFIX}{u:.2f}'
	case_path = f'{PATH}/{case_folder_name}'

2. Create the ``case_path``.

.. code-block:: python

	os.mkdir(case_path)

3. Copy the ``.prm`` file and the ``.msh`` file from the current ``PATH`` to the ``case_path``.

.. code-block:: python

	shutil.copy(f'{PATH}/{PRM_FILE}', f'{case_path}/{PRM_FILE}')
	shutil.copy(f'{PATH}/{MESH_FILE}', f'{case_path}/{MESH_FILE}')

1. Enter the ``case_path`` and write the parameter file with the rendered template.

.. code-block:: python

	os.chdir(case_path)

	write_prm = open(PRM_FILE, 'w')
	write_prm.write(parameters)
	write_prm.close()

1. Never forget to step back from the directory, in order to create another template and another folder for the next case.

.. code-block:: python

	os.chdir('../')

And voilà! The final current directory should look like this:

.. code-block:: text

	+---automatic_launch
	|   +---cylinder_u_1.00
	|   |       cylinder-structured.msh
	|   |       cylinder.prm
	|   |
	|   +---cylinder_u_1.95
	|   |       cylinder-structured.msh
	|   |       cylinder.prm
	|   |
	|   +---cylinder_u_2.42
	|   |       cylinder-structured.msh
	|   |       cylinder.prm
	|   |
	|   \---cylinder_u_10.00
	|   |       cylinder-structured.msh
	|   |       cylinder.prm

""""""""""""""""""""""""""""""""""""""
On Digital Alliance of Canada clusters
""""""""""""""""""""""""""""""""""""""
If you want to generate the same folders, but on a cluster, the same script applies.

The main differences are:

1. Specify the shell script that will launch a job on the cluster.

.. code-block:: python

	SHELL_FILE = 'launch_lethe.sh'

2. Copy the ``.sh`` from the current ``PATH`` to the ``case_path``.

.. code-block:: python

	shutil.copy(f'{PATH}/{SHELL_FILE}', f'{case_path}/{SHELL_FILE}')

This last step will allow to launch one job for each case.

-----------------------------------
Launch automatically multiple cases
-----------------------------------

""""""""""""""""""""""""""""""""""
Locally
""""""""""""""""""""""""""""""""""


""""""""""""""""""""""""""""""""""""""
On Digital Alliance of Canada clusters
""""""""""""""""""""""""""""""""""""""
