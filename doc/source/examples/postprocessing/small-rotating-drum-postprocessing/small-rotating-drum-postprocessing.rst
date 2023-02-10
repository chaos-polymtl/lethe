
==========================================
Small scale rotating drum
==========================================

This is an example of how to post-process results obtained in the `Small scale rotating drum example <../../dem/rotating-drum/small-scale-rotating-drum.html>`_ using `lethe_pyvista_tools <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_, a Python module is based on `PyVista <https://docs.pyvista.org/>`_, built to facilitate the reading of Lethe results using `Python <https://www.python.org/>`_. 

.. important::
  This example uses data from the `Small scale rotating drum example <../../dem/rotating-drum/small-scale-rotating-drum.html>`_.

.. warning::
  For `lethe_pyvista_tools <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_ to work, along with `Python 3 <https://www.python.org/downloads/>`_, the following libraries are needed: `os <https://docs.python.org/3/library/os.html>`_, `NumPy <https://numpy.org/>`_, `PyVista <https://docs.pyvista.org/>`_, and `tqdm <https://tqdm.github.io/>`_. If any of the modules are missing, use `pip <https://pypi.org/project/pip/>`_ to install it running ``pip3 install $NAME_OF_THE_MODULE`` on the terminal.

Features
----------------------------------
- DEM simulation
- Post-processing using `Python <https://www.python.org/>`_, `PyVista <https://docs.pyvista.org/>`_, `lethe_pyvista_tools <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_, and `ParaView <https://www.paraview.org/>`_.


Files used in this example
----------------------------

- Parameters file for particle insertion: ``/examples/postprocessing/small-scale-rotating-drum-postprocessing/packing-rotating-drum.prm``
- Parameters file for drum rotation: ``/examples/postprocessing/small-scale-rotating-drum-postprocessing/small-rotating-drum-dem.prm``
- Python module for Lethe data post-processing: ``/contrib/postprocessing/lethe_pyvista_tools.py``
- Python script using module for rotating drum post-processing: ``/examples/postprocessing/small-scale-rotating-drum-postprocessing/example_small_rotating_drum.py``



Description of the case
-----------------------

This example post-processes Lethe-DEM data using `Python <https://www.python.org/>`_, `PyVista <https://docs.pyvista.org/>`_, `lethe_pyvista_tools <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_, and `ParaView <https://www.paraview.org/>`_. First, we follow the `Small scale rotating drum example <../../dem/rotating-drum/small-scale-rotating-drum.html>`_ to obtain the DEM files.

.. note::
  It is not necessary to use all mentioned tools, but they are used in this example to show different possibilities of data access according to user's need.


Python code
---------------

Module importing
~~~~~~~~~~~~~~~~~

The module `lethe_pyvista_tools <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_ was conceived to optimize the reading and post-treatment of Lethe data using Python. It is based on `PyVista <https://docs.pyvista.org/>`_, a versatile module that can be used to manipulate Lethe results.

First of all, we import the module to our Python script. There are two ways to do so. In the case of this example, we use the `sys <https://docs.python.org/3/library/sys.html>`_ module to import it directly from Lethe's directory:

.. code-block::

  import sys
  path_to_module = '../../../contrib/postprocessing/'
  sys.path.append(path_to_module)
  from lethe_pyvista_tools import *

where `sys.path <https://docs.python.org/3/library/sys.html#sys.path:~:text=in%20version%203.10.-,sys.path%C2%B6,-A%20list%20of>`_ is a list of strings that specifies the search path for modules. However, one can simply copy the `lethe_pyvista_tools.py <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_ file to the same folder as the Python post-processing script is and simply import it, such as:
 
.. code-block::

  from lethe_pyvista_tools import *

The ``*`` means that we want to import all members of lethe_pyvista_tools. 

Constructing the object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following step is to create an object to receive the data. In the present case, the object is called ``particles``.

.. code-block::

  particles = lethe_pyvista_tools(case_path = ".", prm_file_name = "small-rotating-drum-dem.prm")

Here, the constructor ``lethe_pyvista_tools`` receives two arguments: ``case_path`` and ``prm_file_name``. In the above code line, ``"."`` means that the ``case_path`` is the path where we currently are, but it can be any path where the case is. The ``prm_file_name`` argument takes ``"small-rotating-drum-dem.prm"``.

.. tip::
  Together with the object ``particles``, ``lethe_pyvista_tools`` creates a dictionary with all parameters in the ``.prm`` file. To access the parameter, we can use ``particles.prm_dict['$NAME_OF_THE_PARAMETER']``. In the present case for example, the diameter of the particles can be easily printed using ``print(particles.prm_dict['diameter'])``. This can be useful for post-processing routines with multiple simulations.

Reading Lethe results to Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next, we read the results into Python:

.. code-block::
  
  particles.read_lethe_to_pyvista(pvd_name = "out.pvd", first = 40)

``read_lethe_to_pyvista`` takes the name of the ``.pvd`` file generated by the simulation. The method can take 3 other arguments: ``first``, ``last``, and ``interval``, standing for the first and last time-steps to be read and the interval between the time-steps, respectively. By default, ``first = 0``, ``interval = 1``, and ``last`` is the last time-step of the data.

In the present case, since we want to use data after the full packing of particles only (after the 40th time-step), the parameter ``first`` is set to ``40``.

The ``read_lethe_to_pyvista`` reading function assigns the datasets of each time-step to the object ``particles``. Each time-step corresponds to a PyVista dataset, . adds each time-step dataset to the object ``particles``. To access the datase


Results
---------

The following movie displays the rolling regime inside the rotating drum obtained with a rotational velocity of 1 rad/s. 

.. raw:: html

  <iframe width="560" height="315" src="https://www.youtube.com/embed/qxO4MD_zg2w" title="Rotating drum - mixing study" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>


Possibilities for extension
----------------------------

- Use two types of particles with different radius to prove the Brazil-Nut effect.
- Perform an unresolved CFD-DEM simulation for wet granular flows to see the impact of the hydrodynamics of the fluid over the particles dynamics.


 
