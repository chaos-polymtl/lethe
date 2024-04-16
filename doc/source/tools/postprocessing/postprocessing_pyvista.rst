====================================
Post-processing Lethe with PyVista
====================================

Lethe has a `post-processing module <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_ written in Python based on Lethe users specific needs. `lethe_pyvista_tools <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_ was conceived to optimize the reading and post-treatment of Lethe data using Python.

The module is powered by `PyVista <https://docs.pyvista.org/>`_, a powerful 3D plotting and mesh analysis tool. It is a *pythonic* interface to deal with Visualization Toolkit (VTK) data. The module also uses the powerful `multiprocessing <https://docs.python.org/3/library/multiprocessing.html>`_ Python library to increase post-processing speed by running tasks in multiple processors.



.. warning::

  For `lethe_pyvista_tools <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_ to work, along with `Python 3 <https://www.python.org/downloads/>`_, the following libraries are needed: `os <https://docs.python.org/3/library/os.html>`_, `NumPy <https://numpy.org/>`_, `PyVista <https://docs.pyvista.org/>`_, `tqdm <https://tqdm.github.io/>`_, `matplotlib <https://matplotlib.org/stable/index.html>`_, `SciPy <https://scipy.org/>`_, and `scikit-learn <https://scikit-learn.org/stable/index.html>`_. If any of the modules are missing, use `pip <https://pypi.org/project/pip/>`_ to install it by running ``pip3 install $NAME_OF_THE_MODULE`` on the terminal. Alternatively, use the ``requirements.txt`` file and install them all at once running ``pip install -r requirements.txt`` located at ``contrib/postprocessing/``.


------------------------------
Importing lethe_pyvista_tools
------------------------------

There are two ways to have access to ``lethe_pyvista_tools``: you can install the module and have it on your machine just like ``os``, ``sys``, or other PyPi modules; or, you can use it without installing. Here, we explain how to do both.


Installing
~~~~~~~~~~

Installing ``lethe_pyvista_tools`` is very simple. In the `contrib/postprocessing <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_ of your Lethe version you will find all source files, including ``setup.py`` and ``requirements.txt``. Those two files take care of configuring your installation and assuring that everything is set correctly. To install it, navigate to `contrib/postprocessing <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_ on your terminal and run:



.. code-block::

    pip3 install .

or

.. code-block::

    pip install .

.. note::

    If you do not have ``pip`` on your machine, run ``sudo apt install python3-pip``.

.. tip::

    Uninstalling ``lethe_pyvista_tools`` is exactly like doing so with any other library: ``pip3 uninstall lethe_pyvista_tools`` or ``pip3 uninstall lethe_pyvista_tools``.


Importing without installing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case, we use the `sys <https://docs.python.org/3/library/sys.html>`_ module to import it directly from Lethe's directory into a python session:


.. code-block::

  import sys
  path_to_module = '$LETHE_PATH/contrib/postprocessing/'
  sys.path.append(path_to_module)
  from lethe_pyvista_tools import *
  import matplotlib.pyplot as plt

where `sys.path <https://docs.python.org/3/library/sys.html#sys.path:~:text=in%20version%203.10.-,sys.path%C2%B6,-A%20list%20of>`_ is a list of strings that specifies the search path for modules. Another option is to simply copy the `lethe_pyvista_tools.py <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_ file to the same folder as the Python post-processing script, and import it as follows:


.. code-block::

  from lethe_pyvista_tools import *

One third and very convenient way to always import the module without copying it or even adding the ``sys.path.append(path_to_module)`` is permanently adding the path to the module (i.e., ``/contrib/postprocessing/``) to your `PYTHONPATH <https://docs.python.org/3/library/sys_path_init.html#:~:text=The%20PYTHONPATH%20environment%20variable%20is,all%20installed%20Python%20versions%2Fenvironments.>`_.


The ``*`` means that we want to import all members of lethe_pyvista_tools.



Using lethe_pyvista_tools
~~~~~~~~~~~~~~~~~~~~~~~~~~

To get quick-started, follow the hand-on :doc:`Small Scale Rotating Drum Post-processing example <../../examples/dem/small-scale-rotating-drum-postprocessing/small-scale-rotating-drum-postprocessing>`. It has a detailed explanation on how to use the module. You can also start with a "raw" template file, check `this example file <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing/example.py>`_.

