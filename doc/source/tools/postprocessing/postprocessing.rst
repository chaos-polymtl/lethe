====================================
Post-processing Lethe with PyVista
====================================

Lethe has a `post-processing module <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_ written in Python based on Lethe users specific needs.

The module is powered by `PyVista <https://docs.pyvista.org/>`_, a powerful 3D plotting and mesh analysis tool. It is a "pythonic" interface to deal wih Visualization Toolkit (VTK) data.

To get quick-started, follow the `Small Scale Rotating Drum Post-processing example <../../examples/dem/small-scale-rotating-drum-post-processing>`_.

.. warning::

  For `lethe_pyvista_tools <https://github.com/lethe-cfd/lethe/tree/master/contrib/postprocessing>`_ to work, along with `Python 3 <https://www.python.org/downloads/>`_, the following libraries are needed: `os <https://docs.python.org/3/library/os.html>`_, `NumPy <https://numpy.org/>`_, `PyVista <https://docs.pyvista.org/>`_, `tqdm <https://tqdm.github.io/>`_, `matplotlib <https://matplotlib.org/stable/index.html>`_, and `SciPy <https://scipy.org/>`_, and `scikit-learn <https://scikit-learn.org/stable/index.html>`_. If any of the modules are missing, use `pip <https://pypi.org/project/pip/>`_ to install it running ``pip3 install $NAME_OF_THE_MODULE`` on the terminal.