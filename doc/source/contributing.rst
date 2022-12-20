.. tip::
    Check the `utils <https://github.com/lethe-cfd/lethe-utils>`_ repository for additional cases, scripts, geo files and meshes.

############
Contributing
############

Pull requests
=============

Title
-----

Each pull request should only either fix an issue, add a single feature or restructure the code for one aspect, plus a test, if needed, to certify it. If there are many aspects to change, it is better to create multiples branches and create multiple pull requests. The title should start with an action verb and only contain one idea.

Review process
--------------

* At least 2 reviewers for minor changes, like a bug fix or a small feature
* At least 3 reviewers for bigger changes like bigger features and architectural reconfiguration.
* The reviewer should be notified before opening the pull request. The reviewers should be selected from already existing contributor to the code. (See section: Insights\Contributors).

Review responsibility
---------------------

* Give a review of the code implementation and general functionality of the code.
* Give a review on the code description and comments.
* Verify that the pull request meets the requirement described above regarding testing and the format of the pull request.
* Give the review in a timely manner.

Documentation
=============

Build
-----

Setup
^^^^^

To build the doc on your personal machine, you'll need the following requirements:

* Python v3.x
* make
* Sphinx v4.x

To install required system packages:

.. code-block:: shell
  :class: copy-button

  sudo apt-get install python3 build-essential

To install Sphinx globally:

.. code-block:: shell
  :class: copy-button

  pip install 'sphinx==4.*'

To install additional packages:

.. code-block:: shell
  :class: copy-button

  pip install sphinx-copybutton

Then, install the Sphinx requirements:

.. code-block:: shell
  :class: copy-button

  pip install -r doc/requirements.txt

Build HTML
^^^^^^^^^^

To build standalone HTML files like the CI would, enter the following commands:

.. code-block:: shell
  :class: copy-button

  cd doc
  make html

The generated files should be in the ``build/html`` directory. Open ``index.html`` in a browser to view the rendered documents.

Equations
---------

Sphinx can render equations using the MathJax backend.

Examples
^^^^^^^^

.. code-block:: RST

    .. math::
        i^2=-1

gives:

.. math::

   i^2=-1


.. code-block:: RST

    .. math::
        df=\frac{\partial f}{\partial t}

gives:

.. math::

   df=\frac{\partial f}{\partial t}

.. code-block:: RST

    .. math::
        \rho\left[\frac{\partial \bar{u}}{\partial t} + \bar{u}\cdot\bar{\nabla} \bar{u} \right] = - \bar{\nabla} \bar{p} + \mu \bar{\nabla}^2 \bar{u} + \rho \bar{g}

gives: 

.. math::
    \rho\left[\frac{\partial \bar{u}}{\partial t} + \bar{u}\cdot\bar{\nabla} \bar{u} \right] = - \bar{\nabla} \bar{p} + \mu \bar{\nabla}^2 \bar{u} + \rho \bar{g}

Code
----

Python syntax highlight: 

.. code-block:: RST

    .. code-block:: python
        your code

C++ syntax highlight: 

.. code-block:: RST

    .. code-block:: cpp
        your code

The code block in text mode is commonly used in the examples to show sections of a parameter file:

.. code-block:: RST

    .. code-block:: text
        your code

.. warning::
    Keep in mind the proper indentation of the sections of the parameter file. We recommend the use of the ``prmindent`` script located in the ``contrib/utilities`` folder before copying sections of a parameter file to the documentation. Do NOT use tabs in these blocks of code, as they will not be recognized, leading to the wrong indentation.

Examples
^^^^^^^^

.. code-block:: RST

    .. code-block:: python
        for i in range(5):
        print(i)

gives:

.. code-block:: python

    for i in range(5):
    print(i)

.. code-block:: RST

    .. code-block:: cpp
        for (int i = 0 ; i<5 ; i++) {
            std::cout << i << std::endl;
        }

.. code-block:: cpp

    for (int i = 0 ; i<5 ; i++) {
        std::cout << i << std::endl;
    }

Tables
------

.. code-block:: RST

    .. list-table::
        :header-rows: 1

        * - My
            - Beautiful
            - Table
        * - tables
            - are
            - rendered
        * - with
            - automatic
            - strip

Gives:

.. list-table::
   :header-rows: 1
   :align: center

   * - My
     - Beautiful
     - Table
   * - tables
     - are
     - rendered
   * - with
     - automatic
     - strip
