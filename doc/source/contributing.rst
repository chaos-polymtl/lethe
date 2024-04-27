.. tip::
    Check the `utils <https://github.com/chaos-polymtl/lethe-utils>`_ repository for additional cases, scripts, geo files and meshes.

############
Contributing
############

Pull Requests
=============

Title
-----

Each pull request should only either fix an issue, add a single feature or restructure the code for one aspect, plus a test, if needed, to certify it. If there are many aspects to change, it is better to create multiples branches and create multiple pull requests. The title should start with an action verb and only contain one idea.

Format
------

Lethe uses clang-format to have a uniform indentation across all source files. You should always run the indentation script in the contrib folder before creating a PR: ``contrib/utilities/indent-all``. If there are compatibility issues, you can run before that script either ``./contrib/utilities/download_clang_format`` or  ``./contrib/utilities/compile_clang_format``. For parameter files, there is a specific script that can be run as follows:  ``prmindent -i name_of_file.prm``.


Review Process
--------------

* At least 2 reviewers for minor changes, like a bug fix or a small feature
* At least 3 reviewers for bigger changes like bigger features and architectural reconfiguration.
* The reviewer should be notified before opening the pull request. The reviewers should be selected from already existing contributor to the code. (See section: Insights\Contributors).

Reviewers Responsibility
------------------------

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

  pip install 'sphinx==7.*'

To install additional packages:

.. code-block:: shell
  :class: copy-button

  pip install sphinx-copybutton

and:

.. code-block:: shell
  :class: copy-button

  sudo apt-get install graphviz

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

General Rules and Format
------------------------

Throughout the documentation, you may notice that the different pages follow a certain format to ensure uniformity and help users to navigate more fluidly. Here are the different elements that must be considered when contributing to the documentation of Lethe:

- Titles and subtitles must be capitalized following the `Chicago Manual of Style <https://www.chicagomanualofstyle.org/book/ed17/frontmatter/toc.html>`_:

  - First and last words of the title/subtitle must be capitalized.
  - Capitalize:

    - adjectives (e.g., Small, Large)
    - adverbs (e.g., Warmly, Rapidly)
    - nouns (e.g., Ball, Cylinder)
    - pronouns (e.g., They, She, He)
    - subordinating conjunction when fewer than 5 letters (e.g., When, Once)
    - verbs (e.g., Melt, Create)

  - **Do not** capitalize:

    - articles (e.g., a, an, the)
    - coordinating conjunctions (e.g., and, but, for)
    - prepositions (e.g., at, by, to)
    - second word after a hyphenated prefix (e.g., Mid-, Non-) in compound modifiers (e.g., Mid-year, Non-linear)
    - words with less than four letters

  .. tip::

    If you are unsure of the capitalization of your title or subtitle, you can use `online tools <https://capitalizemytitle.com/style/Chicago/>`_ to help you out.

- Examples generally contain the following subsections in the listed order:

  - **Features** lists features of the example.
  - **Files Used in This Examples** lists files used in the example in alphabetical order.
  - **Description of the Case** describes the system studied in the example.
  - **Parameter File** describes the different parameter subsections involved in the example. Each subsection of interest of the parameter file begins with its name as a sub-heading.
  - **Running the Simulation** displays the command used to run the example and gives an scale of the running duration.
  - **Results** or **Results and Discussion** displays results of the simulation and comments on them.
  - **Possibility for Extension** lists different interesting ways to extend the example.
  - **References** lists references used in the example in IEEE referencing style format.

    .. seealso::

      More information on the IEEE referencing style can be found in the `IEEE Reference Guide <https://ieeeauthorcenter.ieee.org/wp-content/uploads/IEEE-Reference-Guide.pdf>`_.

      Here is an example of how references should appear in:

      - In-text citation:

        The *Lethe: An open-source parallel high-order adaptative CFD solver for incompressible flows* article by Blais *et al.* `[1] <https://doi.org/10.1016/j.softx.2020.100579>`_ is used as an example.

      - The ``References`` list:

        `[1] <https://doi.org/10.1016/j.softx.2020.100579>`_ B. Blais *et al.*, “Lethe: An open-source parallel high-order adaptative CFD solver for incompressible flows,” *SoftwareX*, vol. 12, p. 100579, Jul. 2020, doi: 10.1016/j.softx.2020.100579.

        Following the format:

        .. container::

          [#] A. A. Author, "Name of the paper," *Abbreviated Title of the Journal*, vol. x, no. x, pp. xxx-xxx, Abbreviated month, year, doi: xxx.

        .. Important::

          When citing, the "*et al.*" notation is used in:

          - In-text citation if there are three or more authors for a given reference.
          - The ``References`` list if there are more than six authors for a given reference.

  .. note::

    These subheadings can take the singular or plural form depending on the example. Feel free to adapt them and add more layers to structure your own examples.

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
