=====================
Parameter File Checks
=====================

The ``lethe-parameter-check`` application validates a parameter file without
running a simulation. It can also describe all the parameters of a solver
family in the JSON format, for tools such as editors or graphical interfaces
that need to know the parameters of Lethe without reimplementing them.

The option is given before the other arguments:

.. code-block:: shell

   lethe-parameter-check input_file solver_family
   lethe-parameter-check --schema solver_family dimension

where ``solver_family`` is ``lethe-fluid`` (which also covers
``lethe-fluid-block``, ``lethe-fluid-nitsche`` and ``lethe-fluid-sharp``),
``lethe-particles`` or ``lethe-fluid-particles``.

---------------------------
Check a Parameter File
---------------------------

.. code-block:: shell

   lethe-parameter-check cylinder.prm lethe-fluid

The parameter file is parsed exactly as the solvers of the family parse it. The
application exits with an error if the file is not valid, and it reports the
use of deprecated parameters or parameter values. Files read during the
parsing are given relative to the parameter file, so the check should be run
from the folder of the parameter file.

-------------------------------------
Print the Schema of the Parameters
-------------------------------------

.. code-block:: shell

   lethe-parameter-check --schema lethe-fluid 3 > lethe-fluid-3d.json

The schema lists every parameter declared by a solver family for a dimension,
in the order of their declaration. It is generated from the declarations
themselves, so it is always up to date with the version of Lethe that prints it.
Its top-level fields are:

- ``lethe_version`` and ``dealii_version``: the versions of Lethe and deal.II
  the schema was printed with.
- ``solver_family`` and ``dimension``: the solver family and the dimension the
  parameters were declared for.
- ``parameters``: the subsections and the parameters, nested as in a parameter
  file and named exactly as in a parameter file. Each parameter is described by
  its ``value`` and ``default_value`` (the same in the schema), its
  ``documentation`` and its ``pattern_description``, which gives the type of
  the parameter and, for a selection, the list of its choices (e.g.
  ``[Selection gmsh|dealii|lethe ]``).

For example, the ``type`` parameter of the ``mesh`` subsection is described as:

.. code-block:: text

   "mesh": {
       "type": {
           "value": "dealii",
           "default_value": "dealii",
           "documentation": "Type of mesh. Choices are <gmsh|dealii|lethe>.",
           "pattern": "...",
           "pattern_description": "[Selection gmsh|dealii|lethe ]"
       },
       ...
   }

The subsections whose number is set by the parameter file, such as
``bc 0``, ``bc 1``, ... in ``boundary conditions`` and ``manifold 0``,
``manifold 1``, ... in ``manifolds``, appear once, as ``bc 0`` and
``manifold 0``. They are the templates of all the subsections of their kind.

The schema only describes the parameters one by one. Whether a parameter is used
can depend on the value of another parameter (e.g. the ``time step`` is not
used by a steady simulation), which only the solvers know. Checking a parameter
file remains the way to validate a combination of parameters.
