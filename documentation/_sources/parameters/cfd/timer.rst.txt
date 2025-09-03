=====
Timer
=====

The timer section controls the frequency at which the timing of the simulation is output. Setting ``type=iteration`` outputs the timing of all elements of the simulation at the end of each iteration. Setting its value to ``end`` outputs the timing of all elements only at the end of the simulation. The timing can be disabled completely, which is necessary for the unit and functional tests.

.. code-block:: text

  subsection timer
    # Clock monitoring methods. Choices are none, iteration or end
    set type = none
  end


The timer output in Lethe is an important mechanism to monitor the relative cost of the different functions. At every time step, the timing functions will report the time spent in each function.

The following table shows an example of the output of the timer:

.. code-block:: text

    +---------------------------------------------------------+------------+------------+
    | Total wallclock time elapsed since start                |      1.11s |            |
    |                                                         |            |            |
    | Section                                     | no. calls |  wall time | % of total |
    +---------------------------------------------+-----------+------------+------------+
    | Assemble RHS                                |         9 |     0.275s |        25% |
    | Assemble matrix                             |         7 |     0.293s |        26% |
    | Calculate and output norms after Newton its |         8 |    0.0135s |       1.2% |
    | Output VTU                                  |         2 |    0.0491s |       4.4% |
    | Read mesh and manifolds                     |         1 |   0.00733s |      0.66% |
    | Setup DOFs                                  |         1 |    0.0144s |       1.3% |
    | Setup ILU                                   |         7 |    0.0801s |       7.2% |
    | Solve linear system                         |         7 |     0.369s |        33% |
    +---------------------------------------------+-----------+------------+------------+

For every block of functions, Lethe reports the number of calls, the wall time spent and the percentage of time spent in the function. It is important to monitor this when running a simulation. The most important ones from a user perspective in most of the simulations are:

* ``Assemble RHS`` refers to the assembly of the right-hand side of the equation, that is the residual of the equation.

* ``Assemble matrix`` refers to the time spent assembling the matrix associated with the equation.

* ``Output VTU`` refers to the output of the ``.vtu``, ``.pvtu`` and ``.pvd`` files.

* ``Read mesh and manifolds`` refers to the time it takes to create a ``dealii`` triangulation or read a ``gmsh`` one, and attach corresponding manifolds.

* ``Setup ILU`` describe the time spent to set-up the preconditioner. Other options are ``Setup AMG`` and ``Setup GMG``.

* ``Solve linear system`` refers to the time spent solving the linear system of equations, without the time required to assemble the preconditioner.

Depending on the type of simulation, there may be other categories that appear in the timer. This is useful to monitor the functions responsible for taking up most of the simulation time. 
