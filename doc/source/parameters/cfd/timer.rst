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
    | Total wallclock time elapsed since start                |      16.1s |            |
    |                                                         |            |            |
    | Section                                     | no. calls |  wall time | % of total |
    +---------------------------------------------+-----------+------------+------------+
    | Assemble RHS                                |         9 |      5.72s |        36% |
    | Assemble matrix                             |         7 |      7.59s |        47% |
    | Calculate and output norms after Newton its |         8 |     0.318s |         2% |
    | Distribute constraints after linear solve   |         7 |   0.00887s |         0% |
    | Output VTU                                  |         2 |     0.765s |       4.8% |
    | Read mesh and manifolds                     |         1 |     0.296s |       1.8% |
    | Set initial conditions                      |         1 |     0.243s |       1.5% |
    | Setup DOFs                                  |         1 |     0.727s |       4.5% |
    | Setup ILU                                   |         7 |     0.104s |      0.65% |
    | Solve linear system                         |         7 |     0.267s |       1.7% |
    +---------------------------------------------+-----------+------------+------------+

For every block of functions, Lethe reports the number of calls, the wall time spent and the percentage of time spent in the function. It is important to monitor this when running a simulation. The most important ones from a user perspective in most of the simulations are:

* ``Assemble RHS`` refers to the assembly of the right-hand side of the equation, that is the residual of the equation.

* ``Assemble matrix`` refers to the time spent assembling the matrix associated with the equation.

* ``Output VTU`` refers to the output of the ``.vtu``, ``.pvtu`` and ``.pvd`` files.

* ``Read mesh and manifolds`` refers to the time it takes to create a ``dealii`` triangulation or read a ``gmsh`` one, and attach corresponding manifolds.

* ``Setup ILU`` describe the time spent to set-up the preconditioner. Other options are ``Setup AMG`` and ``Setup GMG``.

* ``Solve linear system`` refers to the time spent solving the linear system of equations, without the time required to assemble the preconditioner.

Depending on the type of simulation, there may be other categories that appear in the timer. This is useful to monitor the functions responsible for taking up most of the simulation time. 
