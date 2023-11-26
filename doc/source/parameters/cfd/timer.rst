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

   +---------------------------------------------+------------+------------+
   | Total wallclock time elapsed since start    |      11.5s |            |
   |                                             |            |            |
   | Section                         | no. calls |  wall time | % of total |
   +---------------------------------+-----------+------------+------------+
   | assemble_rhs                    |         4 |      2.28s |        20% |
   | assemble_system                 |         4 |      7.78s |        68% |
   | output                          |         2 |      1.21s |        11% |
   | setup_ILU                       |         4 |     0.164s |       1.4% |
   | solve_linear_system             |         4 |    0.0812s |      0.71% |
   +---------------------------------+-----------+------------+------------+

For every block of functions, Lethe reports the number of calls, the wall time spent and the percentage of time spent in the function. It is important to monitor this when running a simulation.

* ``assemble_rhs`` refers to the assembly of the right-hand side of the equation, that is the residual of the equation.

* ``assemble_system`` refers to the time spent assembling the matrix associated with the equation.

* ``output`` refers to the output of the ``.vtu``, ``.pvtu`` and ``.pvd`` files.

* ``setup_ILU`` or ``setup_AMG`` describe the time spent to set-up the preconditioner.

* ``solve_linear_system`` refers to the time spent solving the linear system of equations, without the time required to assemble the preconditioner.


Depending on the type of simulation, there may be other categories that appear in the timer. This is useful to monitor the functions responsible for taking up most of the simulation time. 
