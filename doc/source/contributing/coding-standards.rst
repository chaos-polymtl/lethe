================
Coding Standards
================

Lethe follows the `coding conventions of deal.II <https://dealii.org/developer/doxygen/deal.II/CodingConventions.html>`_, with the additions described on this page. Formatting is not discussed here, since it is handled automatically by clang-format (see :doc:`pull-requests`). These rules apply to all new and modified code and are verified during the review. When in doubt, imitate the code that surrounds your change.

Variable and Function Names
---------------------------

* Variables, functions and files use ``snake_case``.
* Names are complete words that describe the physical or numerical meaning of the quantity: ``kinematic_viscosity`` or ``pressure_scaling_tolerance``, not ``nu`` or ``tol``. Short mathematical names are only used for local indices that follow the deal.II conventions, such as ``q`` for the quadrature points and ``i`` and ``j`` for the shape functions in assembly loops.
* Quantities that count something start with ``n_``: ``n_q_points``, ``n_dofs``.
* A variable that stores a parameter mirrors the name of the parameter in the parameter file. For example, ``set number of fluids`` is stored in ``number_of_fluids``. This makes it easy to find where a parameter is used in the code.
* Booleans read as a statement: ``is_``, ``has_`` or ``enable_`` variables, ``density_is_constant()`` or ``is_non_newtonian()`` functions.
* Functions start with a verb that describes what they do: ``setup_dofs``, ``assemble_system_matrix``, ``get_`` and ``set_`` accessors.
* Member variables have no prefix or suffix (no ``m_`` or trailing ``_``). In class templates, members of a base class are accessed through ``this->``.
* Constants have a name, are declared ``constexpr`` or ``const`` close to where they are used and also use ``snake_case``. Do not use magic numbers in the code. ``ALL_CAPS`` is reserved for preprocessor macros.

.. code-block:: cpp

  // Avoid
  unsigned int nf = prm.get_integer("number of fluids");
  for (unsigned int k = 0; k < nf; ++k)

  // Prefer
  number_of_fluids = prm.get_integer("number of fluids");
  for (unsigned int i_fluid = 0; i_fluid < number_of_fluids; ++i_fluid)

Class and Type Names
--------------------

* Classes, structures and type aliases use ``UpperCamelCase`` without underscores: ``SimulationControl``, ``FluidDynamicsMatrixBased``, ``GlobalVectorType``. Acronyms keep their capital letters: ``CFDDEMSolver``, ``GLSNavierStokesAssemblerCore``, ``CLSAssemblerBDF``.
* An abstract base class ends with ``Base`` and the name of a derived class states what it specializes. For example, ``CahnHilliardFilterBase`` is implemented by ``CahnHilliardTanhFilter`` and ``CahnHilliardClipFilter``.
* A file is named after the main class it contains, in ``snake_case``. ``FluidDynamicsMatrixBased`` is declared in ``include/solvers/fluid_dynamics_matrix_based.h`` and implemented in ``source/solvers/fluid_dynamics_matrix_based.cc``. The header guard is ``lethe_<file_name>_h``, here ``lethe_fluid_dynamics_matrix_based_h``.
* Namespaces use ``UpperCamelCase``: ``Parameters``, ``BoundaryConditions``, ``DEM``.
* A choice between options is an ``enum class`` with an ``UpperCamelCase`` name and ``snake_case`` values, such as ``Parameters::SimulationControl::TimeSteppingMethod::bdf2``. Do not encode options with integers, booleans or strings.
* Template parameters follow deal.II: ``dim`` and ``spacedim`` for the dimensions, ``UpperCamelCase`` for types (``VectorType``, ``Number``).
* Exceptions declared with the ``DeclException`` macros use ``UpperCamelCase`` and are named after the error they report: ``FileDoesNotExist``, ``NumberOfFluidsError``.

Defensive Programming
---------------------

Lethe is used by many people who did not write the code they run, and simulations often run for hours or days on clusters. A simulation that stops immediately with a clear error message costs its user a few minutes. A simulation that silently runs with a wrong value can produce no results at all or results that look plausible but are wrong, which costs much more. Lethe must therefore be written defensively: every assumption made by a piece of code is either enforced by the compiler or checked explicitly, and errors are reported as early as possible, with a message that explains how to fix them.

In practice:

* **Validate the parameters when they are read.** Restrict the range of an entry with its pattern in ``declare_parameters`` (e.g., ``Patterns::Double(0.)`` for a quantity that cannot be negative) and check the consistency between parameters with ``AssertThrow`` in ``parse_parameters``. Errors are then reported before the mesh is even read.
* **Check the preconditions of functions**: sizes with ``AssertDimension``, indices with ``AssertIndexRange`` and values that must be finite with ``AssertIsFinite``.
* **Handle every case.** A ``switch`` or an ``if``/``else if`` chain on an ``enum class`` handles every value. A combination of features that is not supported throws an exception instead of silently doing nothing.
* **Guard divisions** by quantities that can vanish (cell measure, density, time step, norm of a gradient) and give tolerances a name. Never compare floating-point numbers with ``==``.
* **Let the compiler help you.** Declare ``const`` every variable that does not change after its initialization, pass large objects by ``const &`` and initialize every variable where it is declared.
* **Write useful error messages.** A message states what is wrong, the offending value and how to fix it. Exceptions raised in more than one place are declared in ``include/core/exceptions.h``: reuse them and add yours there. For example:

.. code-block:: cpp

  DeclException2(
    FileDoesNotExist,
    std::string,
    std::string,
    << "The " << arg2 << " <" << arg1 << "> does not exist. "
    << "Verify that its name is spelled correctly. Relative paths are "
    << "interpreted from the directory in which the application was launched, "
    << "not from the directory of the parameter file.");

``Assert`` or ``AssertThrow``
-----------------------------

deal.II provides two families of assertion macros, which differ in when they are checked:

* ``Assert`` (as well as ``AssertDimension``, ``AssertIndexRange`` and ``AssertIsFinite``) is only checked in Debug mode. In Release mode, it is removed by the compiler and costs nothing.
* ``AssertThrow`` is checked in both Debug and Release modes, at the cost of a test every time it is evaluated.

Choose between them according to where the check is located:

* **In loops and other performance-critical sections** (loops over cells, faces, quadrature points, degrees of freedom or particles, matrix-free operator evaluations, contact detection), use ``Assert``. Since it is free in Release mode, use it generously to check sizes, indices and internal invariants.
* **Outside performance-critical sections** (constructors, setup functions, parameter parsing, reading of meshes and files, checks done once per time step or per linear solve), use ``AssertThrow`` for robustness. Its cost is negligible there, and the check then also protects Release builds, which are the ones used for production runs. This also applies to internal consistency checks.
* **Errors that a user can trigger** (a value in the parameter file, a mesh, a restart file, an unsupported combination of features, a file that cannot be read) are always checked with ``AssertThrow``. Do this check once, outside the loops, during parsing or setup: an ``Assert`` located in a loop catches nothing in a Release build.

.. list-table::
   :header-rows: 1
   :align: center

   * - Check
     - Location
     - Macro
   * - Entry of the parameter file out of range or inconsistent with another entry
     - ``parse_parameters``
     - ``AssertThrow``
   * - Boundary condition or feature not supported by a solver
     - Setup
     - ``AssertThrow``
   * - Size of a solution vector once the degrees of freedom are distributed
     - ``setup_dofs``
     - ``AssertThrow``
   * - Index of a shape function or of a quadrature point
     - Assembly loop
     - ``AssertIndexRange``
   * - Size of a local matrix or vector
     - Assembly loop
     - ``AssertDimension``
   * - Value computed at a quadrature point is finite
     - Assembly loop
     - ``AssertIsFinite``

.. code-block:: cpp

  // Parameter parsing: a user error, checked once and in every build
  AssertThrow(number_of_fluids <= max_fluids,
              NumberOfFluidsError(number_of_fluids));

  // Assembly: an internal invariant in a loop, checked in Debug mode only
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          AssertIndexRange(i, local_rhs.size());
          // ...
        }
    }

.. warning::

  An ``AssertThrow`` must fail on all the processes at the same time. Checks on parameters already do, since every process reads the same parameter file. A condition that depends on locally owned data can, however, fail on a single process. That process stops while the others wait for it in the next collective communication, and the simulation can hang until it reaches its wall time. Reduce such a condition over all the processes before checking it:

  .. code-block:: cpp

    const bool negative_density_found =
      Utilities::MPI::logical_or(local_negative_density_found,
                                 mpi_communicator);
    AssertThrow(!negative_density_found,
                ExcMessage("A negative density was found. Verify the "
                           "physical properties of the fluids."));

Smart Pointers
--------------

Memory in Lethe is managed with smart pointers, so that every object has a clearly identified owner and is released automatically:

* Do not use ``new`` and ``delete``. Create objects with ``std::make_unique`` or ``std::make_shared``. A raw pointer, when one is needed to interface with an external library, never owns the object it points to.
* Use ``std::unique_ptr`` when an object has a single owner. This is the default choice, including for objects that are only created when they are needed:

  .. code-block:: cpp

    std::unique_ptr<FEValues<dim>> fe_values_temperature;
    if (needs_temperature)
      fe_values_temperature = std::make_unique<FEValues<dim>>(
        *this->mapping,
        dof_handler_temperature->get_fe(),
        temperature_quadrature,
        update_values | update_quadrature_points);

* Use ``std::shared_ptr`` only when the ownership is really shared, that is when several objects hold the same object and none of them is guaranteed to outlive the others. In Lethe, this is the case of the triangulation, the mapping and the ``SimulationControl``, which are shared by the solver, its auxiliary physics and the physical property models:

  .. code-block:: cpp

    simulation_control = std::make_shared<SimulationControlTransient>(
      simulation_parameters.simulation_control);

* Only pass a smart pointer to a function that shares the ownership of the object, such as a constructor that stores it. Lethe passes them as ``const std::shared_ptr<T> &``. A function that only uses the object takes it by reference (``const T &``).
* Avoid cycles of ``std::shared_ptr`` (two objects that own each other): neither of them is ever destroyed.
