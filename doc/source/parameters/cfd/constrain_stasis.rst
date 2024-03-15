=================
Constrain Stasis
=================

This subsection is used to define temperature-dependent solid domains within a defined fluid.
Homogenous constraints are applied on velocity and pressure degrees of freedom (DoFs) of cells found within the prescribed temperature range to mimic a solid.

.. note::
  Pressure DoFs of solid cells that are adjacent to fluid cells (i.e., cells with temperature values outside the specified range) are left unconstrained in order to adequately solve fluid flows.

.. attention::
    In order to use this feature, the ``heat transfer`` physic must be enabled (``true``) in the :doc:`./multiphysics` subsection.

The subsection with default parameters goes as follows:

.. code-block:: text

    subsection constrain stasis
      set enable                = false
      set number of constraints = 0
      subsection constraint 0
        set fluid id                 = 0
        set phase fraction tolerance = 1e-4
        set min temperature          = -999.0
        set max temperature          = 0.0
      end
    end

* The ``enable`` parameter is set to ``true`` when at least one temperature-dependent stasis constraint should be applied.

* The ``number of constraints`` parameter is an integer representing the number of constraints that will be applied. It is used in multiphase (VOF) simulations to apply different constraints to each fluid. Only one constraint per fluid can be imposed. Each constraint comes with its own subsection (starting with number ``0``) containing its own set of parameters as detailed below.

  * The ``fluid id`` parameter is an integer representing the fluid on which the current constraint should be applied.

  * The ``phase fraction tolerance`` parameter is an absolute tolerance on the filtered phase fraction :math:`(\phi')` used in conjunction with VOF simulations (see :doc:`volume_of_fluid` subsection) to select the cells on which the constraint is applied. For example, if a ``phase fraction tolerance`` of :math:`10^{-4}` is specified for a constraint on a fluid of ``fluid id = 1``, cells with :math:`\phi' \in [0.9999,1.0001]` are considered as the cells of interest. Note that if the phase fraction filter is disabled, the tolerance will be applied on the solved phase fraction field :math:`(\phi)`.

  * The ``min temperature`` parameter is a double representing the minimum temperature value for a cell to be considered a solid :math:`[\Theta]`.

  * The ``max temperature`` parameter is a double representing the maximum temperature value for a cell to be considered a solid :math:`[\Theta]`.

