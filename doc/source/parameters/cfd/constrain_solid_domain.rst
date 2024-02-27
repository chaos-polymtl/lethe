=======================
Constrain Solid Domain
=======================

This subsection is used to define temperature-dependant solid domains within a defined fluid.
Homogenous constraints are applied on velocity and pressure degrees of freedom of cells found within the prescribed temperature range to mimic a solid.

.. attention::
    In order to use this feature, the ``heat transfer`` physic must be enabled (``true``) in the :doc:`./multiphysics` subsection.

The subsection with default parameters goes as follows:

.. code-block:: text

    subsection constrain solid domain
      set enable                = false
      set number of constraints = 1
      subsection constraint 0
        set fluid id        = 0
        set min temperature = -999.0
        set max temperature = 0.0
      end
    end

* The ``enable`` parameter is set to ``true`` when at least one temperature-dependant solid domain constraint should be applied.

* The ``number of constraints`` parameter is an integer representing the number of constraints that will be applied. It is used in multiphase (VOF) simulations to apply different constraints to each fluid. Only one constraint per fluid can be imposed. Each constraint comes with its own subsection (starting with number ``0``) containing its own set of parameters as detailed below.

  .. warning::
      Currently, only single fluid simulations can handle this feature. Consequently, only one constraint should be defined.

  * The ``fluid id`` parameter is an integer representing the fluid on which the current constraint should be applied.

  * The ``min temperature`` parameter is a double representing the minimum temperature value for a cell to be considered a solid :math:`[\Theta]`.

  * The ``max temperature`` parameter is a double representing the maximum temperature value for a cell to be considered a solid :math:`[\Theta]`.

