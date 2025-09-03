===============
Velocity Source
===============

This subsection allows you to define velocity-dependent source terms. Two families of velocity source terms are supported. They enable the simulation of  problems in a Lagrangian reference frame for which the Coriolis and centrifugal forces must be added or of problems including porous media.


Rotating Frame
~~~~~~~~~~~~~~

Rotating frame simulations use the following parameters:

.. code-block:: text

  subsection velocity source
    set rotating frame type    = none
    set omega_x                = 0
    set omega_y                = 0
    set omega_z                = 0
  end

* The ``rotating frame type`` parameter specifies the type of reference frame that is selected. The options are ``none`` for an Eulerian reference frame or ``srf`` for a single rotating frame. By selecting ``srf`` the additional contributions of the centrifugal force and Coriolis force induced by the rotating movement are taken into account in the Navier-Stokes equations.

* The ``omega_x`` parameter is a double representing the :math:`x` component of the angular velocity of the reference frame.

* The ``omega_y`` parameter is a double representing the :math:`y` component of the angular velocity of the reference frame.

* The ``omega_z`` parameter is a double representing the :math:`z` component of the angular velocity of the reference frame.


Darcy Penalization
~~~~~~~~~~~~~~~~~~

A Darcy-like source term can be added to the simulation using the following parameters:

.. code-block:: text

  subsection velocity source
    set Darcy type = none
  end

* The ``Darcy type`` parameter specifies the type of Darcy penalization term to be applied to the Navier-Stokes equations. The options are ``none`` or ``phase_change``. The ``phase_change`` model uses the values of the ``Darcy penalty liquid``  and ``Darcy penalty solid`` set-up within the ``phase change`` subsection of the :doc:`physical_properties`.

.. caution::
  The phase change Darcy model does not currently have a Cahn-Hilliard implementation.