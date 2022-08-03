Velocity Source
~~~~~~~~~~~~~~~~

This subsection allows you to change the reference frame of the calculated velocity in a given system. This is especially useful when it is easier to solve a problem in a Lagrangian reference frame.

The default values are:

.. code-block:: text

    subsection velocity source
        set type         = none
        set omega_x      = 0
        set omega_y      = 0
        set omega_z      = 0
    end

* The ``type`` parameter specifies the type of reference frame that is selected. The options are ``none`` for an Eulerian reference frame or ``srf`` for a single rotating frame. By selecting ``srf`` the additional contributions of the centrifugal force and Coriolis force induced by the rotating movement are taken into account in the Navier-Stokes equations.

* The ``omega_x`` parameter is a double representing the :math:`x` component of the angular velocity of the reference frame.

* The ``omega_y`` parameter is a double representing the :math:`y` component of the angular velocity of the reference frame.

* The ``omega_z`` parameter is a double representing the :math:`z` component of the angular velocity of the reference frame.
