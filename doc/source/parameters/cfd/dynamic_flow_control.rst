Dynamic flow control
~~~~~~~~~~~~~~~~~~~~

The purpose of this subsection is to enable dynamic flow control. It is important when we want to target a volumetric flow on a specific boundary. To control the volumetric flow, the code calculates a :math:`\beta`  coefficient at each time step that is used to keep the volumetric flow value inside an acceptable threshold. If the chosen wall boundary on which a volumetric flow is imposed is an inlet, the volumetric flow rate targeted must be negative since the outward normal vector is used for the calculation. Thereby, if the chosen wall boundary is an outlet, volumetric flow rate must be positive.

The default parameters are:

.. code-block:: text

  subsection flow control
    set enable = false
    set volumetric flow rate = 0
    set boundary id = 0
    set flow direction = 0
    set initial beta = 0
  end

* The ``enable`` parameter is set to ``true`` to enable an imposed volumetric flow rate.

* The ``volumetric flow rate`` parameter specifies the flow rate (:math:`m^{dim}/s`). 

.. warning::

  This value should be negative if imposed on an the inlet boundary, and positive if imposed on an outlet boundary.

.. note::

  If you are simulating a problem with periodic boundary conditions at inlet and outlet walls, it is important to control the volumetric flow rate specifying the boundary id that corresponds to the inlet only (where the flow would have a negative value).

* The ``boundary id`` parameter is the wall boundary where the flow rate intended is controlled (inlet or outlet).

* The ``flow direction`` parameter sets the absolute direction. It should be ``0`` if the flow is in the :math:`x` direction, ``1`` in the :math:`y` direction and ``2`` in the :math:`z` direction.

* The ``initial beta`` parameter sets an initial :math:`\beta` coefficient that may speed up the convergence of the flow rate targeted. Some tests should be done to find a good initial :math:`\beta` coefficient.

.. tip:: 

  A good method to find a reasonable initial beta is to test two or three different initial beta parameters, write down the given flow rate at the first time step in the simulation and do a regression. The correlation is linear and giving a proper value will greatly speed up the convergence.

