Dynamic flow control
~~~~~~~~~~~~~~~~~~~~

This subsection's purpose is to enable a dynamic flow control. It enables setting a volumetric flow on a specific boundary toward to normal of the wall. 
It calculates a Beta coefficient at each time step and may be added to a source term. If the chosen wall boundary on which a volumetric flow is imposed is at inlet, 
the volumetric flow rate targeted must be negative since outward normal vector is used. Thereby, if the chosen wall boundary is at outlet, 
volumetric flow rate must be positive.

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

  This value should be negative if imposed on an inlet boundary, and positive if imposed on an outlet boundary.

.. note::

  When there is a periodic boundary, the wall boundary should be at inlet only and if the chosen wall boundary is at inlet, flow rate must be negative and vice versa.

* The ``boundary id`` parameter is the wall boundary where the flow rate intended is controlled (inlet or outlet).

* The ``flow direction`` parameter sets the absolute direction. It should be ``0`` if the flow is in the :math:`x` direction, ``1`` in the :math:`y` direction and ``2`` in the :math:`z` direction.

* The ``initial beta`` parameter sets an initial :math:`\beta` coefficient that may speed up the convergence of the flow rate targeted. Some tests should be done to find a good initial :math:`\beta` coefficient.

