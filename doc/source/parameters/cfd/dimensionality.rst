==============
Dimensionality
==============
This subsection is used to automatically rescale the physical properties of a simulation for cases in which 
the fundamental dimensions (length, time, mass, temperature) that are used are not in SI units but the physical properties (e.g. thermal conductivity) are known in SI units.

.. code-block:: text

  subsection dimensionality
    set length      = 1 #meter
    set mass        = 1 #kilogramm
    set time        = 1 #second
    set temperature = 1 #Kelvin
  end

``length``, ``mass``, ``time`` and ``temperature`` are the fundamental dimensions that may be required. By default, it is assumed that SI units are used.
For example, if a mesh is made in millimeters, such that a height of :math:`0.100`m is :math:`100` in the mesh file, but the user still wishes to specify the viscosity in :math:`m^2/s`, the user can set ``length=0.001`` to automatically rescale the physical properties.
