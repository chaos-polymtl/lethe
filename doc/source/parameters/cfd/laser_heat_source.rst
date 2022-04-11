Laser heat source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a laser heat source is present in a simulation, it can be added in this section. The default parameters are:

.. code-block:: text

   subsection laser parameters
	set enable                 = false
	set concentration factor   = 2.0
	set power                  = 100.0
	set absorptivity           = 0.5
	set penentration depth     = 0.0
	set beam radius            = 0.0
	set start time             = 0.0
	set end time               = 1.0
	set beam orientation       = z
	subsection path
	  set Function expression  =  0.0; 0.0
	end
   end

* The ``enable`` parameter is set to true if the problem has a laser heat source term and enables its calculation.

* Laser ``concentration factor`` parameter indicates the definition of the beam radius. In almost all the articles, it is assumed equal to 2.0.

* The ``power`` parameter sets the power of the laser in W.

* The ``absorptivity`` parameter is defined as the fraction of the amount of incident radiation that is absorbed by the surface, and it is measured using diffuse reﬂectance spectroscopy (DRS). Generally, a constant value in the range of 0.3-0.8 (for welding processes with titanium) are used in the literature. However, recent studies show that it varies with powder particle size distribution, and the angle of incidence that changes due to the dynamic melt pool surface [1].

* The ``penentration depth`` parameter determines the penetration depth of the laser in the simulation domain in the direction of emission.

* The ``beam radius`` parameter defines the radius of the laser beam.

* The ``start time`` and ``end time`` parameters define the operation time window of the laser.

* The ``beam orientation`` parameter shows the orientation of the laser beam. For instance, if a laser beam is emitted perpendicular on a plane in x-y coordinates, the orientation of the laser beam will be in the z direction.

.. note:: 
    In two-dimensional simulations, the laser beam orientation is always in the perpendicular direction to the simulation domain.


* In the ``path`` subsection, the laser scanning path is defined using a ``Function expression``.


-----------
References
-----------
[1] Zhang, Z., Huang, Y., Kasinathan, A.R., Shahabad, S.I., Ali, U., Mahmoodkhani, Y. and Toyserkani, E., 2019. 3-Dimensional heat transfer modeling for laser powder-bed fusion additive manufacturing with volumetric heat sources based on varied thermal conductivity and absorptivity. Optics & Laser Technology, 109, pp.297-312.

