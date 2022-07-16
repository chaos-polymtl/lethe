Laser heat source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a laser heat source is present in a simulation, it can be added in this section. The default parameters are:

.. code-block:: text

   subsection laser parameters
	set enable                 = false
	set concentration factor   = 2.0
	set power                  = 100.0
	set absorptivity           = 0.5
	set penetration depth      = 0.0
	set beam radius            = 0.0
	set start time             = 0.0
	set end time               = 1.0
	set beam orientation       = z-
	subsection path
	  set Function expression  =  0.0; 0.0
	end
   end

* The ``enable`` parameter is set to true if the problem has a laser heat source term and enables its calculation.

* Laser ``concentration factor`` parameter indicates the definition of the beam radius. In almost all the articles, it is assumed equal to 2.0.

* The ``power`` parameter sets the power of the laser [:math:`ML^2T^{-3}`].

* The ``absorptivity`` parameter is defined as the fraction of the amount of incident radiation that is absorbed by the surface, and it is measured using diffuse reﬂectance spectroscopy (DRS). Generally, a constant value in the range of 0.3-0.8 (for welding processes with titanium) ise used in the literature. However, recent studies show that it varies with powder particle size distribution and the angle of incidence that changes due to the dynamic melt pool surface `[1] <https://doi.org/10.1016/j.optlastec.2018.08.012>`_.

* The ``penetration depth`` parameter determines the penetration depth of the laser in the simulation domain in the direction of emission.

* The ``beam radius`` parameter defines the radius of the laser beam.

* The ``start time`` and ``end time`` parameters define the operation time window of the laser.

* The ``beam orientation`` parameter shows the orientation and direction of the laser beam. For instance, if a laser beam is emitted perpendicular to a plane in x-y coordinates, the orientation of the laser beam will be in the z-direction. Negative (-) or positive (+) defines the direction of the laser beam. For instance if the laser beam is emitted in the negative z direction, the value of ``beam orientation`` will be ``z-``.

.. note:: 
    In two-dimensional simulations, the laser beam orientation cannot be in the z-direction.


* In the ``path`` subsection, the laser scanning path is defined using a ``Function expression``.

The exponential decaying model `[2] <https://doi.org/10.1016/j.matdes.2018.01.022>`_ is used to simulate the laser heat source. In the exponential decaying model, the laser heat flux is calculated using the following equation:

    .. math:: 
        q(x,y,z) = \frac{\eta \alpha P}{\pi r^2 \mu} \exp{(-\eta \frac{r^2}{R^2})} \exp{(- \frac{|z|}{\mu})}


where :math:`\eta`, :math:`\alpha`, :math:`P`, :math:`R`, :math:`\mu`, :math:`r`, \:math:`z` denote concentration factor, absorptivity, laser power, beam radius, penetration depth, radial distance from the laser focal point, and axial distance from the laser focal point, respectively.

-----------
References
-----------
`[1] <https://doi.org/10.1016/j.optlastec.2018.08.012>`_ Zhang, Z., Huang, Y., Kasinathan, A.R., Shahabad, S.I., Ali, U., Mahmoodkhani, Y. and Toyserkani, E., 2019. 3-Dimensional heat transfer modeling for laser powder-bed fusion additive manufacturing with volumetric heat sources based on varied thermal conductivity and absorptivity. Optics & Laser Technology, 109, pp.297-312.

`[2] <https://doi.org/10.1016/j.matdes.2018.01.022>`_ Liu, S., Zhu, H., Peng, G., Yin, J. and Zeng, X., 2018. Microstructure prediction of selective laser melting AlSi10Mg using finite element analysis. Materials & Design, 142, pp.319-328.

