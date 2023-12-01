=================
Laser Heat Source
=================

If a laser heat source is present in a simulation, it can be added in this section. The default parameters are:

.. code-block:: text

  subsection laser parameters
    set enable               = false
    set type                 = exponential_decay
    set concentration factor = 2.0
    set power                = 0.0
    set absorptivity         = 0.5
    set penetration depth    = 1.0
    set beam radius          = 0.0
    set start time           = 0.0
    set end time             = 1.0
    set beam orientation     = z-

    subsection path
      set Function expression = 0.0; 0.0
    end

    subsection free surface radiation
      set enable                    = false
      set emissivity                = 0.6
      set Tinf                      = 0.0
      set Stefan-Boltzmann constant = 5.6703e-8
    end
  end


* The ``enable`` parameter is set to ``true`` if the problem has a laser heat source term and enables its calculation.

* The ``type`` parameter is set to ``exponential_decay`` (default) if we assume that the laser behaves as a volumetric source.  If the laser is assumed to be a surface flux, the ``type`` can be set at ``heat_flux_vof_interface`` and used in conjunction with the :doc:`VOF auxiliary physic <./volume_of_fluid>`. The different models are detailed :ref:`below <LaserTypes>`.

* Laser ``concentration factor`` parameter indicates the definition of the beam radius. In almost all the articles, it is assumed equal to :math:`2.0`.

* The ``power`` parameter sets the power of the laser :math:`[ML^2T^{-3}]`.

* The ``absorptivity`` parameter is defined as the fraction of incident radiation that is absorbed by the surface, and it is measured using diffuse reﬂectance spectroscopy (DRS). Generally, a constant value in the range of :math:`0.3`-:math:`0.8` (for welding processes with titanium) is used in the literature. However, recent studies show that it varies with powder particle size distribution and the angle of incidence that changes due to the dynamic melt pool surface `[1] <https://doi.org/10.1016/j.optlastec.2018.08.012>`_.

* The ``penetration depth`` parameter determines the penetration depth of the laser in the simulation domain in the direction of emission.

  .. attention::
    The ``penetration depth`` value should be greater than :math:`0` and it is only taken into account if the laser ``type`` is set to ``exponential_decay``.

* The ``beam radius`` parameter defines the radius of the laser beam.

* The ``start time`` and ``end time`` parameters define the operation time window of the laser.

* The ``beam orientation`` parameter shows the orientation and direction of the laser beam. For instance, if a laser beam is emitted perpendicular to a plane in :math:`x`-:math:`y` coordinates, the orientation of the laser beam will be in the z-direction. Negative (-) or positive (+) defines the direction of the laser beam. For instance if the laser beam is emitted in the negative :math:`z` direction, the value of ``beam orientation`` will be ``z-``.

  .. attention::
      In two-dimensional simulations, the laser beam orientation cannot be in the z-direction.


* In the ``path`` subsection, the laser scanning path is defined using a ``Function expression``.

* ``subsection free surface radiation``: In additive manufacturing simulations, radiation at the interface between the air and the metal is a significant cooling mechanism. When this interface (i.e., free surface) is resolved by the :doc:`volume_of_fluid` solver, the ``free surface radiation`` subsection defines the parameters to impose this radiation cooling following the Stefan-Boltzmann law of radiation:

  .. math::
      q_\text{rad} = \epsilon \sigma (T^4 - T_\text{inf}^4)

  * ``enable``: controls if the radiation cooling is enabled. The radiation sink is modulated by the filtered phase fraction gradient norm, :math:`|\nabla \phi'|`, in such way that the flux is applied at the interface between the fluids.

    .. warning::
        To apply this radiation cooling, the ``VOF`` parameter must be set to ``true`` in the :doc:`multiphysics` subsection.

  * ``emissivity``, ``Tinf``, and ``Stefan-Boltzmann constant`` are respectively the emissivity :math:`\epsilon` of the surface, the environment temperature :math:`T_\text{inf}`, and the Stefan-Boltzmann constant :math:`\sigma`.

.. _LaserTypes:

Laser types
^^^^^^^^^^^^^

* When the ``type`` parameter is set to ``exponential_decay``, the exponential model from Zhang *et al.* `[2] <https://doi.org/10.1016/j.matdes.2018.01.022>`_ is used to simulate the laser heat source:

  .. math::
      q(x,y,z) = \frac{\eta \alpha P}{\pi r^2 \mu} \exp{\left(-\eta \frac{r^2}{R^2}\right)} \exp{\left(- \frac{|z|}{\mu}\right)}

  where :math:`\eta`, :math:`\alpha`, :math:`P`, :math:`R`, :math:`\mu`, :math:`r`, and :math:`z` denote the concentration factor, absorptivity, laser power, beam radius, penetration depth, radial distance from the laser focal point, and axial distance from the laser focal point, respectively.

  When the ``exponential_decay`` is used in conjunction with the :doc:`VOF auxiliary physic <./volume_of_fluid>` the equation takes the following form:

  .. math::
      q(x,y,z) = \frac{\phi' \eta \alpha P}{\pi r^2 \mu} \exp{\left(-\eta \frac{r^2}{R^2}\right)} \exp{\left(- \frac{|z|}{\mu}\right)}

  where :math:`\phi'` is the filtered phase fraction.

  .. attention::
    In this case, the heat affects the fluid initialized as ``fluid 1``.

* When ``type`` is set to ``heat_flux_vof_interface``, it **must be used in conjunction with the** :doc:`VOF auxiliary physic <./volume_of_fluid>`. This model is used to apply the heat flux, given by the expression below, only at the interface.

  .. math::
      q(x,y,z) = \frac{|\nabla \phi'| \eta \alpha P}{\pi r^2} \exp{\left(-\eta \frac{r^2}{R^2}\right)}

  where :math:`r` is the radial distance from the laser's axis and :math:`|\nabla \phi'|` is the :math:`L^2` norm of the filtered phase fraction gradient.


-----------
References
-----------
`[1] <https://doi.org/10.1016/j.optlastec.2018.08.012>`_ Z. Zhang *et al.*, “3-Dimensional heat transfer modeling for laser powder-bed fusion additive manufacturing with volumetric heat sources based on varied thermal conductivity and absorptivity,” *Opt. Laser Technol.*, vol. 109, pp. 297–312, Jan. 2019, doi: 10.1016/j.optlastec.2018.08.012.

`[2] <https://doi.org/10.1016/j.matdes.2018.01.022>`_ 	S. Liu, H. Zhu, G. Peng, J. Yin, and X. Zeng, “Microstructure prediction of selective laser melting AlSi10Mg using finite element analysis,” *Mater. Des.*, vol. 142, pp. 319–328, Mar. 2018, doi: 10.1016/j.matdes.2018.01.022.
