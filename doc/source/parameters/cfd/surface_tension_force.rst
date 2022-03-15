Surface tension force
-----------------------------

Surface tension is the tendency of a liquid to maintain the minimum possible surface area.

When running a multiphase flow simulation, where at least one phase is liquid, using the volume of fluid method (VOF), enabling surface tension force ensures an accurate interface between the phases.

This subsection includes parameters related to the surface tension force.

.. code-block:: text

  subsection surface tension force
   set surface tension coefficient 	    = 0.1
   set phase fraction gradient filter 	= 0.0005
   set curvature filter			        = 0.0005
   set output auxiliary fields 		    = true
  end

.. warning::
   Do not forget to ``set surface tension force = true`` in the :doc:`multiphysics` subsection of the ``.prm``.   
  
   
* ``surface tension coefficient``: surface tension coefficient in :math:`Nm^{-1}`.


* ``phase fraction gradient filter``: this parameter (:math:`\eta_n \geq 0`) applies a `projection step <https://onlinelibrary.wiley.com/doi/full/10.1002/fld.2643>`_) to damp high frequency errors, that are magnified by differentiation, in the phase fraction gradient (:math:`\bf{\psi}`).
    .. math:: 
        \int_\Omega {\bf{v}} . {\bf{\psi}} + \eta_n \nabla {\bf{v}} . \nabla {\bf{\psi}} d\Omega = \int_\Omega {\bf{v}} . \nabla {\phi} d\Omega

where :math:`\bf{v}`, :math:`\bf{\psi}`, :math:`\eta_n \geq 0`, and :math:`\phi` denote a piecewise continuous vector-valued test function, filtered phase fraction gradient, phase fraction gradient filter value, and phase fraction, respectively.


.. tip::

  This parameter must be a small value larger than 0. We recommend the following procedure to choose a proper value for ``phase fraction gradient filter`` value: 1, Enable ``output auxiliary fields`` 2, Choose a small value for :math:`\eta_n = h/10`, where :math:`h` is the smallest mesh size 3, Run the simulation and check whether the filtered phase fraction gradient field is smooth and without oscillation 4, If the filtered phase fraction gradient field shows oscillations, increase the value :math:`\eta_n` to a larger value (:math:`\eta_n = h/5`, for example), and repeat this process until reaching a smooth filtered phase fraction gradient field without oscillations.

* ``curvature filter``: this parameter (:math:`\eta_k`) applies a `projection step <https://onlinelibrary.wiley.com/doi/full/10.1002/fld.2643>`_) to damp high frequency errors, that are magnified by differentiation, in the curvature (:math:`k`).
    .. math:: 
        \int_\Omega v k + \eta_k \nabla v . \nabla k d\Omega = \int_\Omega \nabla v . \frac{\bf{\psi}}{|\bf{\psi}|} d\Omega

where :math:`v`, :math:`k`, and :math:`\eta_k \geq 0` denote a test function, filtered curvature, and curvature filter value, respectively.

.. tip::

  The same procedure as choosing a value for ``phase fraction gradient filter`` should be done for ``curvature filter`` value.

* ``output auxiliary fields``: the solver writes filtered phase fraction gradient and filtered curvature if this parameter is enabled.

.. seealso::

  The :doc:`../../examples/multiphysics/rising-bubble-VOF/rising-bubble-VOF` example using VOF represents the application of surface tension force.
