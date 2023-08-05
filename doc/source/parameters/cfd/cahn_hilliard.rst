===============
Cahn-Hilliard
===============

In this subsection, the parameters for multiphase flow simulation using the Cahn-Hilliard equations are specified. 

In this method, the two fluids considered are given an index of :math:`0` and :math:`1` respectively. The amount of fluid at any given quadrature point is represented by a phase order parameter :math:`\phi` between :math:`-1` and :math:`1`. The interface between the two fluids is naturally diffuse and the motion of the interface are driven by differences in the local chemical potential :math:`\eta`.

.. note::
    A Canh-Hilliard theory documentation explaining the origins of the diffuse interface will be added in an upcoming update.

The equations solved are as follows:

.. math::
        \partial_t\phi + (\mathbf{a} \cdot \nabla) \phi - \nabla \cdot (M(\phi)\nabla \eta) = 0 \\
        
         \eta - 4W(\phi^3 - \phi) + \epsilon^2 \Delta \phi = 0

where:

* :math:`\mathbf{a}` corresponds to the velocity field; this vector field is used when the problem is driven by convection.
* :math:`M(\phi)` is the mobility function; two cases are considered: a constant mobility model, i.e, :math:`M = D` with :math:`D` the mobility constant, or a quartic mobility model with :math:`M(\phi) = D(1-\phi^2)^2`.
* :math:`W` is the potential well height.
* :math:`\epsilon` is a parameter linked to the thickness of the interface.

.. note::

  At the moment, a maximum of two fluids is supported. By convention, air is usually the ``fluid 0`` and the other fluid of interest is the ``fluid 1``.    See :doc:`initial_conditions` for the definition of the Cahn-Hilliard initial conditions and :ref:`Physical properties - Two phase simulations<two phase simulations>` for the definition of the physical properties of both fluids.  Do not forget to ``set cahn hilliard = true`` in the :doc:`multiphysics` subsection of the ``.prm`` file.


The default values of the Cahn-Hilliard parameters are given in the text box below:

.. code-block:: text

  subsection cahn hilliard
  
    set well height = 1

    subsection epsilon
      set method = automatic
      set value  = 1
    end
  end
  
* ``well height``: defines the :math:`W` parameter in the equations above. Its value should be increased if the solution is excessively diffused.

* ``epsilon``: defines the :math:`\epsilon` parameter. It can either be user-defined or determined automatically for each cell. For the latter, epsilon is equal to two times the characteristic length of the cell. The choices are ``automatic`` (default) or ``manual``.

.. attention::
     The ``mobility model`` and ``mobility constant`` must be set in the :doc:`physical_properties` section.
