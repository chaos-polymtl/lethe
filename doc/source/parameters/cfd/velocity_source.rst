..
  SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

===============
Velocity Source
===============

This subsection allows you to define velocity-dependent source terms. Two families of velocity source terms are supported: rotating frames and permeability models. They respectively enable you to simulate problems in a Lagrangian reference frame where Coriolis and centrifugal forces must be taken into account, and problems that include solid-liquid phase-changes.


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


Permeability Models
~~~~~~~~~~~~~~~~~~~

A permeability source term can be added to the simulation using the following parameters:

.. code-block:: text

  subsection velocity source
    set permeability model                    = none
    set enable Darcy multiply by density      = false
    set Carman-Kozeny fluid with phase change = fluid 0
    set Carman-Kozeny permeability area       = 1e-3
    set Carman-Kozeny division tolerance      = 1e-3
  end

* The ``permeability model`` parameter specifies the type of penalization term to be applied to the Navier-Stokes equations. The options are ``none``, ``darcy phase change``, or ``carman-kozeny phase change``.

  .. caution::
    The phase change Darcy and Carman-Kozeny models do not currently have a Cahn-Hilliard implementation.

Darcy Phase Change
++++++++++++++++++

The ``darcy phase change`` model uses the values of the ``Darcy penalty liquid``  and ``Darcy penalty solid`` set-up within the ``phase change`` subsection of the :doc:`physical_properties`. The added source term corresponds to:

  .. math::
     \boldsymbol{F}_\mathrm{Darcy} = K \boldsymbol{u}

  where

  - :math:`K= \alpha_\mathrm{l}K_\mathrm{l} + (1-\alpha_\mathrm{l})K_\mathrm{s} \, [\mathsf{T^{-1}}]` is the liquid fraction :math:`\left(\alpha_\mathrm{l} \right)` weighted Darcy penalty with :math:`K_\mathrm{l}` and :math:`K_\mathrm{s}` respectively the Darcy penalty in the liquid and the solid phases;
  - :math:`\boldsymbol{u} \, [\mathsf{LT^{-1}}]` is the velocity.



* The ``enable Darcy multiply by density`` parameter enables the multiplication by the density within the Darcy force term (:math:`\boldsymbol{F}_\mathrm{Darcy}`). This is for dimensional consistency when solving the pressure (:math:`p`) rather than the kinematic pressure (:math:`p^* = p / \rho`, with :math:`\rho` the density of the fluid) in the momentum equation. This parameter should be used when coupling of the :doc:`CLS equation<../../theory/multiphase/cfd/cls>` with the :doc:`incompressible Navier-Stokes equations<../../theory/multiphysics/fluid_dynamics/navier-stokes>`.

    .. math::
      \boldsymbol{F}_\mathrm{Darcy} = \left(\rho K\right)_\mathrm{eff} \boldsymbol{u}

    where

    * :math:`\left(\rho K\right)_\mathrm{eff} = (1-\phi)\rho_0 K_\mathrm{0} + \phi \rho_1 K_\mathrm{1}` is the effective Darcy penalization coefficient. It corresponds to the weighted product of the Darcy penalty :math:`\left(K_i \, [\mathsf{T^{-1}}]\right)` and the density :math:`\left(\rho_i \, [\mathsf{M L^{-3}}]\right)` by the phase indicator (:math:`\phi`). If the ``phase_change`` model is enabled in a fluid, the Darcy penalty in this fluid is computed with the penalty values in the liquid (:math:`K_{i\mathrm{,l}}`) and solid (:math:`K_{i\mathrm{,s}}`), and the liquid fraction (:math:`\alpha_\mathrm{l}`): :math:`K_i = \alpha_\mathrm{l}K_{i\mathrm{,l}} + (1-\alpha_\mathrm{l})K_{i\mathrm{,s}}`, and;
    * :math:`\boldsymbol{u}` is the velocity vector.


Carman-Kozeny Phase Change
++++++++++++++++++++++++++

The ``carman-kozeny phase change`` uses the ``Carman-Kozeny permeability area``, the ``Carman-Kozeny division tolerance``, the ``kinematic viscosity`` specified in the the :doc:`physical_properties` and the liquid fraction (:math:`\alpha_\mathrm{l}`) to compute it's penalty. The added source term corresponds to [#zenz2024]_:

  .. math::
     \boldsymbol{F}_\mathrm{Carman-Kozeny} = \frac{\nu}{A_\mathrm{perm}} \left[ \frac{(1-\alpha_\mathrm{l})^2}{(\alpha_\mathrm{l})^3 + \delta}\right] \boldsymbol{u}

  where

  - :math:`\nu \, [\mathsf{L^2T^{-1}}]` is the kinematic viscosity;
  - :math:`A_\mathrm{perm} \, [\mathsf{L^2}]` is the permeability area;
  - :math:`\delta` a tolerance to avoid division by zero in the solid, and;
  - :math:`\boldsymbol{u} \, [\mathsf{LT^{-1}}]` is the velocity.

  .. note::
    Here, the kinematic viscosity is used instead of the dynamic viscosity alike [#zenz2024]_, since in the single fluid formulation of the Navier-Stokes equations, we solve for the kinematic pressure (:math:`p^* = p / \rho`, with :math:`\rho` the density of the fluid) rather than the pressure (:math:`p`).

The figure below gives an idea of how the different parameters influence the resulting force term. It can be seen that for different values of :math:`\delta`, the evolution of the Caraman-Kozeny force magnitude (:math:`\lVert \boldsymbol{F}_\mathrm{Carman-Kozeny} \rVert`) varies mostly within the range :math:`\alpha_\mathrm{l} \in [0,0.5]`.

+----------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/carman_kozeny_penalty_evolution.svg                                                              |
|     :align: center                                                                                                   |
|     :width: 650                                                                                                      |
|     :name: Carman-Kozeny penalty evolution                                                                           |
|                                                                                                                      |
|  Evolution of the Caraman-Kozeny force magnitude (:math:`\lVert \boldsymbol{F}_\mathrm{Carman-Kozeny} \rVert`) in    |
|  function of the liquid fraction, the permeability  area, and the tolerance parameter for :math:`\nu = 1`,           |
|  :math:`\rho =1`, and :math:`\lVert \boldsymbol{u} \rVert = 1`. The red curve highlights the iso-contour             |
|  :math:`\lVert \boldsymbol{F}_\mathrm{Carman-Kozeny} \rVert = 1`. The orange horizontal lines delimit the range of   |
|  typical :math:`A_\mathrm{perm}` that could be used, in this case :math:`[10^{-3},10^{-6}]`.                         |
+----------------------------------------------------------------------------------------------------------------------+


For :doc:`CLS<../../theory/multiphase/cfd/cls>` simulations, the source term takes the following form:

  .. math::
     \boldsymbol{F}_\mathrm{Carman-Kozeny} =  \left[ \sum_{i=0}^1 \frac{w_i \mu_i}{A_{i,\mathrm{perm}}}  \left[ \frac{(1-\alpha_{i,\mathrm{l}})^2}{(\alpha_{i,\mathrm{l}})^3 + \delta}\right] \right] \boldsymbol{u}

  where

  - :math:`w_i = \begin{cases}  1-\phi \quad &\mathrm{if} \quad  i = 0\\ \phi \quad &\mathrm{if} \quad  i = 1\\ \end{cases}` is the phase indicator weight, and;
  - :math:`\mu \, [\mathsf{ML^{-1}T^{-1}}]` is the dynamic viscosity.

* The ``Carman-Kozeny fluid with phase change`` specifies on which fluid(s) the permeability model should be applied on. The options are ``fluid 0``, ``fluid 1``, or ``both``.

  .. note::
    This only affects the ``carman-kozeny phase change`` permeability model. For the ``darcy phase change`` model, the penalization is computed according to the ``Darcy penalty liquid``  and ``Darcy penalty solid`` as described above.

  .. attention::
    When ``Carman-Kozeny fluid with phase change`` is set to ``fluid 1`` or ``both``, ensure that the ``cls`` physics is enabled in the :doc:`multiphysics` subsection.

* The ``Carman-Kozeny permeability area`` parameter corresponds to :math:`A_\mathrm{perm}` in :math:`\boldsymbol{F}_\mathrm{Carman-Kozeny}`. It represents the permeability area of the pseudo-porous bed (or solid phase) that is simulated. Typically the value of :math:`A_\mathrm{perm}` is chosen in function of :math:`\mu`, such that :math:`\frac{\mu}{A_\mathrm{perm}} \in [10^{3}, 10^{6}]`.

  .. caution::
    When ``Carman-Kozeny fluid with phase change = both``, two values of ``Carman-Kozeny permeability area`` separated by a comma must be specified. The first value corresponds to ``fluid 0``, and the second to ``fluid 1``.

* The ``Carman-Kozeny division tolerance`` parameter avoids division by zero in :math:`\boldsymbol{F}_\mathrm{Carman-Kozeny}` when in the solid phase (:math:`\alpha_\mathrm{l} = 0`). Typically, :math:`\delta \in [10^{-3},10^{-6}]`.

  .. caution::
    When ``Carman-Kozeny fluid with phase change = both``, two values of ``Carman-Kozeny division tolerance`` separated by a comma must be specified. The first value corresponds to ``fluid 0``, and the second to ``fluid 1``.


References
~~~~~~~~~~

.. [#zenz2024] \C. Zenz *et al.*, “A compressible multiphase Mass-of-Fluid model for the simulation of laser-based manufacturing processes,” *Computers & Fluids*, vol. 268, p. 106109, Jan. 2024, doi: `10.1016/j.compfluid.2023.106109 <https://doi.org/10.1016/j.compfluid.2023.106109>`_\.
