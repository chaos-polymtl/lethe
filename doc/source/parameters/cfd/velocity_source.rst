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
    set permeability model               = none
    set enable Darcy multiply by density = false
    set Carman-Kozeny permeability area  = 1e-3
    set Carman-Kozeny division tolerance = 1e-4
  end

* The ``permeability model`` parameter specifies the type of penalization term to be applied to the Navier-Stokes equations. The options are ``none``, ``darcy_phase_change``, or ``carman_kozeny_phase_change``.

  - The ``darcy_phase_change`` model uses the values of the ``Darcy penalty liquid``  and ``Darcy penalty solid`` set-up within the ``phase change`` subsection of the :doc:`physical_properties`. The added source term corresponds to:

    .. math::
       \boldsymbol{F}_\mathrm{Darcy} = K \boldsymbol{u}

    where

    - :math:`K= \alpha_\mathrm{l}K_\mathrm{l} + (1-\alpha_\mathrm{l})K_\mathrm{s}` is the liquid fraction :math:`\left(\alpha_\mathrm{l} \right)` weighted Darcy penalty with :math:`K_\mathrm{l}` and :math:`K_\mathrm{s}` respectively the Darcy penalty in the liquid and the solid phases;
    - :math:`\boldsymbol{u}` is the velocity.

  - The ``carman_kozeny_phase_change`` uses the ``Carman-Kozeny permeability area``, the ``Carman-Kozeny division tolerance``, the ``kinematic viscosity`` specified in the the :doc:`physical_properties` and the liquid fraction (:math:`\alpha_\mathrm{l}`) to compute it's penalty. The added source term corresponds to:

    .. math::
       \boldsymbol{F}_\mathrm{Carman-Kozeny} = \frac{\nu}{A_\mathrm{perm}} \left[ \frac{(1-\alpha_\mathrm{l})^2}{(\alpha_\mathrm{l})^3 + \delta}\right] \boldsymbol{u}

    where

    - :math:`\nu` is the kinematic viscosity;
    - :math:`A_\mathrm{perm}` is the permeability area;
    - :math:`\delta` a tolerance to avoid division by zero in the solid, and;
    - :math:`\boldsymbol{u}` is the velocity.

    The figure below gives an idea of how the different parameters influence the resulting force term.

    +----------------------------------------------------------------------------------------------------------------------+
    |  .. figure:: images/carman_kozeny_penalty_evolution.svg                                                              |
    |     :align: center                                                                                                   |
    |     :width: 650                                                                                                      |
    |     :name: Carman-Kozeny penalty evolution                                                                           |
    |                                                                                                                      |
    |  Evolution of the Caraman-Kozeny penalty (force magnitude) in function of the liquid fraction, the permeability      |
    |  area, and the tolerance parameter for :math:`\nu = 1` and :math:`\lVert \boldsymbol{u} \rVert = 1`.                 |
    +----------------------------------------------------------------------------------------------------------------------+

    For :doc:`CLS<../../theory/multiphase/cfd/cls>` simulations, the source term takes the following form:

    .. math::
       \boldsymbol{F}_\mathrm{Carman-Kozeny} =  \left[ \frac{1}{A_\mathrm{perm}} \sum_{i=0}^1 w_i \mu_i  \left[ \frac{(1-\alpha_{i,\mathrm{l}})^2}{(\alpha_{i,\mathrm{l}})^3 + \delta}\right] \right] \boldsymbol{u}

    where

    - :math:`w_i = \begin{cases}  1-\phi \quad &\mathrm{if} \quad  i = 0\\ \phi \quad &\mathrm{if} \quad  i = 1\\ \end{cases}` is the phase indicator weight, and;
    - :math:`\mu` is the dynamic viscosity.

    .. note::
      Here, the dynamic viscosity is used, since in the CLS formulation of the Navier-Stokes equations, we solve for the pressure (:math:`p`) rather than the kinematic pressure (:math:`p^* = p / \rho`, with :math:`\rho` the density of the fluid).

  .. caution::
    The phase change Darcy and Carman-Kozeny models do not currently have a Cahn-Hilliard implementation.

* The ``enable Darcy multiply by density`` parameter enables the multiplication by the density within the Darcy force term (:math:`\boldsymbol{F}_\mathrm{Darcy}`). This is for dimensional consistency when solving the pressure (:math:`p`) rather than the kinematic pressure (:math:`p^* = p / \rho`, with :math:`\rho` the density of the fluid) in the momentum equation. This parameter should be used when coupling of the :doc:`CLS equation<../../theory/multiphase/cfd/cls>` with the :doc:`incompressible Navier-Stokes equations<../../theory/multiphysics/fluid_dynamics/navier-stokes>`.

  .. math::
    \boldsymbol{F}_\mathrm{Darcy} = \left(\rho K\right)_\mathrm{eff} \boldsymbol{u}

  where

  * :math:`\left(\rho K\right)_\mathrm{eff} = (1-\phi)\rho_0 K_\mathrm{0} + \phi \rho_1 K_\mathrm{1}` is the effective Darcy penalization coefficient. It corresponds to the weighted product of the Darcy penalty :math:`\left(K_i [=] T^{-1}\right)` and the density :math:`\left(\rho_i[=] M L^{-3}\right)` by the phase indicator (:math:`\phi`). If the ``phase_change`` model is enabled in a fluid, the Darcy penalty in this fluid is computed with the penalty values in the liquid (:math:`K_{i\mathrm{,l}}`) and solid (:math:`K_{i\mathrm{,s}}`), and the liquid fraction (:math:`\alpha_\mathrm{l}`): :math:`K_i = \alpha_\mathrm{l}K_{i\mathrm{,l}} + (1-\alpha_\mathrm{l})K_{i\mathrm{,s}}`, and;
  * :math:`\boldsymbol{u}` is the velocity vector.

* The ``Carman-Kozeny permeability area`` parameter corresponds to :math:`A_\mathrm{perm}` in :math:`\boldsymbol{F}_\mathrm{Carman-Kozeny}`. It represents the permeability area of the pseudo-porous bed (or solid phase) that is simulated.

* The ``Carman-Kozeny division tolerance`` parameter avoids division by zero in :math:`\boldsymbol{F}_\mathrm{Carman-Kozeny}` when in the solid phase (:math:`\alpha_\mathrm{l} = 0`).