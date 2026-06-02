..
  SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

===============
Velocity Source
===============

This subsection allows you to define velocity-dependent source terms. Two families of velocity source terms are supported. They enable the simulation of  problems in a Lagrangian reference frame for which the Coriolis and centrifugal forces must be added or of problems including porous media.


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


Darcy Penalization
~~~~~~~~~~~~~~~~~~

A Darcy-like source term can be added to the simulation using the following parameters:

.. code-block:: text

  subsection velocity source
    set Darcy type                       = none
    set enable Darcy multiply by density = false
  end

* The ``Darcy type`` parameter specifies the type of Darcy penalization term to be applied to the Navier-Stokes equations. The options are ``none`` or ``phase_change``. The ``phase_change`` model uses the values of the ``Darcy penalty liquid``  and ``Darcy penalty solid`` set-up within the ``phase change`` subsection of the :doc:`physical_properties`.

  .. caution::
    The phase change Darcy model does not currently have a Cahn-Hilliard implementation.

* The ``enable Darcy multiply by density`` parameter enables the multiplication by the density within the Darcy force term (:math:`\boldsymbol{F}_\mathrm{Darcy}`). This is for dimensional consistency when solving the pressure (:math:`p`) rather than the kinematic pressure (:math:`p^* = p / \rho`, with :math:`\rho` the density of the fluid) in the momentum equation. This parameter should be used when coupling of the :doc:`CLS equation<../../theory/multiphase/cfd/cls>` with the :doc:`incompressible Navier-Stokes equations<../../theory/multiphysics/fluid_dynamics/navier-stokes>`.

  .. math::
    \boldsymbol{F}_\mathrm{Darcy} = \left(\rho K\right)_\mathrm{eff} \boldsymbol{u}

  where

  * :math:`\left(\rho K\right)_\mathrm{eff} = (1-\phi)\rho_0 K_\mathrm{0} + \phi \rho_1 K_\mathrm{1}` is the effective Darcy penalization coefficient. It corresponds to the weighted product of the Darcy penalty :math:`\left(K_i [=] T^{-1}\right)` and the density :math:`\left(\rho_i[=] M L^{-3}\right)` by the phase indicator (:math:`\phi`). If the ``phase_change`` model is enabled in a fluid, the Darcy penalty in this fluid is computed with the penalty values in the liquid (:math:`K_{i\mathrm{,l}}`) and solid (:math:`K_{i\mathrm{,s}}`), and the liquid fraction (:math:`\alpha_\mathrm{l}`): :math:`K_i = \alpha_\mathrm{l}K_{i\mathrm{,l}} + (1-\alpha_\mathrm{l})K_{i\mathrm{,s}}`, and;
  * :math:`\boldsymbol{u}` is the velocity vector.