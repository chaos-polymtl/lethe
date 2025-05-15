====================================
Discrete Element Method (DEM)
====================================

In this guide, we summarize the theory behind DEM. For further details, we refer the reader to the article by Blais *et al.*  [#blais2019]_ and the one by Golshan *et al.* [#golshan2023]_


.. math::
    m_i\frac{d\mathbf{v_i}}{dt} &= \sum_{j\in \mathcal C_i} (\mathbf{F}_{ij}^\mathrm{n} + \mathbf{F}_{ij}^\mathrm{t}) + m_i\mathbf{g} + \mathbf{F}_i^\mathrm{ext} \\
    I_i\frac{d\mathbf{\omega_i}}{dt} &= \sum_{j\in \mathcal C_i} (\mathbf{M}_{ij}^\mathrm{t} + \mathbf{M}_{ij}^\mathrm{r}) +  \mathbf{M}_i^\mathrm{ext}

Where:

* :math:`m_i` mass of the particule i;
* :math:`\mathbf{v_i}` velocity of the particule i;
* :math:`\mathcal C_i` particles in contact list;
* :math:`\mathbf{F}_{ij}^\mathrm{n}` normal contact force due to the contact between particles i and j;
* :math:`\mathbf{F}_{ij}^\mathrm{t}` tangential contact force due to the contact between particles i and j;
* :math:`\mathbf{g}` gravitationnal acceleration;
* :math:`\mathbf{F}_i^\mathrm{ext}` external forces;
* :math:`I_i` moment of inertia of particle i;
* :math:`\mathbf{\omega_i}` angular velocity of particle i;
* :math:`\mathbf{M}_{ij}^\mathrm{t}` tangential friction torque due to the contact between particles i and j;
* :math:`\mathbf{M}_{ij}^\mathrm{r}` rolling friction torque due to the contact between particles i and j;
* :math:`\mathbf{M}_i^\mathrm{ext}` external torques;


--------------------------------
Contact Force and Torque Models
--------------------------------

The normal and tangential contact forces use linear or nonlinear viscoelastic models and are calculated as followed:

.. math::
    \mathbf{F}_{ij}^\mathrm{n} &= -(k_\mathrm{n}\delta_{\mathrm{n}})\mathbf{n}_{ij}-(\eta_\mathrm{n}\mathbf{v}_{rn}) \\
    \mathbf{F}_{ij}^\mathrm{t} &= -(k_\mathrm{t}\mathbf{\delta}_\mathrm{t})-(\eta_\mathrm{t}\mathbf{v}_{rt})

Where:

* :math:`\delta_{\mathrm{n}}` normal overlap;
* :math:`\mathbf{\delta_\mathrm{t}}` tangential displacement vector;
* :math:`\mathbf{n}_{ij}` contact normal vector;
* :math:`k_\mathrm{n}, k_\mathrm{t}` spring constants;
* :math:`\eta_\mathrm{n}, \eta_\mathrm{t}` damping model constants;
* :math:`\mathbf{v}_{rn}` relative velocity in the normal direction;
* :math:`\mathbf{v}_{tn}` relative velocity in the tangential direction;


.. figure:: images/collision_particles.png
    :width: 400
    :align: center
    :alt: particle-particle_collision

    Representation of a typical particle-particle contact. [#golshan2023]_

The contact normal vector :math:`\mathbf{n}_{ij}` is computed as:

.. math::
    \mathbf{n}_{ij}=\frac{\mathbf{x}_{j}-\mathbf{x}_{i}}{\left|\mathbf{x}_{j}-\mathbf{x}_{i}\right|}

The normal overlap (:math:`\delta_{\mathrm{n}}`) is the contact distance between the particles i and j. In the case of a collision between a particle and a wall, the wall is considered as j. The tangential displacement (:math:`\delta_\mathrm{t}`) depends on the contact history and is updated during a contact.
The normal and tangential displacements are calculated as follow:

.. math::
    \delta_{\mathrm{n}} =& \:R_i + R_j - |\mathbf{x}_{j} - \mathbf{x}_{i}| \\
    \mathbf{\delta}_{ij}^{\mathrm{t,new}} &= \mathbf{\delta}_{ij}^{\mathrm{t,old}}+\mathbf{v}_{ij,\mathrm{t}}dt

~~~~~~~~~~~~~~~~~~~~~
Relative Velocities
~~~~~~~~~~~~~~~~~~~~~
The relative velocities are calculated to update the tangential displacement and for their contribution in the force models:

.. math::
    \mathbf{v}_{ij} &= \mathbf{v}_i-\mathbf{v}_j+\left(R_i\mathbf{\omega}_i+R_j\mathbf{\omega}_j\right)\times\mathbf{n}_{ij} \\
    \mathbf{v}_{ij,\mathrm{n}} &= \left(\mathbf{v}_{ij}.\mathbf{n}_{ij}\right)\mathbf{n}_{ij} \\
    \mathbf{v}_{ij,\mathrm{t}} &= \mathbf{v}_{ij}-\mathbf{v}_{ij,\mathrm{n}}

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Spring and damping constants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The spring and damping constants for the linear and nonlinear viscoelastic models are calculated as follows:

.. list-table:: Spring and damping Models used in Lethe.
   :widths: 40 30 30
   :header-rows: 1

   * - Parameters
     - Linear model definitions
     - Nonlinear viscoelastic model definitions [#garg2012]_
   * - Normal spring constant
     - :math:`k_\mathrm{n} = \frac{16}{15}\sqrt{R_{\mathrm{e}}}Y_{e}\left(\frac{15m_{e}V^2}{16\sqrt{R_{\mathrm{e}}}Y_{e}}\right)^{0.2}`
     - :math:`k_\mathrm{n} = \frac{4}{3}Y_{e}\sqrt{R_{\mathrm{e}}\delta_{\mathrm{n}}}`
   * - Normal damping model constant
     - :math:`\eta_\mathrm{n} = -2\beta\sqrt{m_{e} k_\mathrm{n}}`
     - :math:`\eta_\mathrm{n} = -2\sqrt{\frac{5}{6}}\beta\sqrt{S_\mathrm{n}m_{e}}`
   * - Tangential spring constant
     - :math:`k_\mathrm{t} = 0.4 k_\mathrm{n}`
     - :math:`k_\mathrm{t} = 8G_{e}\sqrt{R_{\mathrm{e}}\delta_{\mathrm{n}}}`
   * - Tangential damping model constant
     - :math:`\eta_\mathrm{t} = -2\beta\sqrt{m_{e} k_\mathrm{t}}`
     - :math:`\eta_\mathrm{t} = -2\sqrt{\frac{5}{6}}\beta\sqrt{S_\mathrm{t}m_{e}}`

Where:

* :math:`R_{\mathrm{e}}` effective radius;
* :math:`Y_\mathrm{e}` effective Young's modulus;
* :math:`m_\mathrm{e}` effective mass;
* :math:`V` characteristic impact velocity, this parameters is set to 1.0;
* :math:`e` coefficient of restitution;
* :math:`G_\mathrm{e}` effective shear modulus;

These parameters are computed as follows:

.. math::
    \frac{1}{m_\mathrm{e}} &= \frac{1}{m_i}+\frac{1}{m_j} \\
    \frac{1}{R_{\mathrm{e}}} &= \frac{1}{R_i}+\frac{1}{R_j} \\
    \frac{1}{G_\mathrm{e}} &= \frac{2(2-\nu_i)(1+\nu_i)}{Y_i}+\frac{2(2-\nu_j)(1+\nu_j)}{Y_j} \\
    \frac{1}{Y_\mathrm{e}} &= \frac{\left(1-\nu_i^2\right)}{Y_i}+\frac{\left(1-\nu_j^2\right)}{Y_j} \\
    \beta &= \frac{\ln{e}}{\sqrt{\ln^2{e}+\pi^2}} \\
    S_\mathrm{n} &= 2Y_{e}\sqrt{R_{\mathrm{e}}\delta_{\mathrm{n}}} \\
    S_\mathrm{t} &= 8G_{e}\sqrt{R_{\mathrm{e}}\delta_{\mathrm{n}}}

Where:

* :math:`\nu_i, \nu_j` poisson coefficient of particle i or j;

~~~~~~~~~~~~~~~~~~~~
Coulomb's limit
~~~~~~~~~~~~~~~~~~~~

Coulomb's criterion is breached when the following condition is broken during a collision:

.. math::
    |\mathbf{F}_{ij}^{\mathrm{t}}| \geq \mu |\mathbf{F}_{ij}^\mathrm{n}|


A breach means the collision is having gross sliding and tangential force needs to be limited to the Coulomb's limit.
To do so, the tangential displacement :math:`\mathbf{\delta_\mathrm{t}}` is first limited and then the tangential force is recalculated.

When using nonlinear viscoelastic contact model, the tangential displacement is computed from tangential spring force :

.. math::
    \mathbf{\delta_\mathrm{t}} &= \frac{\mathbf{\tilde{F}_{ij}}}{-k_\mathrm{t}} \\
    \mathbf{\tilde{F}_{ij}} &= \mathbf{\hat{F}_{ij}} + \eta_\mathrm{t}\mathbf{v}_{ij,\mathrm{t}} \\
    \mathbf{\hat{F}_{ij}^\mathrm{t}} &= \mu |\mathbf{F}_{ij}^\mathrm{n}| \frac{\mathbf{F}_{ij}^\mathrm{t}}{|\mathbf{F}_{ij}^\mathrm{t}|}

Regarding the particle-wall contacts, the applied models are the same as for particle-particle contacts.

.. note::
    When using a cohesive force model, Coulomb's criterion needs to be modified. For further information on cohesive force models, see `Cohesive force models`_ .

~~~~~~~~~~~~~~~~~~~~~~~~~
Tangential torque
~~~~~~~~~~~~~~~~~~~~~~~~~

Tangential torque is the torque generated by the tangential force. It can be calculated through:

.. math::
    \mathbf{M}_{\mathrm{t},ij} = R_{i}\mathbf{n}_{ij} \times \mathbf{F}_{\mathrm{t},ij}

.. note::
    As of now, the ``lethe-particles`` solver only uses spherical particles, thus the normal force does not generate a torque on the particle during a collision.

~~~~~~~~~~~~~~~~~~~~~~~~~
Rolling friction models
~~~~~~~~~~~~~~~~~~~~~~~~~

Rolling friction may be computed through a constant torque model, a viscous torque model or an elastic-plastic spring-dashpot torque model. It is also possible to ignore the rolling resistance by using the no-resistance model. The corresponding models are described by the following equations:

.. list-table:: Rolling Friction Models used in Lethe.
   :width: 80%
   :widths: 40 40
   :header-rows: 1
   :align: center

   * - Models
     - Equations
   * - Constant resistance
     - :math:`\mathbf{M}_{\mathrm{r},ij} = -\mu_\mathrm{r}R_{\mathrm{e}}|\mathbf{F}_{\mathrm{n},ij}| \mathbf{\hat{\omega}}_{ij}`
   * - Viscous resistance
     - :math:`\mathbf{M}_{\mathrm{r},ij} = -\mu_\mathrm{r}R_{\mathrm{e}}|\mathbf{F}_{\mathrm{n},ij}||\mathbf{V}_{\omega}| \mathbf{\hat{\omega}}_{ij}`
   * - Elastic-plastic spring-dashpot resistance
     - :math:`\mathbf{M}_{\mathrm{r},ij} = \mathbf{M}_{\mathrm{t+\Delta t}}^{k} + \mathbf{M}_{\mathrm{t+\Delta t}}^{d}`
   * - No resistance
     - :math:`\mathbf{M}_{\mathrm{r},ij} = 0`

Where:

* :math:`\mu_\mathrm{r}` rolling friction coefficient;
* :math:`\hat{\omega}_{ij}` relative angular velocity;
* :math:`V_{\omega}` contact point relative velocity caused by the angular velocities;
* :math:`\mathbf{M}_{\mathrm{r,t+\Delta t}}^{k}` elastic resistance torque at the end of the current time step;
* :math:`\mathbf{M}_{\mathrm{r,t+\Delta t}}^{d}` viscous damping resistance torque at the end of the current time step;

The parameters for the constant and viscous models are computed as follows:

.. math::
    \mathbf{\hat{\omega}}_{ij} &= \frac{\omega_{i} - \omega_{j}}{|\omega_{i} - \omega_{j}|} \\
    \mathbf{V}_{\omega} &= \left( \omega_{i} \times R_{i}\mathbf{n}_{ij}-\omega_{j} \times R_{j}\mathbf{n}_{ji} \right).

For the elastic-plastic spring-dashpot model, :math:`\mathbf{M}_{\mathrm{r}}^{k}` and :math:`\mathbf{M}_{\mathrm{r}}^{d}` are computed using the following algorithm:

.. math::
    \mathbf{\omega}_{ji} &= \mathbf{\omega}_{i}- \mathbf{\omega}_{j}\\
    \mathbf{\omega}_{ji,\mathrm{plane}} &= \mathbf{\omega}_{ij}- \left( \mathbf{\omega}_{ij}\cdot\mathbf{n}_{ij} \right) \mathbf{n}_{ij}\\
    \mathbf{\Delta\theta} &= \Delta t \; \mathbf{\omega}_{ij,\mathrm{plane}}\\
    k_\mathrm{r} &= 2.25 k_\mathrm{n} \left( \mu_\mathrm{r} R_\mathrm{e} \right)^2\\
    \mathbf{\Delta M}_{\mathrm{r},t}^k &= -k_\mathrm{r}\mathbf{\Delta\theta}\\
    \mathbf{M}_{\mathrm{r},t+\Delta t}^\mathrm{k} &= \mathbf{M}_{\mathrm{r},t}^\mathrm{k}+ \mathbf{\Delta M}_{\mathrm{r},t}^\mathrm{k} \\
    M\mathrm{^{m}_{r}} &= \mu_\mathrm{r} R_\mathrm{e} |\mathbf{F}_{\mathrm{n},ij}|\\
    \mathbf{M}_{\mathrm{r},t+\Delta t}^{\mathrm{k}} &= \begin{cases}
         \mathbf{M}_{\mathrm{r},t+\Delta t}^{\mathrm{k}}, & |\mathbf{M}_{\mathrm{r},t+\Delta t}^{\mathrm{k}}| <  M\mathrm{^{m}_{r}} \\
         \frac{ \mathbf{M}_{\mathrm{r},t+\Delta t}^{\mathrm{k}}}{| \mathbf{M}_{\mathrm{r},t+\Delta t}^{\mathrm{k}}|} M\mathrm{^{m}_{r}}, & \text{else}
    \end{cases}\\
    I_\mathrm{e} &= \left( \frac{1}{I_i + m_iR_i^2} + \frac{1}{I_j + m_jR_j^2} \right)\\
    C_r^{crit} &= 2 \sqrt{I_\mathrm{e} k_\mathrm{r}} \\
    C_r &= \eta_r C_r^{crit}\\
    \mathbf{M}_{\mathrm{r},t+\Delta t}^{\mathrm{e}} &= \begin{cases}
         -C_r \mathbf{\omega}_{ij,\mathrm{plane}} , & |\mathbf{M}_{\mathrm{r},t+\Delta t}^{\mathrm{k}}| <  M\mathrm{^{m}_{r}} \\
         -f C_r \mathbf{\omega}_{ij,\mathrm{plane}}, & \text{else}
    \end{cases}

Where:

* :math:`\mathbf{\omega}_{ji}` relative angular velocity between particle j and i;
* :math:`\mathbf{\omega}_{ji,t}` relative angular velocity between particle j and i perpendicular to the normal contact vector (vector in the contact plane);
* :math:`\mathbf{\Delta\theta}` incremental relative rotation between particle j and i;
* :math:`k_\mathrm{r}` rolling stiffness;
* :math:`\mathbf{\Delta M}_{r,t}^\mathrm{k}` incremental elastic rolling resistance torque;
* :math:`M\mathrm{^{m}_{r}}` limiting spring torque which is achieved at a full angular mobilisation;
* :math:`I_\mathrm{e}` effective inertia;
* :math:`C_\mathrm{r}^{\mathrm{crit}}` rolling critical viscous damping constant;
* :math:`C_\mathrm{r}` rolling viscous damping constant;
* :math:`f` full mobilisation model parameter;

:math:`\mathbf{M}_{t}^{\mathrm{k}}` starts at :math:`\mathbf{0}` at the beginning of a contact and is set back to :math:`\mathbf{0}` when the contact ends. :math:`\mathbf{M}_{\mathrm{r},ij}` is applied on particle i. The rolling resistance torque applied on particle j can be found using Newton's Third Law.

For further details on all three rolling resistance model, we refer the reader to the article by Ai *et al.*  [#ai2011]_

-----------------------
Cohesive force models
-----------------------

Lethe supports two cohesive force models: the Johnson-Kendall-Roberts (JKR) and the Derjaguin-Muller-Toporov (DMT). Both models describe attractive forces due to van der Waals effects. Choosing the right model can be based on the Tabor parameter :math:`\mathbf{\tau}` which represents the ratio between the normal elastic deformation caused by adhesion and the distance at which adhesion forces occur. [#grierson2005]_

This parameter can be described as:


.. math::
    \mathbf{\tau} = \left( \frac{R_{\mathrm{e}} \gamma_{\mathrm{e}}^2}{Y_\mathrm{e}^2 z_{\mathrm{o}}^3}\right)^{1/3}

Where :math:`\mathbf{z_{\mathrm{o}}}` is the equilibrium separation of the surfaces and :math:`\mathbf{\gamma}_\mathrm{e}` the effective surface energy. The DMT model is applicable for low :math:`\mathbf{\tau}` values (:math:`\mathbf{\tau} < 1`) while the JKR model is more appropriate for high :math:`\mathbf{\tau}` values (:math:`\mathbf{\tau} > 1`) . In essence, the DMT model is preferred for small, hard particles (high :math:`Y`) and the JKR model for large, soft particles.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Johnson-Kendall-Roberts force model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Johnson-Kendall-Roberts (JKR) model describes attractive forces due to van der Waals effects. [#coetzee2023]_
This model modifies the Hertz formulation by defining a larger contact path radius (:math:`\mathbf{a}`) and by taking into account the effective surface energy (:math:`\mathbf{\gamma}_{e}`).
The model is defined by:

.. math::
    a^{3} = \frac{3 R_{\mathrm{e}}}{4 Y_{\mathrm{e}}} \left[|\mathbf{F_{n}^{JKR}}| + 3\pi\gamma_{\mathrm{e}}R_{\mathrm{e}}  + \sqrt{6 |\mathbf{F_{n}^{JKR}}| \pi\gamma_{\mathrm{e}}R_{\mathrm{e}} + (3\pi\gamma_{\mathrm{e}}R_{\mathrm{e}})^2 }\right]

Where :math:`\mathbf{F_{n}^\mathrm{JKR}}` corresponds to the normal spring force and attractive force combined and :math:`\mathbf{\gamma_{\mathrm{e}}}` is the effective surface energy.
Note that if the effective surface energy is equal to zero, the JKR model reverts to Hertz model.

The effective surface energy can be computed as:

.. math::
    \gamma_{\mathrm{e}} = \gamma_{i} + \gamma_{j} - 2\gamma_{i,j}

Where :math:`\gamma_{i}` and :math:`\gamma_{j}` are the surface energy of each material (particle or wall) and where :math:`\gamma_{i,j}` is the interface energy which is equal to zero when both surfaces are the same material.
The interface energy term is approximated using [#israelachvili–289]_:

.. math::
    \gamma_{i,j} \approx \left( \sqrt{\gamma_{i}} - \sqrt{\gamma_{j}}  \right)^{2}

To compute the :math:`\mathbf{F_{n}^{JKR}}`, the contact patch radius needs to be determined. The contact patch radius can be related to the normal overlap as follows:

.. math::
    \delta_{\mathrm{n}} = \frac{ a^{2} }{ R_{\mathrm{e}} } -  \sqrt{ \frac{2 \pi \gamma_{\mathrm{e}} a }{ Y_\mathrm{e} }}

This equation can be rewritten as a fourth-order polynomial function with two complex and two real roots.

.. math::
    0 = a^{4} - 2R_{\mathrm{e}}\delta_{\mathrm{n}}a^{2} - 2\pi\gamma_{\mathrm{e}}R_{\mathrm{e}}^{2}a + R_{\mathrm{e}}^{2}\delta_{\mathrm{n}}^{2}

Since we are always solving for the same real root, a straightforward procedure, described by Parteli et al. can be used [#parteli2014]_:

.. math::
    c_\mathrm{0} &= R_{\mathrm{e}}^{2}\delta_{\mathrm{n}}^{2} \\
    c_\mathrm{1} &= \frac{-2\pi\gamma_{\mathrm{e}}R_{\mathrm{e}}^{2}}{Y_{\mathrm{e}}}\\
    c_\mathrm{2} &= -2R_{\mathrm{e}}\delta_{\mathrm{n}}\\
    P &= -\frac{c_{\mathrm{2}}^{2}}{12} - c_{\mathrm{0}} \\
    Q &= - \frac{c_{\mathrm{2}}^{3}}{108} + \frac{c_{\mathrm{0}}c_{\mathrm{2}}}{3} - \frac{c_{\mathrm{1}}^{2}}{8} \\
    U &= \left[ -\frac{ Q }{ 2 } + \sqrt{  \frac{ Q^{2} } {4} + \frac{ P^{3} }{ 27 }  }  \right]^{ \frac{1}{3} } \\
    s &=
    \begin{cases}
    -5c_{\mathrm{2}}/6 + U - \frac{P}{3U} &{if}\: P \neq 0 \\
    -5c_{\mathrm{2}}/6 + Q^{\frac{1}{3}}  &{if}\: P = 0
    \end{cases}\\
    \omega &= \sqrt{c_{\mathrm{2}} + 2 s} \\
    \lambda &= \frac{c_{\mathrm{1}} }{2 \omega}\\
    a &= \frac{1}{2}\left(\omega + \sqrt{\omega^{2} - 4(c_{\mathrm{2}} + s + \lambda ) } \right)

Finally, the :math:`\mathbf{F_{\mathrm{n}}^{JKR}}` can be computed as follows:

.. math::
    \mathbf{F_{n}^{JKR}} = \left( \frac{4 Y_\mathrm{e} a^{3}}{3 R_{\mathrm{e}}} - \sqrt{8 \pi \gamma_{\mathrm{e}} Y_{\mathrm{e}} a^{3}} \right) \mathbf{n}_{ij}

The normal damping, tangential damping and tangential spring constants need to be computed using the same procedure as the nonlinear model.

A simplified version of the JKR model (SJKR-A) is implemented in Lethe. Please refer to C. J. Coetzee and O. C. Scheffler for more information on the different versions of the JKR model and their specific features [#coetzee2023]_.

A modified Coulomb's limit, based on the work of C. Thornton [#thornton1991]_, is used for the JKR model. Using the usual limit can result in permanent slip since the total normal force can be equal to zero even when there is a substantial overlap between particles.

The modified Coulomb's criterion is breached when the following condition is broken during a collision:

.. math::
    |\mathbf{F}_{ij}^{t}| \geq \mu |\mathbf{F_{n}^{JKR} + 2F_{\mathrm{po}}}|.

Where :math:`\mathbf{F_{\mathrm{po}}}` is the pull-off force, which can be computed as follows:

.. math::
    \mathbf{F_{\mathrm{po}}} = \left(1.5\pi\gamma_{\mathrm{e}}R_{\mathrm{e}}\right) \mathbf{n}_{ij}


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Derjaguin-Muller-Toporov force model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Derjaguin-Muller-Toporov (DMT) model describes attractive forces due to van der Waals effects. This model is more suitable for particles with smaller diameter, lower surface energy and higher Young's modulus. In Lethe, the DMT model is implemented following the work of Meier *et al.* [#meier2019]_. This implementation includes non-contact forces between particles. The model is described by the following equations:

.. math::
    \mathbf{F_{ad}^{DMT}} =    \begin{cases}
        F_{\mathrm{po}} = -2\pi\gamma_{\mathrm{e}}R_{\mathrm{e}}, & \delta_\mathrm{n} \leq \delta_{\mathrm{o}} \\
        \frac{-AR_{\mathrm{e}}}{6 \delta_{n}^2}, & \delta_{\mathrm{o}} < \delta_\mathrm{n} < \delta^* \\
        0, &  \delta^* \leq \delta_{\mathrm{n}}
    \end{cases}

where :math:`A` is the Hamaker constant which is used to quantify the strength of van der Waals forces. :math:`\delta_{\mathrm{o}}` represents the distance at which the van der Waals force curve equals the pull-off force :math:`F_{\mathrm{po}}` and :math:`\delta^*` represents a cut-off radius at which the van der Waals has a relative decline of :math:`C_{\mathrm{FPO}}` [#meier2019]_. They are computed using:

.. math::
    \begin{align}
        \delta_{\mathrm{o}} &= - \sqrt{\frac{ -A R_{\mathrm{e}}}{6 F_{\mathrm{po}}}}\\
        \delta^* &= \frac{\delta_{\mathrm{o}}}{ \sqrt{C_{\mathrm{FPO}}}}
    \end{align}

where :math:`C_{\mathrm{FPO}}` is a user parameter used to determined the cut-off distance at which the non-contact forces are being performed.

The Coulomb's limit threshold for the DMT model is computed in the same way as for the non-linear viscoelastic model. This means that the adhesion force term is not taken into account when computing the norm of the normal force. For further information, see `Coulomb's limit`_ .

--------------------
Integration Methods
--------------------

Two types of integration methods are implemented in Lethe-DEM:

* Explicit Euler method;
* Velocity Verlet method

Explicit Euler method is calculated as:

.. math::
    \mathbf{v}_{i}^{n+1} &= \mathbf{v}_{i}^{n} + \mathbf{a}_{i}^{n}dt \\
    \mathbf{x}_{i}^{n+1} &= \mathbf{x}_{i}^{n} + \mathbf{v}_{i}^{n}dt

And velocity Verlet method is calculated with half-step velocity as:

.. math::
    \mathbf{v}_{i}^{n+\frac{1}{2}} &= \mathbf{v}_{i}^{n} + \mathbf{a}_{i}^{n}\frac{dt}{2} \\
    \mathbf{x}_{i}^{n+1} &= \mathbf{x}_{i}^{n} + \mathbf{v}_{i}^{n+\frac{1}{2}}dt \\
    \mathbf{v}_{i}^{n+1} &= \mathbf{v}_{i}^{n+\frac{1}{2}} + \mathbf{a}_{i}^{n+1}\frac{dt}{2}


-------------
References
-------------

.. [#blais2019] \B. Blais, D. Vidal, F. Bertrand, G. S. Patience and J. Chaouki, “Experimental Methods in Chemical Engineering: Discrete Element Method—DEM,” *Can. J. Chem. Eng.*, vol. 97, pp. 1964-1973, 2019, doi: `10.1002/cjce.23501 <https://doi.org/10.1002/cjce.23501>`_\.

.. [#golshan2023] \S. Golshan, P. Munch, R. Gassmöller, M. Kronbichler, and B. Blais, “Lethe-DEM: an open-source parallel discrete element solver with load balancing,” *Comput. Part. Mech.*, vol. 10, no. 1, pp. 77–96, Feb. 2023, doi: `10.1007/s40571-022-00478-6 <https://doi.org/10.1007/s40571-022-00478-6>`_\.

.. [#garg2012] \R. Garg, J. Galvin-Carney, T. Li, and S. Pannala, “Documentation of open-source MFIX–DEM software for gas-solids flows,” Tingwen Li Dr., p. 10, Accessed: Sep. 2012, Available: https://mfix.netl.doe.gov/doc/mfix-archive/mfix_current_documentation/dem_doc_2012-1.pdf\.

.. [#ai2011] \R J. Ai, J.-F. Chen, J. M. Rotter, and J. Y. Ooi, “Assessment of rolling resistance models in discrete element simulations,” Powder Technology, vol. 206, no. 3, pp. 269–282, Jan. 2011, doi: 10.1016/j.powtec.2010.09.030.

.. [#grierson2005] \D. S. Grierson, E. E. Flater, and R. W. Carpick, “Accounting for the JKR–DMT transition in adhesion and friction measurements with atomic force microscopy,” *Journal of Adhesion Science and Technology*, vol. 19, no. 3–5, pp. 291–311, Jan. 2005, doi: `10.1163/1568561054352685 <https://doi.org/10.1163/1568561054352685>`_\.

.. [#coetzee2023] \C. J. Coetzee and O. C. Scheffler, “Review: The Calibration of DEM Parameters for the Bulk Modelling of Cohesive Materials,” *Processes*, vol. 11, no. 1, Art. no. 1, Jan. 2023, doi: `10.3390/pr11010005 <https://doi.org/10.3390/pr11010005>`_\.

.. [#israelachvili2011] \J. N. Israelachvili, “Chapter 13 - Van der Waals Forces between Particles and Surfaces,” in *Intermolecular and Surface Forces*, 3rd ed., J. N. Israelachvili, Ed., Boston: Academic Press, 2011, pp. 253–289, doi: `10.1016/B978-0-12-391927-4.10013-1 <https://doi.org/10.1016/B978-0-12-391927-4.10013-1>`_\.

.. [#parteli2014] \E. J. R. Parteli, J. Schmidt, C. Blümel, K.-E. Wirth, W. Peukert, and T. Pöschel, “Attractive particle interaction forces and packing density of fine glass powders,” *Sci Rep*, vol. 4, no. 1, Art. no. 1, Sep. 2014, doi: `10.1038/srep06227 <https://doi.org/10.1038/srep06227>`_\.

.. [#violano2018] \G. Violano, G. P. Demelio, and L. Afferrante, “On the DMT Adhesion Theory: From the First Studies to the Modern Applications in Rough Contacts.” *Procedia Structural Integrity*, vol. 12, pp. 58–70, Jan. 2018, doi: `0.1016/j.prostr.2018.11.106 <https://doi.org/10.1016/j.prostr.2018.11.106.>`_\.

.. [#thornton1991] \C. Thornton, “ Interparticle sliding in the presence of adhesion,” *Journal of Physics D: Applied Physics*, vol. 24, no. 11, pp. 1942–1946, 1991, doi: `10.1088/0022-3727/24/11/007 <https://doi.org/10.1088/0022-3727/24/11/007>`_\.

.. [#meier2019] \C. Meier, R. Weissbach, J. Weinberg, W. A. Wall, and A. John Hart, “Modeling and characterization of cohesion in fine metal powders with a focus on additive manufacturing process simulations,” *Powder Technology*, vol. 343, pp. 855–866, Feb. 2019, doi: `10.1016/j.powtec.2018.11.072 <https://doi.org/10.1016/j.powtec.2018.11.072>`_\.

