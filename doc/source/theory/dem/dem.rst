Discret Element Method (DEM)
############################
*under construction*

In this guide, we summarize the theory behind DEM. For further details, we refer the reader to the article by Golshan, Munch, Gassmöller, Kronbichler & Blais (2022) `[1] <https://doi.org/10.1007/s40571-022-00478-6>`_

.. math::
    m_i\frac{d\mathbf{v_i}}{dt} &= \sum_{j\in \mathcal C_i} (\mathbf{F}_{ij}^n + \mathbf{F}_{ij}^t) + m_i\mathbf{g} + \mathbf{F}_i^\text{ext} \\
    I_i\frac{d\mathbf{\omega_i}}{dt} &= \sum_{j\in \mathcal C_i} (\mathbf{M}_{ij}^t + \mathbf{M}_{ij}^r) +  \mathbf{M}_i^\text{ext}

Where:

* :math:`m_i` mass of the particule i;
* :math:`\mathbf{v_i}` velocity of the particule i;
* :math:`\mathcal C_i` particles in contact list;
* :math:`\mathbf{F}_{ij}^n` normal contact force due to the contact between particles i and j;
* :math:`\mathbf{F}_{ij}^t` tangential contact force due to the contact between particles i and j;
* :math:`\mathbf{g}` gravitationnal acceleration;
* :math:`\mathbf{F}_i^\text{ext}` external forces;
* :math:`I_i` moment of inertia of particle i;
* :math:`\mathbf{\omega_i}` angular velocity of particle i;
* :math:`\mathbf{M}_{ij}^t` tangential friction torque due to the contact between particles i and j;
* :math:`\mathbf{M}_{ij}^r` rolling friction torque due to the contact between particles i and j;
* :math:`\mathbf{M}_i^\text{ext}` external torques;



Contact force and torque models
--------------------------------

The normal and tangential contact forces use linear or nonlinear viscoelastic models and are calculated as followed:

.. math::
    \mathbf{F}_{ij}^n &= -(k_n\delta_n)\mathbf{n}_{ij}-(\eta_n\mathbf{v}_{rn}) \\
    \mathbf{F}_{ij}^t &= -(k_t\mathbf{\delta}_t)-(\eta_t\mathbf{v}_{rt})

Where:

* :math:`\delta_n` normal overlap;
* :math:`\mathbf{\delta_t}` tangential overlap vector;
* :math:`\mathbf{n}_{ij}` contact normal vector;
* :math:`k_n, k_t` spring constants;
* :math:`\eta_n, \eta_t` damping model constants;
* :math:`\mathbf{v}_{rn}` relative velocity in the normal direction;
* :math:`\mathbf{v}_{tn}` relative velocity in the tangential direction;


.. figure:: images/collision_particles.png
    :width: 400
    :align: center
    :alt: particle-particle_collision

    Representation of a typical particle-particle contact. [1]

The contact normal vector :math:`\mathbf{n}_{ij}` is computed as:

.. math::
    \mathbf{n}_{ij}=\frac{\mathbf{x}_{j}-\mathbf{x}_{i}}{\left|\mathbf{x}_{j}-\mathbf{x}_{i}\right|}

The normal overlap (:math:`\delta_n`) is the contact distance between the particles i and j. The tangential overlap (:math:`\delta_t`) depends on the contact history and is updated during a contact.
The normal and tangential overlaps are calculated as follow:

.. math::
    \delta_n =& \:R_i + R_j - |\mathbf{x}_{j} - \mathbf{x}_{i}| \\
    \mathbf{\delta}_{ij}^{t,\text{new}} &= \mathbf{\delta}_{ij}^{t,\text{old}}+\mathbf{v}_{rt}dt


The relative velocities calculated to allow updating the tangential overlap are described by:

.. math::
    \mathbf{v}_{rn} &= \left(\mathbf{v}_{ij}.\mathbf{n}_{ij}\right)\mathbf{n}_{ij} \\
    \mathbf{v}_{rt} &= \mathbf{v}_{ij}-\mathbf{v}_{rn} \\
    \mathbf{v}_{ij} &= \mathbf{v}_i-\mathbf{v}_j+\left(R_i\mathbf{\omega}_i+R_j\mathbf{\omega}_j\right)\times\mathbf{n}_{ij}



.. list-table:: Spring and Damping Models used in Lethe-DEM.
   :widths: 40 30 30
   :header-rows: 1

   * - Parameters
     - Linear model definitions
     - Nonlinear viscoelastic model definitions
   * - Normal spring constant
     - :math:`k_n = \frac{16}{15}\sqrt{R_{e}}Y_{e}\left(\frac{15m_{e}V^2}{16\sqrt{R_{e}}Y_{e}}\right)^{0.2}`
     - :math:`k_n = \frac{4}{3}Y_{e}\sqrt{R_{e}\delta_n}`
   * - Tangential spring constant
     - :math:`\eta_n = \sqrt{\frac{4m_{e}k_n}{1+\left(\frac{\pi}{\ln{e}}\right)^2}}`
     - :math:`\eta_n = -2\sqrt{\frac{5}{6}}\beta\sqrt{S_nm_{e}}`
   * - Normal damping model constant
     - :math:`k_t = k_n`
     - :math:`k_t = 8G_{e}\sqrt{R_{e}\delta_n}`
   * - Tangential damping model constant
     - :math:`\eta_t = \eta_n`
     - :math:`\eta_t = -2\sqrt{\frac{5}{6}}\beta\sqrt{S_tm_{e}}`

Where:

* :math:`R_e` effective radius;
* :math:`Y_e` effective Young's modulus;
* :math:`m_e` effective mass;
* :math:`V` characteristic impact velocity;
* :math:`e` coefficient of restitution;
* :math:`G_e` effective shear modulus;

The parameters are computed as followed:

.. math::
    \frac{1}{m_{e}} &= \frac{1}{m_i}+\frac{1}{m_j} \\
    \frac{1}{R_{e}} &= \frac{1}{R_i}+\frac{1}{R_j} \\
    \frac{1}{G_{e}} &= \frac{2(2-\nu_i)(1+\nu_i)}{Y_i}+\frac{2(2-\nu_j)(1+\nu_j)}{Y_j} \\
    \frac{1}{Y_{e}} &= \frac{\left(1-\nu_i^2\right)}{Y_i}+\frac{\left(1-\nu_j^2\right)}{Y_j} \\
    \beta &= \frac{\ln{e}}{\sqrt{\ln^2{e}+\pi^2}} \\
    S_n &= 2Y_{e}\sqrt{R_{e}\delta_n} \\
    S_t &= 8G_{e}\sqrt{R_{e}\delta_n}

Where:

* :math:`\nu_i, \nu_j` poisson coefficient of particle i or j;

Rolling friction may be calculated through a constant torque model or a viscous torque model corresponding to those equations:

.. math::
    \mathbf{M}_{ij}^{r} &= -\mu_{r}R_{e}|\mathbf{F}_{ij}^{n}| \mathbf{\hat{\omega}}_{ij} \\
    \mathbf{M}_{ij}^{r} &= -\mu_{r}R_{e}|\mathbf{F}_{ij}^{n}||\mathbf{V}_{\omega}| \mathbf{\hat{\omega}}_{ij}

Where the parameters are:

.. math::
    \mathbf{\hat{\omega}}_{ij} &= \frac{\omega_{i} - \omega_{j}}{|\omega_{i} - \omega_{j}|} \\
    \mathbf{V}_{\omega} &= \left( \omega_{i} \times R_{i}\mathbf{n}_{ij}-\omega_{j} \times R_{j}\mathbf{n}_{ji} \right)

Where:

* :math:`\mu_{r}` rolling friction coefficient;

Tangential torque is calculated through:

.. math::
    \mathbf{M}_{ij}^{t} = R_{i}\mathbf{n}_{ij} \times \mathbf{F}_{ij}^{c}

Coulomb's criterion is violated when this condition is not respected during a collision:

.. math::
    |\mathbf{F}_{ij}^{t}| \geq \mu |\mathbf{F}_{ij}^{n}|


A violation means the collision is having gross sliding and tangential force needs to be limited to the Coulomb limit.
To do so, the tangential overlap :math:`\mathbf{\delta_t}` is first limited and then the tangential force is recalculated.

The tangential overlap is calculated with the tangential force with no damping force as default nonlinear contact model as follow:

.. math::
    \mathbf{\delta_t} &= \frac{\mathbf{\tilde{F}_{ij}}}{-k_{t}} \\
    \mathbf{\tilde{F}_{ij}} &= \mathbf{\hat{F}_{ij}} + \eta_{t}\mathbf{v}_{rt} \\
    \mathbf{\hat{F}_{ij}^{t}} &= \mu |\mathbf{F}_{ij}^{n}| \frac{\mathbf{F}_{ij}^{t}}{|\mathbf{F}_{ij}^{t}|}

Regarding the particle-wall contacts, applied models are the same than particle-particle contacts with a background triangulation and mapping with walls.

Integration methods
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



References
-------------
`[1] <https://doi.org/10.1007/s40571-022-00478-6>`_ Golshan et al. "Lethe-DEM: An open-source parallel discrete element solver with load balancing." Computational Particle Mechanics (2022) p.1-20

