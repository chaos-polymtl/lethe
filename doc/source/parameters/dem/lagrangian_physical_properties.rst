==============================
Lagrangian Physical Properties
==============================

In this subsection, gravitational acceleration, and the physical properties of the particles and walls are defined. These properties include ``number of particle types``, and for each type, particle ``diameter``, particle ``density``, ``Young's modulus`` of particle and wall, ``Poisson ratio`` of particle and wall, ``restitution coefficient`` of particle and wall, ``friction coefficient`` of particle and wall and ``rolling friction coefficient`` of particle and wall.

.. code-block:: text

  subsection lagrangian physical properties
    # Gravitational acceleration vector
    set g                       = 0.0, 0.0, 0.0

    # Number of particle types
    set number of particle types = 1

    # Entering particle type 0
    subsection particle type 0

      # Choices are uniform, normal, lognormal or custom
      set size distribution type            = uniform

      # If distribution type = uniform
      set diameter                          = 0.001

      # If distribution type = custom
      set custom diameters                  = 0.001 , 0.0005
      set custom volume fractions           = 0.6   , 0.4

      # If distribution type = normal or lognormal
      set average diameter                  = 0.001
      set standard deviation                = 0.0
      set minimum diameter cutoff           = -1.
      set maximum diameter cutoff           = -1.

      # If distribution type = normal, lognormal or custom
      set distribution prn seed             = 1

      # For every distribution types
      set number of particles               = 0
      set density particles                 = 1000
      set young modulus particles           = 1000000
      set poisson ratio particles           = 0.3
      set restitution coefficient particles = 0.1
      set friction coefficient particles    = 0.1
      set rolling friction particles        = 0.1
      set rolling viscous damping particles = 0.1
      set surface energy particles          = 0.0
      set Hamaker constant particles        = 4.e-19
      set thermal conductivity particles    = 1
      set specific heat particles           = 1000
      set microhardness particles           = 1.e9
      set surface slope particles           = 0.1
      set surface roughness particles       = 1.e-9
      set thermal accommodation particles   = 0.7
      set real young modulus wall           = 0.
    end

    # Wall properties
    set young modulus wall           = 1000000
    set poisson ratio wall           = 0.3
    set restitution coefficient wall = 0.1
    set friction coefficient wall    = 0.1
    set rolling friction wall        = 0.1
    set rolling viscous damping wall = 0.1
    set surface energy wall          = 0.0
    set Hamaker constant wall        = 4.e-19
    set thermal conductivity wall    = 100
    set microhardness wall           = 1.e9
    set surface slope wall           = 0.1
    set surface roughness wall       = 1.e-10
    set thermal accommodation wall   = 0.7
    set real young modulus wall      = 0.

    # Interstitial gas properties
    set thermal conductivity gas     = 0.01
    set specific heat gas            = 1000
    set dynamic viscosity gas        = 1.e-5
    set specific heats ratio gas     = 1
    set molecular mean free path gas = 68.e-9
  end

* The ``g`` parameter defines the gravitational acceleration in `x`, `y`, and `z` directions. The deprecated version of this parameter is the 3 parameters ``gx``, ``gy``, and ``gz``.

* The ``number of particle types`` parameter specifies the number of particle types in a simulation. Particles with different sizes, size distributions, and physical properties have to be defined as separate particle types.

* For each particle type, we have to define a separate subsection (for instance, ``subsection particle type 0``) to specify its physical properties.

.. note::
    If the particles in a simulation are monodispersed and have the same physical properties, the ``number of particle types`` should be equal to zero. For polydispersed systems, the ``number of particle types`` is selected equal to the number of particles types in the simulation. For each particle type, a separate subsection ``particle type n`` should be defined (n starts from zero to ``number of particle types`` - 1) which contains all the physical properties related to that particle type.

* The ``size distribution type`` parameter specifies the size distribution for each particle type. For each particle type, four ``size distribution type`` can be defined: ``uniform``, ``normal``, ``lognormal`` and ``custom``.

  - For the ``uniform`` size distribution, the diameter of the particles is constant.
  - For the ``normal`` size distribution, the particle diameters are randomly sampled from a normal distribution with an average diameter and a standard deviation.
  - For the ``lognormal`` size distribution, the particle diameters are randomly sampled from a lognormal distribution with an average diameter and a standard deviation.
  - For the ``custom`` size distribution, particle diameters are sampled from a list of diameters with a corresponding list of probabilities.

.. note::
    In the ``custom`` size distribution, the probability values are based on the volume fraction taken by all the particles of the associated diameter, not to the total number of particles. For example, if a probability is equal to ``0.5`` , this means that half of the total volume of inserted particles will be occupied by particle with the associated diameter value.

* The ``diameter`` parameter defines the diameter of the particles in a ``uniform`` distribution.

* For a ``normal`` distribution, the ``average diameter`` and the ``standard deviation`` parameters defines the average (:math:`{\mu_d}`) and the standard deviation (:math:`{\sigma_d}`) of the particle size distribution *weighted by number*. The ``minimum cutoff`` and ``maximum cutoff`` parameters can be used to limit the lower and upper value of the diameter sampled from the distribution. If set to ``-1``, those two bounds are set to :math:`\mu_d \pm2.5 \sigma_d`. The number-weighted probability density function (:math:`f_{x}^{N}(d)`)  of a normal distribution is defined as:

.. math::
    f_{x}^{N}(d)=\frac{1}{\sigma_d\sqrt{2\pi}}\exp\left(-\frac{(d-\mu_d)^2}{2\sigma_{d}^2}\right).

* For a ``lognormal`` distribution, the ``average diameter`` and the ``standard deviation`` parameters defines the average (:math:`{\mu_d}`) and the standard deviation (:math:`{\sigma_d}`) of the particle size distribution *weighted by number*. The ``minimum cutoff`` and ``maximum cutoff`` parameters can be used to limit the lower and upper value of the diameter sampled from the distribution. If set to ``-1``, those two bounds are set to :math:`\exp(\mu \pm2.5 \sigma)`. The number-weighted probability density function (:math:`f_{x}^{N}(d)`)  of a lognormal distribution is defined as:

.. math::
    f_{x}^{N}(d)=\frac{1}{\sigma\sqrt{2\pi}}\exp\left(-\frac{(d-\mu)^2}{2\sigma^2}\right).

where :math:`{\mu}` and :math:`{\sigma}` are the mean and standard deviation of the underlying normal distribution which can be computed using:

.. math::
    \sigma = \sqrt{\ln \left(1 + \left(\frac{\sigma_d}{\mu_d}\right)^2\right)},

.. math::
    \mu = \ln(\mu_d) - 0.5 \sigma^2.

* For a ``custom`` distribution, the ``custom diameters`` parameter defines the different diameter values used when generating particles. The ``custom volume fractions`` parameter defines the probabilities corresponding to each diameter value previously declared based on volume fraction. Both list must have the same length.

* For a ``normal`` or a ``custom`` distribution, the ``distribution prn seed`` parameter defines the pseudo-random number (PRN) generator with which the diameters values are getting generated.

* The ``number of particles`` parameter defines the number of particles for each type.

* The ``density particles`` defines the density of particles for each type.

* The ``young modulus particles`` defines the Young's modulus for particles in each type.

* The ``poisson ratio particles`` defines the Poisson's ratio for particles in each type.

* The ``restitution coefficient particles`` defines the restitution coefficient for particles in each type.

* The ``friction coefficient particles`` defines the friction coefficient for particles in each type.

* The ``rolling friction particles`` defines the rolling friction coefficient of particles for each type.

* The ``rolling viscous damping particles``` defines the rolling viscous damping coefficient of the particles for the elasto-plastic spring-dashpot rolling friction model.

* The ``surface energy particles`` defines the surface energy of particles for each type. This parameter is used with the JKR and DMT force model.

* The ``Hamaker constant particles`` defines the Hamaker constant of particles for each type. This parameter is used with the DMT force model.

* The ``young modulus wall`` defines the Young's modulus of the walls.

* The ``poisson ratio wall`` defines the Poisson's ratio of the walls.

* The ``restitution coefficient wall`` defines the restitution coefficient of the walls.

* The ``friction coefficient wall`` defines the friction coefficient of the walls.

* The ``rolling friction wall`` defines the rolling friction coefficient of the walls.

* The ``rolling viscous damping wall`` defines the rolling viscous damping coefficient of the walls for the elasto-plastic spring-dashpot rolling friction model.

* The ``surface energy wall`` defines the surface energy of the walls. This parameter is used with the JKR and DMT force model.

* The ``Hamaker constant wall`` defines the Hamaker constant of the walls. This parameter is used with the DMT force model.

.. note::
    The following DEM parameters are used for multiphysic DEM simulations. All parameters should be specified in a consistent set of units (ideally SI).

* The ``thermal conductivity particles`` defines the thermal conductivity of particles for each type.

* The ``specific heat particles`` defines the specific heat of particles for each type.

* The ``microhardness particles`` defines the microhardness of particles for each type.

* The ``surface slope particles`` defines the surface slope of particles for each type. It is a non-dimensional parameter related to roughness and more precisely to the angle of the asperities on the surface. A higher surface slope entails a smaller microcontact resistance, as there are more microcontacts.

* The ``surface roughness particles`` defines the surface roughness of particles for each type.

* The ``thermal accommodation particles`` defines the thermal accommodation coefficient of particles for each type. The thermal accommodation coefficient characterizes the quality of thermal energy exchange between gas molecules and a solid surface.

* The ``real young modulus particles`` defines the real Young's modulus of particles for each type. It is used in multiphysic DEM to correct the thermal contact radius. This is useful when the Young's modulus in the simulation is lowered to increase the time-step. An artificially low Young's modulus would lead to an overestimated thermal contact radius. The real Young's modulus can only be given a value that is higher than the Young's modulus. Otherwise, the regular Young's modulus will be used in the calculations.

* The ``thermal conductivity gas`` defines the thermal conductivity of the interstitial gas.

* The ``specific heat gas`` defines the specific heat capacity of the interstitial gas.

* The ``dynamic viscosity gas`` defines the dynamic viscosity of the interstitial gas.

* The ``specific heats ratio gas`` defines the specific heats ratio of the interstitial gas.

* The ``molecular mean free path gas`` defines the molecular mean free path of the interstitial gas. It is the average distance a gas molecule will travel between collisions with other gas molecules.

* The ``thermal conductivity wall`` defines the thermal conductivity of the wall.

* The ``microhardness wall`` defines the microhardness of the wall.

* The ``surface slope wall`` defines the surface slope of the wall.

* The ``surface roughness wall`` defines the surface roughness of the wall.

* The ``thermal accommodation wall`` defines the thermal accommodation coefficient of the wall.

* The ``real young modulus wall`` defines the real Young's modulus of the wall.
