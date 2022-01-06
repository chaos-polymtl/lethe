Model Parameters
-------------------
In this subsection, DEM simulation parameters are defined. These parameters include ``contact detection method`` and its subsequent information (``dynamic contact search size coefficient`` **or** ``contact detection frequency`` for ``dynamic`` **or** ``constant`` contact detection method), ``repartition method`` (``once``, ``frequent`` or ``dynamic``) for ``load balancing`` (in parallel simulations), ``particle weight`` which defines the particle weight compared to a default ``cell weight`` of 1000, ``neighborhood threshold`` (which defines the contact neighbor list size: ``neighborhood threshold`` * particle diameter), ``particle particle contact force method``, ``particle wall contact force method`` and ``integration method``. 

.. note::
    In ``constant`` contact search, Lethe-DEM searches for particle-particle and particle-wall collisions at a constant frequency (every ``contact detection frequency`` iterations). Normally, ``contact detection frequency`` should be selected in the range of 5-50. Small values of ``contact detection frequency`` lead to long simulation times, while large values of ``contact detection frequency`` may lead to late detection of collisions. Late detection of collisions can result in very large particles velocities (popcorn jump of particles in a simulation) or particles leaving the simulation domain.

.. note::
    In ``dynamic`` contact search, Lethe-DEM searches for particle-particle and particle-wall collisions automatically. In ``dynamic`` contact search, Lethe-DEM predicts (using the velocities of particles at each iteration) and stores (adds up) the displacements of each particle in the simulation. If the maximum displacement of particles exceeds the smallest contact search criterion (explained in the following), then the iteration is recognized as a contact detection iteration.

Smallest contact search criterion is defined as:
 
.. math::
    \phi=\min({d_c^{min}r_p^{max},\epsilon(\alpha-1)r_p^{max}})

where :math:`{\phi}`, :math:`{d_c^{min}}`, :math:`{r_p^{max}}`, :math:`{\epsilon}`, and :math:`{\alpha}` denote smallest contact search criterion, minimum cell size (in the triangulation), maximum particle radius (in polydisperse simulations), ``dynamic contact search size coefficient``, and ``neighborhood threshold``.

* ``dynamic contact search size coefficient``, as illustrated in the equation above, is a safety factor to ensure the late detection of particles will not happen in the simulations with ``dynamic`` contact search; and its value should be defined generally in the range of 0.5-1. 0.5 is a rather conservative value for ``dynamic contact search size coefficient``.


To optimize the computational performance of Lethe-DEM, the contact search is performed in two steps (broad and fine contact searches). In broad contact search, neighbor particles in adjacent cells are stored in a list. Then this list of adjacent particles is refined once more in the fine search. More information about the contact search algorithm in Lethe-DEM can be found `here <https://arxiv.org/abs/2106.09576>`_ . In the fine search, the particles (from the output list of broad search) which are located in a spherical region around each particle are stored in another list and sent to force calculation class. The diameter of this spherical region :math:`{D}` around each particle is defined as:

.. math::
    D={\alpha}d_p^{max}

where :math:`{d_p^{max}}` denotes the maximum particle size in the simulation.

* According to the definition of the spherical region in the fine search and contact force (equation above), the value of the ``neighborhood threshold`` :math:`{\alpha}` must be larger than 1. It is generally defined in the range of 1.3-1.5.

.. note::
    Load-balancing updates the distribution of the subdomains between the processes in parallel simulation to achieve better computational performance (less simulation time). Three load-balancing methods are available in Lethe-DEM: ``once``, ``frequent``, or ``dynamic``. Read here for more information about different load-balancing methods and their performances in various types of DEM simulations. The total weight of each cell with particles in load-balancing is defined as:
.. math::
    W=1000+W_pn_p

where :math:`{W_p}` is the ``particle weight`` and :math:`{n_p}` is the number of particles in the cell. 1000 is the default weight assigned to one cell.

Selecting ``repartition method = once``, requires defining the step at which the code calls load balancing (``load balance step``). ``Dynamic`` ``repartition method`` requires defining ``load balance frequency``, and in ``dynamic`` ``repartition method``, we should define ``load balance threshold`` and ``dynamic load balance check frequency``. In ``dynamic`` load balancing, the code checks the distribution of particles among the processors, every ``dynamic load balance check frequency`` steps, and if

.. math::
    L_{max}-L_{min}>{\beta}\bar{L}

it calls load-balancing. :math:`{L}` and :math:`{\beta}` denote computational load on a process and ``load balance threshold``, respectively.

.. note::
    Four ``particle particle contact models`` are available in Lethe-DEM (``hertz_mindlin_limit_overlap``, ``hertz_mindlin_limit_force``, ``hertz``, and ``linear``). ``hertz_mindlin_limit_overlap`` and ``hertz_mindlin_limit_force`` are non-linear contact models in which the stiffness and damping forces in both normal and tangential directions are considered. The only difference between these models is in their limiting method of the tangential force during gross sliding (where the tangential force exceeds the coulomb's limit). In ``hertz_mindlin_limit_overlap`` model, Lethe-DEM limits the tangential overlap and with limiting the overlap, the tangential force is limited; while in ``hertz_mindlin_limit_force`` model, the tangential force is limited directly without limiting the tangential overlap. ``hertz`` is another non-linear model in which the damping force is not considered in the tangential direction and the tangential force is limited in gross sliding. Lethe-DEM also has a ``linear`` contact model (the stiffness and damping forces are linear functions of overlap and relative velocity, respectively).

.. note::
    Lethe-DEM has two linear and non-linear particle-wall models.

.. note::
    ``euler`` (1st order) and ``velocity-verlet`` (2nd order) are the available integration methods in Lethe-DEM.

.. note::
    Three rolling resistance models are available in Lethe-DEM: ``no_resistance``, ``constant_resistance``, ``viscous_resistance``.


.. code-block:: text

 subsection model parameters
  # Contact detection method
  # Choices are constant|dynamic
  set contact detection method                          = dynamic

  # Depending on the contact detection method, contact search size coefficient (safety factor multiplier for dynamic contact search) or contact detection frequency should be defined for dynamic and constant contact search methods, respectively.
  set dynamic contact search size coefficient           = 0.9
  set contact detection frequency                       = 20

  # Load balancing method and its subsequent information
  set load balance method                               = frequent
  set load balance frequency                            = 200000
  set load balance particle weight                      = 10000

  # Particle-particle contact neighborhood size
  set neighborhood threshold                            = 1.6

  # Particle-particle contact force model
  # Choices are linear|hertz_mindlin_limit_overlap|hertz_mindlin_limit_force|hertz
  set particle particle contact force method            = hertz_mindlin_limit_overlap

  # Particle-wall contact force model
  # Choices are linear|nonlinear
  set particle wall contact force method                = nonlinear

  # Integration method
  # Choices are euler|velocity_verlet
  set integration method                                = velocity_verlet

  # Rolling resistance method
  # choices are no_resistance|constant_resistance|viscous_resistance
  set rolling resistance torque method                  = no_resistance
 end

