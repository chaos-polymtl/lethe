================
Model Parameters
================

In this subsection, contact detection, force models, time integration, load balancing and dynamic contact disabling parameters are defined.

.. code-block:: text

  subsection model parameters
    subsection contact detection
      # Contact detection method
      # Choices are constant|dynamic
      set contact detection method                = dynamic

      # Particle-particle contact neighborhood size
      set neighborhood threshold                  = 1.3

      set dynamic contact search size coefficient = 0.8
      set frequency                               = 10
    end

    subsection load balancing
      # Choices are none|once|frequent|dynamic|dynamic_with_disabling_contacts
      set load balance method     = none
      set particle weight         = 10000  # Every method, except none
      set step                    = 100000 # if method = once
      set frequency               = 100000 # if method = frequent
      set dynamic check frequency = 10000  # if method = dynamic
      set threshold               = 0.5    # if method = dynamic
    end

    # Particle-particle contact force model
    # Choices are linear|hertz_mindlin_limit_overlap|hertz_mindlin_limit_force|hertz
    set particle particle contact force method = hertz_mindlin_limit_overlap

    # Integration method
    # Choices are euler|velocity_verlet
    set integration method                     = velocity_verlet

    # Integration method
    # Choices are euler|velocity_verlet
    set integration method                     = velocity_verlet

    # Rolling resistance method
    # Choices are no_resistance|constant_resistance|viscous_resistance
    set rolling resistance torque method       = constant_resistance

    subsection dynamic disabling contacts
      set enable dynamic disabling contacts = false
      set granular temperature threshold    = 1e-4
      set solid fraction threshold          = 0.4
    end
  end


--------------------
Contact Detection
--------------------

Particle-particle contact search is a costly operation in DEM simulations. Contact detection parameters must be optimized to ensure a fast and physically accurate DEM simulation.

-  ``neighborhood threshold``  defines the spherical region around each particle which is used to generate the contact list. This parameter should generally be set between 1.3 and 1.5. It must be larger than 1 for contacts to be adequately taken into account.

Lethe defines two contact detection methods: ``dynamic`` and ``constant``

``contact detection method = dynamic``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lethe rebuilds the contact lists automatically. In this mode, Lethe stores the displacements of each particle in the simulation since the last contact detection. If the maximum displacement of a particle exceeds the smallest contact search criterion, then the iteration is a contact search iteration and the contact list is rebuilt. The smallest contact search criterion is the minimum of the smallest cell size in the triangulation or the radius of the spherical region in fine search, and it is defined as:
 
.. math::
  \phi=\min({d_c^{min}-r_p^{max},\epsilon(\alpha-1)r_p^{max}})

where :math:`{\phi}`, :math:`{d_c^{min}}`, :math:`{r_p^{max}}`, :math:`{\epsilon}`, and :math:`{\alpha}` denote smallest contact search criterion, minimum cell size (in the triangulation), maximum particle radius (in polydisperse simulations), ``dynamic contact search size coefficient``, and ``neighborhood threshold``.

* ``dynamic contact search size coefficient`` is a safety factor to ensure the late detection of particles will not happen in the simulations with ``dynamic`` contact search; and its value should be defined generally in the range of 0.5-1. 0.5 is a rather conservative value.


``contact detection method = constant``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Contact search will be carried out at constant frequency

* ``frequency`` is the frequency at which the contact list is renewed. It should be a value between 5 and 50 iterations. Small values of ``frequency`` lead to long simulation times, while large values of ``frequency`` may lead to late detection of collisions. Late detection of collisions can result in very large particles velocities (popcorn jump of particles in a simulation) or particles leaving the simulation domain.

-------------------------------
Contact and Integration Methods
-------------------------------

All contact force models are described in the :doc:`../../theory/dem/dem` section of the theory guide.


* ``integration`` controls the integration method  used. Lethe supports ``euler`` (1st order) and ``velocity-verlet`` (2nd order) time-integrators. The velocity-verlet should be used at all times. 

* ``particle particle contact force method`` control the particle-particle contact force model used. Four models are available in Lethe: ``hertz_mindlin_limit_overlap``, ``hertz_mindlin_limit_force``, ``hertz``, and ``linear``. 
  
* ``particle wall contact force method`` controls the particle-wall contact force model used. Two models are available: ``linear`` and ``non-linear``.

* ``rolling resistance method`` controls the rolling resistance model used. Three rolling resistance models are available: ``no_resistance``, ``constant_resistance``, ``viscous_resistance``


-----------------------
Load Balancing
-----------------------

Load-balancing updates the distribution of the subdomains between the processes in parallel simulation to achieve better computational performance (less simulation time). Three load-balancing methods are available in Lethe: ``once``, ``frequent``, or ``dynamic``. 

The total weight of each cell with particles in load-balancing is defined as:

.. math::
    W=1000+W_pn_p

where :math:`{W_p}` is the ``particle weight`` and :math:`{n_p}` is the number of particles in the cell. 1000 is the default weight assigned to one cell.

* ``particle weight`` must be defined for every ``load balance method``.

``load balance method = once``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Load balancing will be done only once.

* ``step`` the iteration number at which the load balancing will be carried out.

``load balance method = frequent``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Load balancing will be done at a given frequency

* ``frequency`` frequency (in iterations) of the load balancing.

``load balance method = dynamic``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Load balancing will be done when the computational load amongst core is too uneven. If 

.. math::
    L_{max}-L_{min}>{\beta}\bar{L}

load balancing will be executed. :math:`{L}` and :math:`{\beta}` denote computational load on a process and ``threshold``, respectively.

* ``dynamic check frequency`` frequency (in iterations) at which the load check on all processes is performed.
* ``threshold`` is the maximal load unbalance tolerated by the load balancing.

---------------------------
Dynamic Disabling Contacts
---------------------------

The dynamic disabling controls the disabling contact mechanism for performance enhancement. This feature dynamically searches for cells with low particle motion (granular temperature), disabling the computation of contacts for particles within these cells.

* ``enable dynamic disabling contacts`` enables the feature.

* ``granular temperature threshold`` is the threshold of the granular temperature below which the contacts are disabled.
* ``solid fraction threshold`` is the minimum solid fraction of the cell in which the contacts may be disabled.

Some parameters in the load balance section may be used to improve the performance of the dynamic disabling contacts feature using the dynamic load balancing.
Note: The ``load balance method`` may be set to ``dynamic_with_disabling_contacts`` and factors of the weight of the cells by mobility status may be adjusted using the ``active weight factor`` and ``inactive weight factor`` parameters. There is factor only for active and inactive status, mobile factor is always 1. 
