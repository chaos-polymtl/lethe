===============
Post-processing
===============

This subsection controls the post-processing other than the forces and torque on the boundary conditions. Default values are

.. code-block:: text

  subsection post-processing
    set verbosity                        = quiet
    set output frequency                 = 1

    #---------------------------------------------------
    # Fluid dynamic post-processing
    #---------------------------------------------------
    # Kinetic energy calculation
    set calculate kinetic energy         = false
    set kinetic energy name              = kinetic_energy

    # Average velocities calculation
    set calculate average velocities      = false
    set initial time for average velocity = 0.0

    # Average temperature calculation
    set calculate average temperature and heat flux        = false
    set initial time for average temperature and heat flux = 0.0

    # Pressure drop calculation
    set calculate pressure drop          = false
    set pressure drop name               = pressure_drop
    set inlet boundary id                = 0
    set outlet boundary id               = 1

    # Flow rate at boundaries calculation
    set calculate flow rate              = false
    set flow rate name                   = flow_rate

    # Enstrophy calculation
    set calculate enstrophy              = false
    set enstrophy name                   = enstrophy

    # Viscous dissipation
    set calculate viscous dissipation    = false
    set viscous dissipation name         = viscous_dissipation

    # Pressure power
    set calculate pressure power         = false
    set pressure power name              = pressure_power

    # Others
    set smoothed output fields           = false

    #---------------------------------------------------
    # Physical properties post-processing
    #---------------------------------------------------
    set calculate apparent viscosity     = false
    set apparent viscosity name          = apparent_viscosity

    #---------------------------------------------------
    # Multiphysics post-processing
    #---------------------------------------------------
    # Tracer postprocessing
    set calculate tracer statistics      = false
    set tracer statistics name           = tracer_statistics
    set calculate tracer flow rate       = false
    set tracer flow rate name            = tracer_flow_rate

    # Thermal postprocessing
    set postprocessed fluid              = both
    set calculate temperature statistics = false
    set temperature statistics name      = temperature_statistics
    set calculate heat flux              = false
    set heat flux name                   = heat_flux

    # Multiphase postprocessing
    set calculate barycenter             = false
    set barycenter name                  = barycenter_information
    set calculate mass conservation      = true
    set mass conservation name           = mass_conservation_information

    # Other Cahn-Hilliard postprocessing
    set calculate phase statistics       = false
    set phase statistics name            = phase_statistics
    set calculate phase energy           = false
    set phase energy name                = phase_energy

    #---------------------------------------------------
    # Multiphase post-processing
    #---------------------------------------------------
    # CFD-DEM postprocessing
    set calculate volume phases          = false
    set phase volumes name               = phase_volumes
    
  end

* ``verbosity``: enables the display of the post-processing values in the terminal. This does not affect the printing of output files. Choices are: ``quiet`` (default, no output) or ``verbose`` (output at every iteration).

* ``output frequency``: frequency at which the enabled post-processing is outputted in the respective file. For ``output frequency = 1`` (default value), results are outputted at each iteration.

* ``calculate kinetic energy``: controls if calculation of kinetic energy is enabled. 
    * ``kinetic energy name``: output filename for kinetic energy calculations.
    * The kinetic energy :math:`{E}_k` is calculated as 

    .. math::
      {E}_k =  \frac{1}{2 \Omega} \int_{\Omega} \mathbf{u} \cdot \mathbf{u} \ \mathrm{d} \Omega

    with :math:`\Omega` representing the volume of the domain and :math:`\mathbf{u}` the velocity.
    

* ``calculate average velocities``: controls if calculation of time-averaged velocities is enabled.
    * ``initial time for average velocity``: initial time used for the average velocities calculations.

* ``calculate average temperature and heat flux``: controls if calculation of time-averaged temperature and time-averaged heat flux is enabled.
    * ``initial time for average temperature``: initial time used for the average temperature calculations.

* ``calculate pressure drop``: controls if calculation of the pressure drop from the inlet boundary to the outlet boundary is enabled.
    * ``inlet boundary id`` and ``outlet boundary id``: define the IDs for inlet and outlet boundaries, respectively. 
    * ``pressure drop name``: output filename for pressure drop calculations.
    * The pressure drop :math:`\Delta p` and total pressure drop :math:`\Delta p_\text{total}` are calculated as:

    .. math::
      \Delta p =  \frac{ \int_{\Gamma_\text{inlet}} p \mathrm{d} \Gamma}{\int_{\Gamma_\text{inlet}} 1 \mathrm{d} \Gamma} - \frac{ \int_{\Gamma_\text{outlet}} p \mathrm{d} \Gamma}{\int_{\Gamma_\text{outlet}} 1 \mathrm{d} \Gamma}

    .. math::
      \Delta p_\text{total} =  \frac{ \int_{\Gamma_\text{inlet}} (p + \frac{1}{2} \mathbf{u} \cdot \mathbf{u}) \mathrm{d} \Gamma}{\int_{\Gamma_\text{inlet}} \mathrm{d} \Gamma} - \frac{ \int_{\Gamma_\text{outlet}} (p + \frac{1}{2} \mathbf{u} \cdot \mathbf{u}) \mathrm{d} \Gamma}{\int_{\Gamma_\text{outlet}} \mathrm{d} \Gamma}

    with :math:`\Gamma` representing the boundary, :math:`\mathbf{u}` the velocity  and :math:`p` the pressure.

* ``calculate flow rate``: controls if calculation of the volumetric flow rates at every boundary is enabled.
    * ``flow rate name``: output filename for flow rate calculations.
    * The flow rate :math:`Q` is calculated as such, with :math:`\Gamma` representing the boundary, :math:`\mathbf{u}` the velocity and :math:`\mathbf{n}` the vector normal to the surface:

.. math::
    Q =  \int_{\Gamma} \mathbf{n} \cdot \mathbf{u} d \Gamma

* ``calculate enstrophy``: controls if the volume-averaged enstrophy is calculated.
    * ``enstrophy name``: output filename for enstrophy calculations.
    * The enstrophy :math:`\mathcal{E}` is calculated as 

    .. math::
      \mathcal{E} =  \frac{1}{2 \Omega} \int_{\Omega} \mathbf{\omega} \cdot \mathbf{\omega} \mathrm{d} \Omega

    with :math:`\Omega` representing the volume of the domain and :math:`\mathbf{\omega}` the vorticity.

* ``calculate viscous dissipation``: controls if the viscous dissipation is calculated.
    * ``viscous dissipation name``: output filename for the viscous dissipation calculations.
    * The viscous dissipation is calculated as 

    .. math::
       \frac{1}{\Omega} \int_{\Omega} \mathbf{\tau} : \nabla\mathbf{u} \mathrm{d} \Omega

    with :math:`\Omega` representing the volume of the domain and :math:`\mathbf{\tau}` the deviatoric stress tensor.

* ``calculate pressure power``: controls if the pressure power is calculated.
    * ``pressure power name``: output filename for the pressure power calculations.
    * The pressure power is calculated as

    .. math::
       \frac{1}{\Omega} \int_{\Omega}  \nabla p \cdot \mathbf{u} \mathrm{d} \Omega

    with :math:`\Omega` representing the volume of the domain, :math:`\mathbf{u}` the velocity  and :math:`p` the pressure.

* ``smoothed output fields``: controls if the Qcriterion field will be smoothed using an L2-projection over the nodes. The same will shortly be applied to the Vorticity. 

* ``calculate apparent viscosity``: controls if parameter calculation of an apparent viscosity is enabled, when using a non Newtonian flow (see section Physical properties - :ref:`rheological_models`). This is mainly used to define the Reynolds number `a posteriori`. 
    * ``apparent viscosity name``: output filename for apparent viscosity calculations.

* ``calculate tracer statistics``: controls if calculation of tracer statistics is enabled. Statistics include: minimum, maximum, average and standard-deviation.
    .. warning::

        Do not forget to ``set tracer = true`` in the :doc:`multiphysics` subsection of the ``.prm``.

    * ``tracer statistics name``: output filename for tracer statistics calculations.

* ``postprocessed fluid``: fluid domain used for thermal postprocesses. Choices are : ``fluid 0``, ``fluid 1``, or ``both`` (default).
    * For monophasic simulations (``set VOF = false`` in :doc:`multiphysics`), ``both`` and ``fluid 0`` are equivalent and the temperature statistics are computed over the entire domain.
    * For multiphasic simulations (``set VOF = true`` in :doc:`multiphysics`), temperature statistics can be computed over the entire domain (``both``) or inside a given fluid only (``fluid 0`` or ``fluid 1``), with the fluid IDs defined in Physical properties - :ref:`two phase simulations`.

    .. note::

        The output files will have a suffix depending on the ``postprocessed fluid``: ``fluid_0``, ``fluid_1`` and ``all_domain``.

* ``calculate temperature statistics``: controls if calculation of temperature statistics is enabled. Statistics include: minimum, maximum, average and standard-deviation.

    * ``temperature statistics name``: output filename for temperature statistics calculations.

    .. admonition:: Example of temperature statistics table:

        .. code-block:: text

             time  min    max    average std-dev 
            0.0000 0.0000 3.9434  0.1515  0.6943 
            0.2000 2.5183 4.9390  3.3917  0.7229 

* ``calculate heat flux``: controls if calculation of heat flux is enabled. If enabled, these quantities are postprocessed: 

  1. the total heat flux :math:`q_{tot}` for each :ref:`heat transfer bc` boundary condition. The total heat flux on a boundary :math:`\Gamma` is defined as:

  .. math:: 

      q_\text{tot} = \int_\Gamma (\rho C_p \mathbf{u} \mathbf{T} - k \nabla \mathbf{T}) \cdot \mathbf{n}


  The output table is appended with one column per :ref:`heat transfer bc` boundary condition, named ``bc_i`` where ``i`` is the index of the boundary in the parameter file.

  2. the convective heat flux :math:`q_\text{conv}` for each :ref:`heat transfer bc` boundary condition. The convective heat flux on a boundary :math:`\Gamma` is defined as:

  .. math:: 

      q_\text{conv} = \int_\Gamma  h (\mathbf{T}-\mathbf{T}_\infty)

  The output table is appended with one column per :ref:`heat transfer bc` boundary condition, named ``bc_i`` where ``i`` is the index of the boundary in the parameter file.

  3. the thermal energy (:math:`\mathbf{Q} = m c_p \mathbf{T}`) over the domain defined by ``postprocessed fluid``. 

  4. if there is a :doc:`nitsche`, the total heat fluxes on each solid: :math:`q_\text{nitsche} = \beta_\text{heat} \left( \mathbf{T}_\text{nitsche} - \mathbf{T} \right)`

  The output table is appended with one column per solid, named ``nitsche_solid_i`` where ``i`` is the index of the ``nitsche solid`` in the parameter file.

  .. warning ::

      Do not forget to ``set enable heat boundary condition = true`` in the :doc:`nitsche` subsection of the ``.prm``.


  * ``heat flux name``: output filename for heat flux calculations.

    .. admonition:: Example of heat flux table:

        .. code-block:: text

		 time  total_flux_bc_0 convective_flux_bc_0 thermal_energy_fluid flux_nitsche_solid_0 
		0.0000          0.0000               0.0000               0.0000            1000.0000 
		1.0000         -0.9732               0.0000               1.4856               0.9732 

* ``calculate barycenter``: calculates the barycenter of ``fluid 1`` and its velocity in VOF and Cahn-Hilliard simulations. The barycenter :math:`\mathbf{x}_b` and its velocity :math:`\mathbf{v}_b` are defined as:

  .. math::

      \mathbf{x_b} = \frac{\int_{\Omega} \psi \mathbf{x} \mathrm{d}\Omega }{\int_{\Omega} \psi \mathrm{d}\Omega}

  .. math::

      \mathbf{v_b} = \frac{\int_{\Omega} \psi \mathbf{u} \mathrm{d}\Omega }{\int_{\Omega} \psi \mathrm{d}\Omega}

  where :math:`\psi \in [0,1]` is the filtered phase indicator for VOF simulations. 
  
  For Cahn-Hilliard the formula is slightly different since the phase order parameter :math:`\phi` belongs to the :math:`[-1,1]` interval:
  
  .. math::

      \mathbf{x_b} = \frac{\int_{\Omega} 0.5(1-\phi) \mathbf{x} \mathrm{d}\Omega }{\int_{\Omega} 0.5(1-\phi) \mathrm{d}\Omega}

  .. math::

      \mathbf{v_b} = \frac{\int_{\Omega} 0.5(1-\phi) \mathbf{u} \mathrm{d}\Omega }{\int_{\Omega} 0.5(1-\phi) \mathrm{d}\Omega}
      
  where :math:`\phi` is the phase order parameter.
  
  
* ``barycenter name``: name of the output file containing the position and velocity of the barycenter for VOF and Cahn-Hilliard simulations. The default file name is ``barycenter_information``.

* ``calculate mass conservation``: calculates the mass and momentum of both fluids for VOF simulations.

* ``mass conservation name``: name of the output file containing the mass of both fluids for VOF simulations. The default file name is ``mass_conservation_information``.
  
* ``calculate phase statistics``: outputs Cahn-Hilliard phase statistics, including minimum, maximum, average, integral of the phase order parameter, and the volume of each phase.

  .. warning ::

      ``calculate phase statistics = true`` only works with the :doc:`cahn_hilliard` solver.

* ``phase statistics name``: name of the output file containing phase order parameter statistics from Cahn-Hilliard simulations. The default file name is ``phase_statistics``. It is stored in the output folder with in a  ``.dat`` file.

* ``calculate phase energy``: outputs Cahn-Hilliard phase energies, including bulk energy, interface energy and total energy. The energies are computed as follow:

  .. math::

     E_{bulk} = \int_{\Omega} (1-\phi^2)^2 \mathrm{d}\Omega 
      
  .. math::

     E_{interface} = \int_{\Omega} 0.5\epsilon^2|\nabla \phi |^2 \mathrm{d}\Omega 
      
  .. math::

     E_{total} = E_{bulk} + E_{interface}  
    
  where :math:`\epsilon` is the numerical interface thickness. Note that these energies are not homogeneous to physical energies. Nonetheless, they are a convenient way to track the system's evolution.
  
  .. warning ::

      ``calculate phase energy = true`` only works with the :doc:`cahn_hilliard` solver.


* ``phase energy name``: name of the output file containing phase energies from Cahn-Hilliard simulations. The default file name is ``phase_energy``.

* ``calculate phase volumes``: outputs total volume of fluid phase and total volume of solid phase in CFD-DEM simulation. These volumes are computed as follow:

  .. math::

     V_{fluid} = \int_{\Omega} \varepsilon_f \mathrm{d}\Omega 
      
  .. math::

     V_{solid} = \int_{\Omega} (1 - \varepsilon_f) \mathrm{d}\Omega 
      
  where :math:`\varepsilon` is the void fraction.  This is a convenient way to check if the volume of each phase is conserved.
  
  .. warning ::

      ``calculate phase volumes = true`` only works with the ``lethe-fluid-particle`` solver.


* ``phase volumes name``: name of the output file containing phase energies from Cahn-Hilliard simulations. The default file name is ``phase_volumes``.

        
