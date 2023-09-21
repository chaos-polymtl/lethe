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
    set calculate average velocities     = false
    set initial time                     = 0.0

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
    # Tracer statistics
    set calculate tracer statistics      = false
    set tracer statistics name           = tracer_statistics

    # Thermal postprocesses
    set postprocessed fluid              = both
    set calculate temperature statistics = false
    set temperature statistics name      = temperature_statistics
    set calculate heat flux              = false
    set heat flux name                   = heat_flux

    # Multiphase postprocessing
    set calculate barycenter             = false
    set barycenter name                  = barycenter_information

    # Other Cahn-Hilliard postprocessing
    set calculate phase statistics       = false
    set phase statistics name            = phase_statistics
    
  end

* ``verbosity``: enables the display of the post-processing values in the terminal. This does not affect the printing of output files. Choices are: ``quiet`` (default, no output) or ``verbose`` (output at every iteration).

* ``output frequency``: frequency at which the enabled post-processing is outputted in the respective file. For ``output frequency = 1`` (default value), results are outputted at each iteration.

* ``calculate kinetic energy``: controls if calculation of total kinetic energy is enabled. 
    * ``kinetic energy name``: output filename for kinetic energy calculations.

* ``calculate average velocities``: controls if calculation of time-averaged velocities is enabled.
    * ``initial time``: initial time used for the average velocities calculations.

* ``calculate pressure drop``: controls if calculation of the pressure drop from the inlet boundary to the outlet boundary is enabled.
    * ``inlet boundary id`` and ``outlet boundary id``: define the IDs for inlet and outlet boundaries, respectively. 
    * ``pressure drop name``: output filename for pressure drop calculations.
    * The pressure drop :math:`\Delta P` and total pressure drop :math:`\Delta P_\text{total}` are calculated as such, with :math:`\Gamma` representing the boundary, :math:`\pmb{u}` the velocity  and :math:`P` the pressure:

.. math::
    \Delta P =  \frac{ \int_{\Gamma_\text{inlet}} P d \Gamma}{\int_{\Gamma_\text{inlet}} 1 d \Gamma} - \frac{ \int_{\Gamma_\text{outlet}} P d \Gamma}{\int_{\Gamma_\text{outlet}} 1 d \Gamma}

.. math::
    \Delta P_\text{total} =  \frac{ \int_{\Gamma_\text{inlet}} (P + \frac{1}{2} \pmb{u} \cdot \pmb{u}) d \Gamma}{\int_{\Gamma_\text{inlet}} d \Gamma} - \frac{ \int_{\Gamma_\text{outlet}} (P + \frac{1}{2} \pmb{u} \cdot \pmb{u}) d \Gamma}{\int_{\Gamma_\text{outlet}} d \Gamma}

* ``calculate flow rate``: controls if calculation of the volumetric flow rates at every boundary is enabled.
    * ``flow rate name``: output filename for flow rate calculations.
    * The flow rate :math:`Q` is calculated as such, with :math:`\Gamma` representing the boundary, :math:`\pmb{u}` the velocity and :math:`\pmb{n}` the vector normal to the surface:

.. math::
    Q =  \int_{\Gamma} \pmb{n} \cdot \pmb{u} d \Gamma

* ``calculate enstrophy``: controls if calculation of total enstrophy, which corresponds to dissipation effects in the fluid, is enabled. 
    * ``enstrophy name``: output filename for enstrophy calculations.

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
  
* ``calculate phase statistics``: outputs phase statistics from the solution of the Cahn-Hilliard equations, including minimum, maximum, average, and standard deviation of the phase order parameter. This works only with the :doc:`cahn_hilliard` solver.

* ``phase statistics name``: name of the output file containing phase order parameter statistics from Cahn-Hilliard simulations. The default file name is ``phase_statistics``.

        
