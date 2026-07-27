Time-Harmonic Maxwell
=====================

This subsection describes the parameters for the time-harmonic Maxwell solver, a
Discontinuous Petrov-Galerkin (DPG) formulation for time-harmonic electromagnetic
problems. These parameters are declared in the ``time harmonic maxwell`` subsection
of the parameter file.

.. attention::

  The time-harmonic Maxwell solver currently supports 3D problems only.

**Dimensionless formulation**

The solver always works internally in dimensionless form, and rescales the solution
back to physical units only at the end, according to the selected
``electromagnetic scaling type``. The system of equations is non-dimensionalized
using the following conventions:

.. math::

  \begin{align}
     & \tilde{\mathbf{E}} = \frac{1}{E_0}\mathbf{E}_\mathrm{spatial}, &                                                     & \tilde{\mathbf{H}} = \frac{Z_0}{E_0} \mathbf{H}_\mathrm{spatial},
     &                                                                & \varepsilon_r = \frac{1}{\varepsilon_0}\varepsilon, &                                                                                 & \mu_r = \frac{1}{\mu_0}\mu, \nonumber \\
     & \tilde{\nabla} = L \nabla,                                     &                                                     & \sigma_r = \frac{1}{\omega \varepsilon_0}\sigma,                                &
     & \tilde{\omega} =  \frac{ L}{c_0}\omega,                        &                                                     & \tilde{ \mathbf{J} } = \frac{L Z_0}{E_0} \mathbf{J}_\mathrm{spatial}. \nonumber
  \end{align}

where

* :math:`E_0` is the reference electric-field amplitude, set via
  ``electric field amplitude``,
* :math:`Z_0 = \sqrt{\mu_0/\varepsilon_0}` is the impedance of free space,
* :math:`L` is the reference length scale of the mesh,
* :math:`c_0` is the speed of light in vacuum,
* :math:`\varepsilon_0` and :math:`\mu_0` are the permittivity and permeability of
  vacuum,
* :math:`\varepsilon`, :math:`\mu`, and :math:`\sigma` are the (dimensional)
  permittivity, permeability, and electrical conductivity of the medium, and
* :math:`\omega` is the (dimensional) angular frequency set internally from the ``electromagnetic frequency`` parameter in Hz.

Using these conventions, the time-harmonic Maxwell system is written in
dimensionless form (using the :math:`e^{-i \omega t}` convention) as:

.. math::

  \begin{align}
    \nabla \times \tilde{\mathbf{E}} - i \tilde{\omega} \mu_r \tilde{\mathbf{H}}                      & = 0,                  \\
    \nabla \times \tilde{\mathbf{H}} + i \tilde{\omega} \varepsilon_{r,\text{eff}} \tilde{\mathbf{E}} & = \tilde{\mathbf{J}},
  \end{align}


where :math:`\tilde{\mathbf{E}}` and :math:`\tilde{\mathbf{H}}` are the
dimensionless electric and magnetic fields, :math:`\tilde{\omega}` is the
dimensionless angular frequency, :math:`\mu_r` is the relative magnetic
permeability, and :math:`\tilde{\mathbf{J}}` is the dimensionless current
density. 

.. note::

  The effective relative permittivity :math:`\varepsilon_{r,\text{eff}}`
  combines the medium's relative permittivity :math:`\varepsilon_r` with its
  conductive losses :math:`\sigma_r` into a single complex coefficient via the relation :math:`\varepsilon_{r,\text{eff}} = \varepsilon_r + i \sigma_r`.

.. tip::

  If the user wants to use different units than the MKS system because, for example, the problem is defined in centimeters, we recommend to use the dimensionality subsection :doc:`dimensionality <../parameters/dimensionality>` to set the reference length scale :math:`L` to 0.01 m, and leave all the parameters in the PRM file in MKS units. The solver will then automatically convert the parameters to the correct dimensionless form. If not using the dimensionality subsection, the user must ensure that all parameters are consistent with the desired behavior following the conventions described above that define the dimensionless formulation. 

**Parameter file syntax**

The parameters parsed by the solver are grouped as follows:

.. code-block:: text

  subsection time harmonic maxwell

    set electromagnetic frequency    = 1 # Hz
    
    set electromagnetic scaling type = <none|electric field|magnetic field|power>
    set electric field amplitude     = 1
    set magnetic field amplitude     = 1

    set number of waveguide inlets   = 0

    subsection waveguide inlet 0
      set port boundary id = 0
      set waveguide power  = 1 # Watts

      subsection waveguide mode
        set mode type    = <TE|TM>
        set mode order m = 1
        set mode order n = 0
      end

      set corner 0 = 0,0,0
      set corner 1 = 1,0,0
      set corner 2 = 0,1,0
      set corner 3 = 1,1,0
    end

    subsection time coupling strategy
      set type               = <none|iteration|time|threshold>
      set coupling iteration = 1
      set coupling time      = 1
      set coupling threshold = 0.1
    end
  end

.. note::
  The ``number of waveguide inlets`` is at most 10, and the inlets are numbered from 0 to 9. Their properties are defined in the ``waveguide inlet i`` subsections, where ``i`` is the inlet number. At the moment, all inlets need to have the same frequency, which is set in the parent ``time harmonic maxwell`` subsection as the solver does not support multiple frequencies solves at the same time.

* ``electromagnetic frequency``: excitation frequency of the time-harmonic
  wave, in Hz. This is required for any time-harmonic Maxwell problem.

* ``electromagnetic scaling type``: selects how the solution is rescaled back
  to physical units after the linear system is solved. The options are:

  * ``none``: the solution is left in dimensionless form and output as such. 
  * ``electric field``: rescaled so the electric field matches
    ``electric field amplitude`` and the magnetic field mathces :math:`H_0 = E_0 / Z_0`.
  * ``magnetic field``: rescaled so the magnetic field matches
    ``magnetic field amplitude``  and the electric field matches :math:`E_0 = H_0 * Z_0`.
  * ``power``: rescaled so the power delivered through the waveguide inlet
    matches the ``waveguide power`` given in watts by the user to this inlet. The scaling is done internally by isolating the electric field amplitude using the Poynting vector and the waveguide mode properties, and then rescaling the magnetic field accordingly. This is obtained via the following relation:
    .. math::

    \begin{equation}
      \overline{P}_\text{input} = \frac{1}{2} \frac{E_0^2}{Z_0} \int_A \Re{(\mathbf{E} \times \mathbf{H}^*)} \cdot \mathbf{n} \text{d}A.
    \end{equation}

    where :math:`\overline{P}_\text{input}` is the input  defined by the user, :math:`E_0` is the electric field amplitude, :math:`Z_0` is the impedance of free space, and the integral is taken over the waveguide inlet area :math:`A`. The Poynting vector is computed using the electric and magnetic fields of the selected waveguide mode.

  .. attention::

    Setting ``electromagnetic scaling type = power`` requires
    ``number of waveguide inlets`` to be at least 1, since there would otherwise
    be no waveguide power to normalize against.

  .. warning::

    The ``electromagnetic scaling type`` parameter is only used to rescale the
    solution back to physical units after the linear system is solved. It does
    not affect the non-dimensionalization of the equations, which is always done
    using the reference electric-field amplitude :math:`E_0` and the impedance
    of free space :math:`Z_0`. However, when performing multiphysics simulations, the scaling type must not be ``none`` if the electromagnetic fields are to be coupled to other physics, since the other physics will expect dimensional fields.

* ``electric field amplitude``: the reference amplitude :math:`E_0` used as the
    target amplitude when ``electromagnetic scaling type = electric field``.

* ``magnetic field amplitude``: dimensional magnetic-field amplitude used
  as the target amplitude when ``electromagnetic scaling type = magnetic field``.

* ``number of waveguide inlets``: number of ``waveguide inlet`` that needs to be 
    applied to the geometry boundary. Also see the :doc:`boundary conditions <../parameters/boundary_conditions_multiphysics>` section for more information on the waveguide boundary conditions.

.. note::
  At the moment only rectangular waveguide modes are supported, and the user must ensure that the selected mode is valid according to the waveguide geometry. An example of such a geometry is given in the :doc:`waveguide example <../../examples/multiphysics/waveguide>`.

* ``waveguide inlet i``: describes inlet ``i``. Each inlet defines:

  * ``port boundary id``: the boundary id of the mesh on which the inlet is applied.

  * ``waveguide power``: excitation power in Watts of the waveguide inlet. This is used to rescale the solution back to physical units when ``electromagnetic scaling type = power``.

  * ``waveguide mode``: the rectangular waveguide mode:

    * ``mode type``: ``TE`` or ``TM`` which stands for transverse electric or transverse magnetic, respectively. The ``TE`` mode is defined by the following equations:

      .. math::
        \begin{align}
            \mathbf{E}_{TE} & = i\frac{\omega \mu_r}{k_c^2}  \begin{bmatrix} -k_{x_2} \cos(k_{x_1} x_1) \sin(k_{x_2} x_2) \\ k_{x_1} \sin(k_{x_1} x_1) \cos(k_{x_2} x_2) \\ 0 \end{bmatrix} e^{i k_{x_3} x_3},                                    \\
            \mathbf{H}_{TE} & = \begin{bmatrix} - i \frac{k_{x_3} k_{x_1}}{k_c^2} \sin(k_{x_1} x_1) \cos(k_{x_2} x_2)   \\-i \frac{k_{x_3} k_{x_2}}{k_c^2} \cos(k_{x_1} x_1) \sin(k_{x_2} x_2)\\ \cos(k_{x_1} x_1) \cos(k_{x_2} x_2) \end{bmatrix} e^{i k_{x_3} x_3}.
        \end{align}

    The ``TM`` mode is defined by the following equations:
    
    .. math::
      \begin{align}
          \mathbf{E}_{TM} & = \begin{bmatrix} i \frac{k_{x_3} k_{x_1}}{k_c^2} \cos(k_{x_1} x_1) \sin(k_{x_2} x_2) \\  i \frac{k_{x_3} k_{x_2}}{k_c^2} \sin(k_{x_1} x_1) \cos(k_{x_2} x_2) \\ \sin(k_{x_1} x_1) \sin(k_{x_2} x_2) \end{bmatrix} e^{i k_{x_3} x_3}, \\
          \mathbf{H}_{TM} & = i\frac{\omega \epsilon_{r,eff}}{k_c^2} \begin{bmatrix} k_{x_2} \sin(k_{x_1} x_1) \cos(k_{x_2} x_2)    \\-k_{x_1} \cos(k_{x_1} x_1) \sin(k_{x_2} x_2)\\ 0 \end{bmatrix} e^{i k_{x_3} x_3}.
      \end{align}

      In the above equations, :math:`k_{x_1}` and :math:`k_{x_2}` are the wavenumbers in the two transverse directions, and :math:`k_{x_3}` is the wavenumber in the propagation direction. The cutoff wavenumber :math:`k_c` is defined as :math:`k_c^2 = k_{x_1}^2 + k_{x_2}^2`.
        

    * ``mode order m``, ``mode order n``: mode orders in the two transverse
      directions. They are used to compute the wavenumbers :math:`k_{x_1}` and :math:`k_{x_2}` in the above equations, which are given by :math:`k_{x_1} = m \pi / a` and :math:`k_{x_2} = n \pi / b`, where :math:`a` and :math:`b` are the waveguide dimensions in the two transverse directions.



  * ``corner 0`` through ``corner 3``: the four corners of the inlet, as 3D
    coordinates. They must be coplanar to define a valid rectangular waveguide inlet.

    .. important::

      The corners are used the following convention to associate the modes ``m`` and ``n`` to the correct transverse directions. In the reference frame of the inlet, the first transverse direction is defined by the vector from ``corner 0`` to ``corner 1``, and the second transverse direction is defined by the vector from ``corner 0`` to ``corner 2``. The mode orders ``m`` and ``n`` are associated with these two transverse directions, respectively. The following diagram illustrates this convention:

      .. image:: /images/corners_schematic.pdf
         :width: 300
         :alt: Corner Convention

* ``time coupling strategy``: controls how often the electromagnetic fields are
  recomputed during the time evolution of the coupled problem. 

  * ``type``: selects the coupling strategy:

    * ``none``: solved once before the first time step and never recomputed.
    * ``iteration``: recompute every ``coupling iteration`` time iterations.
    * ``time``: recompute every ``coupling time`` units of physical time.
    * ``threshold``: recompute whenever the electromagnetic properties of the
      medium change by more than ``coupling threshold``.

  * ``coupling iteration``: number of time iterations between two
    electromagnetic solves when ``type = iteration``.

  * ``coupling time``: physical time interval between two electromagnetic
    solves when ``type = time``.

  * ``coupling threshold``: maximum change in percent of across all the electromagnetic properties of the medium that triggers a new solve when ``type = threshold``. 

