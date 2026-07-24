==================================
Boundary Conditions - Multiphysics
==================================

This subsection's purpose is defining the boundary conditions associated to multiphysic problems. 

.. _heat transfer bc:

Heat Transfer
^^^^^^^^^^^^^

For heat transfer boundary conditions, the possible ``types`` are ``noflux`` (default), ``temperature`` and ``convection-radiation-flux``.
The default parameters for ``temperature`` and ``convection-radiation-flux`` are shown:

.. code-block:: text

  subsection boundary conditions heat transfer
    set number         = 3
    set time dependent = false

    subsection bc 0
      set id                 = 0
      set type               = periodic
      set periodic id        = 1
      set periodic direction = 0
    end
    subsection bc 1
      set id   = 1
      set type = temperature
      subsection value
        set Function expression = 0
      end
    end
    subsection bc 2
      set id   = 2
      set type = convection-radiation-flux
      subsection h
        set Function expression = 0
      end
      subsection Tinf
        set Function expression = 0
      end
      subsection emissivity
        set Function expression = 0
      end
      subsection heat_flux
        set Function expression = 0
      end
    end
    set Stefan-Boltzmann constant = 0.000000056703
  end


* ``number``: This is the number of boundary conditions of the problem. Periodicity between 2 boundaries counts as 1 condition even if it requires two distinct boundary ids.

* ``time dependent`` specifies if a boundary condition is time-dependent (``true``) or steady (``false``). By default, this parameter is set to ``false``. This is here to improve the computational efficiency for transient cases in which the boundary conditions do not change.

.. warning::
    The ``number`` of boundary conditions must be specified explicitly. This is often a source of error.

.. note::
    The index in ``subsection bc ..`` must be coherent with the ``number`` of boundary conditions set: if ``number = 2``, ``bc 0`` and ``bc 1`` are created but ``bc 2`` does not exist. 

    Likewise, if ``number = 2`` and there is no ``subsection bc 0`` explicitly stated, the boundary is still created, with ``noflux`` by default.

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number ``#`` of the bc. The parameter also accepts a list of boundary ids, to which the same boundary condition will be applied. In this case, id numbers should be separated by commas, for example: ``set id = 0,1,2``.

* ``type``: type of boundary condition being imposed. At the moment, choices are:
    * ``noflux`` (default) so that there is no heat transfer boundary condition,
    * ``temperature`` (Dirichlet BC), to impose a given temperature ``value`` at the boundary,
    * ``periodic`` to impose periodicity between boundaries. ``periodic id`` and ``periodic direction`` specify the id and direction of the matching periodic boundary condition. For example, if boundary id 0 (located at xmin) is matched with boundary id 1 (located at xmax), we would set ``id = 0``, ``periodic id = 1`` and ``periodic direction = 0``;
    * ``convection-radiation-flux`` (Robin BC) for cooling/heating, depending on the environment temperature at the boundary ``Tinf``, with a given heat transfer coefficient ``h`` and ``emissivity`` of the boundary :math:`\mathbf{\epsilon}` following Newton's law of cooling (and heating) and Stefan-Boltzmann law of radiation. It is also possible to impose a given heat flux (:math:`q_0`) by using the parameter ``heat_flux``. This BC can be represented by:

    .. math::
        \frac{ \partial T}{\partial \mathbf{n}} = h (T - T_{inf}) + \epsilon \sigma (T^4 - T_{inf}^4) + q_0

    where :math:`\mathbf{\sigma}` is the Stefan-Boltzmann constant.

    .. important::

      The flux represented by the ``convection-radiation-flux`` BC follow the direction of the normal vector to the boundary, i.e., pointing outwards the boundary. As consequence, a positive value for ``heat_flux``, for example, will result on heat being extracted from the boundary.

.. seealso::

  The :doc:`../../examples/multiphysics/warming-up-a-viscous-fluid/warming-up-a-viscous-fluid` example uses heat transfer boundary conditions.


Tracer
^^^^^^

For tracer boundary conditions, the defaults parameters are:

.. code-block:: text

  subsection boundary conditions tracer
    set number         = 2
    set time dependent = false
    subsection bc 0
      set id                 = 0
      set type               = periodic
      set periodic id        = 1
      set periodic direction = 0
    end
    subsection bc 1
      set id   = 1
      set type = dirichlet
      subsection dirichlet
        set Function expression = 0
      end
    end
  end

* ``number``: This is the number of boundary conditions of the problem. Periodicity between 2 boundaries counts as 1 condition even if it requires two distinct boundary ids.

* ``time dependent`` specifies if a boundary condition is time-dependent (``true``) or steady (``false``). By default, this parameter is set to ``false``. This improves the computational efficiency for transient cases in which the boundary conditions do not change.

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number ``#`` of the bc. The parameter also accepts a list of boundary ids, to which the same boundary condition will be applied. In this case, id numbers should be separated by commas, for example: ``set id = 0,1,2``.

* ``type``: This is the type of boundary condition being imposed:
    * ``dirichlet`` to impose specific values;
    * ``periodic`` to impose periodicity between boundaries. ``periodic id`` and ``periodic direction`` specify the id and direction of the matching periodic boundary condition. For example, if boundary id 0 (located at xmin) is matched with boundary id 1 (located at xmax), we would set ``id = 0``, ``periodic id = 1`` and ``periodic direction = 0``;
    * ``outlet`` to impose an outflow when using the tracer Discontinuous-Galerkin implementation;
    * ``none`` is a do-nothing boundary condition.

CLS
^^^

For CLS boundary conditions (multiphase flow), the possible ``types`` are ``none`` (default) and ``dirichlet``, as shown below.

.. code-block:: text

  subsection boundary conditions CLS
    set number         = 3
    set time dependent = false
    subsection bc 0
      set id                 = 0
      set type               = periodic
      set periodic id        = 1
      set periodic direction = 0
    end
    subsection bc 1
      set id   = 1
      set type = none
    end
    subsection bc 2
      set id   = 2
      set type = dirichlet
      subsection dirichlet
        set Function expression = 0
      end
    end
  end

.. warning::
    The ``number`` of boundary conditions must be specified explicitly. Periodicity between 2 boundaries counts as 1 condition even if it requires two distinct boundary ids.

.. note::
    The index in ``subsection bc ..`` must be coherent with the ``number`` of boundary conditions set: if ``number = 2``, ``bc 0`` and ``bc 1`` are created but ``bc 2`` does not exist. 

    Likewise, if ``number = 2`` and there is no ``subsection bc 0`` explicitly stated, the boundary is still created, with ``none`` by default.

* ``number``: This is the number of boundary conditions of the problem.

* ``time dependent`` specifies if a boundary condition is time-dependent (``true``) or steady (``false``). By default, this parameter is set to ``false``. This improves the computational efficiency for transient cases in which the boundary conditions do not change.

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number ``#`` of the bc. The parameter also accepts a list of boundary ids, to which the same boundary condition will be applied. In this case, id numbers should be separated by commas, for example: ``set id = 0,1,2``.

* ``type``: This is the type of boundary condition being imposed. At the moment, choices are:
    * ``none`` for which nothing happens;
    * ``dirichlet`` for inlet and outlet boundary conditions, to specify which fluid should be at the selected boundary;
    * ``periodic`` to impose periodicity between boundaries. ``periodic id`` and ``periodic direction`` specify the id and direction of the matching periodic boundary condition. For example, if boundary id 0 (located at xmin) is matched with boundary id 1 (located at xmax), we would set ``id = 0``, ``periodic id = 1`` and ``periodic direction = 0``.

    
Cahn-Hilliard
^^^^^^^^^^^^^^

For Cahn-Hilliard boundary conditions, the available ``types`` are ``none`` (default), ``noflux``, ``dirichlet``, ``angle_of_contact``, ``free_angle`` and ``periodic``. The parameters for each type of Cahn-Hilliard boundary conditions are:

.. code-block:: text

  subsection boundary conditions cahn hilliard
    set number         = 5
    set time dependent = false
    subsection bc 0
      set id                 = 0
      set type               = periodic
      set periodic id        = 1
      set periodic direction = 0
    end
    subsection bc 1
      set id   = 1
      set type = none
    end
    subsection bc 2
      set id   = 2
      set type = dirichlet
      subsection phi
        set Function expression = 0
      end
    end
    subsection bc 3
      set id          = 3
      set type        = angle_of_contact
      set angle value = 90 # The angle is given in degrees (°)
    end
    subsection bc 4
      set id   = 4
      set type = free_angle
    end
  end


* ``number``: This is the number of boundary conditions of the problem. Periodicity between 2 boundaries counts as 1 condition even if it requires two distinct boundary ids.

* ``time dependent`` specifies if a boundary condition is time-dependent (``true``) or not (``false``). By default, this parameter is set to ``false``. It is used to improve the computational efficiency of transient cases in which the boundary conditions do not change.

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number ``#`` of the bc. The parameter also accepts a list of boundary ids, to which the same boundary condition will be applied. In this case, id numbers should be separated by commas, for example: ``set id = 0,1,2``.

* ``type``: Type of boundary condition being imposed. At the moment, the choices are:
    * ``none`` (default): for which nothing happens.
    * ``noflux``: no phase leaves the simulation domain.
    * ``dirichlet``: Imposes a given phase order parameter function on the boundary. This function can depend on position (:math:`x,y,z`) and on time (:math:`t`).
    * ``angle_of_contact``: Imposes a given angle of contact ``angle value`` between the two phases at the boundary. It refers to the inner angle of contact, in degrees (°).
    * ``free_angle``: Leaves the angle as a free variable to be solved.
    * ``periodic``: Imposes periodicity between boundaries. ``periodic id`` and ``periodic direction`` specify the id and direction of the matching periodic boundary condition. For example, if boundary id 0 (located at xmin) is matched with boundary id 1 (located at xmax), we would set ``id = 0``, ``periodic id = 1`` and ``periodic direction = 0``.

Time-Harmonic Maxwell
^^^^^^^^^^^^^^^^^^^^^

For Time-Harmonic Maxwell boundary conditions, the possible ``types`` are ``pec`` (perfect electric conductor), ``pmc`` (perfect magnetic conductor), ``silver muller``, ``electric field``, ``magnetic field``, ``impedance boundary`` and ``waveguide port``. The default parameters are shown below.

.. code-block:: text

  subsection boundary conditions time harmonic maxwell
    set number         = 1
    set time dependent = false
    subsection bc 0
      set id   = 0
      set type = silver muller
    end
  end

* ``number``: This is the number of boundary conditions of the problem. Periodicity between 2 boundaries counts as 1 condition even if it requires two distinct boundary ids.

* ``time dependent`` specifies if a boundary condition is time-dependent (``true``) or steady (``false``). By default, this parameter is set to ``false``. This improves the computational efficiency for transient cases in which the boundary conditions do not change.

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number ``#`` of the bc. The parameter also accepts a list of boundary ids, to which the same boundary condition will be applied. In this case, id numbers should be separated by commas, for example: ``set id = 0,1,2``.

* ``type``: This is the type of boundary condition being imposed. At the moment, the choices are:
    * ``pec`` for a perfect electric conductor boundary. The tangential electric field is forced to vanish on the boundary,

      .. math::
        \mathbf{n} \times \mathbf{E} = 0.

      Example:

      .. code-block:: text

        subsection bc 0
          set id   = 0
          set type = pec
        end

    * ``pmc`` for a perfect magnetic conductor boundary. The tangential magnetic field is forced to vanish on the boundary,

      .. math::
        \mathbf{n} \times \mathbf{H} = 0.

      Example:

      .. code-block:: text

        subsection bc 0
          set id   = 0
          set type = pmc
        end

    * ``silver muller`` for a Silver-Müller absorbing boundary condition. In its standard form, the boundary relation is written as

      .. math::
        \mathbf{n} \times \mathbf{H} + \sqrt{\frac{\varepsilon_{r,\mathrm{eff}}}{\mu_r}} \mathbf{n} \times ( \mathbf{E} \times \mathbf{n} ) = 0.

      Example:

      .. code-block:: text

        subsection bc 0
          set id   = 0
          set type = silver muller
        end

    * ``electric field`` to impose a prescribed electric field. The boundary values are taken from the parsed functions in the ``E x/y/z real part`` and ``E x/y/z imag part`` subsections,

      .. math::
        \mathbf{E} = \mathbf{E}^{\mathrm{inc}}.

      Example:

      .. code-block:: text

        subsection bc 0
          set id   = 0
          set type = electric field
          subsection E x real part
            set Function expression = 1
          end
          subsection E x imag part
            set Function expression = 1
          end
          subsection E y real part
            set Function expression = 0
          end
          subsection E y imag part
            set Function expression = 0
          end
          subsection E z real part
            set Function expression = 0
          end
          subsection E z imag part
            set Function expression = 0
          end
        end

    * ``magnetic field`` to impose a prescribed magnetic field. The boundary values are taken from the parsed functions in the ``H x/y/z real part`` and ``H x/y/z imag part`` subsections,

      .. math::
        \mathbf{H} = \mathbf{H}^{\mathrm{inc}} - \mathbf{J}_\mathrm{s}.

      Example:

      .. code-block:: text

        subsection bc 0
          set id   = 0
          set type = magnetic field
          subsection H x real part
            set Function expression = 1
          end
          subsection H x imag part
            set Function expression = 1
          end
          subsection H y real part
            set Function expression = 0
          end
          subsection H y imag part
            set Function expression = 0
          end
          subsection H z real part
            set Function expression = 0
          end
          subsection H z imag part
            set Function expression = 0
          end
        end

    * ``impedance boundary`` to impose an impedance boundary condition. The boundary is governed by the surface admittance :math:`Y_s` supplied through the ``surface admittance real part`` and ``surface admittance imag part`` subsections and the electromagnetic excitation :math:`\mathbf{g}` supplied through the ``excitation x/y/z real part`` and ``excitation x/y/z imag part`` subsections,

      .. math::
        \mathbf{n} \times \mathbf{H} + Y_s \mathbf{n} \times ( \mathbf{E} \times \mathbf{n} ) = \mathbf{g}.

      Example:

      .. code-block:: text

        subsection bc 0
          set id   = 0
          set type = impedance boundary
          subsection surface admittance real part
            set Function expression = 0.1
          end
          subsection surface admittance imag part
            set Function expression = 0
          end
          subsection excitation x real part
            set Function expression = 0
          end
          subsection excitation x imag part
            set Function expression = 0
          end
          subsection excitation y real part
            set Function expression = 0
          end
          subsection excitation y imag part
            set Function expression = 0
          end
          subsection excitation z real part
            set Function expression = 0
          end
          subsection excitation z imag part
            set Function expression = 0
          end
        end

    * ``waveguide port`` to impose a waveguide port boundary condition. The parameters for this boundary conditions are defined in the Time-Harmonic Maxwell section of the input file. In the following, :math:`k_l` is the wave number in the longitudinal direction of the waveguide.

      .. math::
        \mathbf{n} \times \mathbf{H} + \frac{k_{l}}{\omega \mu_r} \mathbf{n} \times ( \mathbf{E} \times \mathbf{n} ) = \mathbf{n} \times \mathbf{H}_{\mathrm{TE}_{mn}} + \frac{k_{l}}{\omega \mu_r} \mathbf{n} \times ( \mathbf{E}_{\mathrm{TE}_{mn}} \times \mathbf{n} )

      Example:

      .. code-block:: text

        subsection bc 0
          set id   = 0
          set type = waveguide port
        end


