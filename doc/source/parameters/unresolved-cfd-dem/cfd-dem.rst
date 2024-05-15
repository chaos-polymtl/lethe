=======
CFD-DEM
=======
This subsection includes parameters related to multiphase flow simulations using the both the gls_vans solver and the cfd-dem_coupling solver within Lethe.

.. code-block:: text

  subsection cfd-dem
    set grad div                      = true
    set void fraction time derivative = true
    set interpolated void fraction    = true
    set vans model                    = modelA
    set drag force                    = true
    set drag model                    = difelice
    set saffman lift force            = false
    set magnus lift force             = false
    set rotational viscous torque     = false
    set vortical viscous torque       = false
    set buoyancy force                = true
    set shear force                   = true
    set pressure force                = true
    set coupling frequency            = 100
    set implicit stabilization        = true
    set grad-div length scale         = 1
    set particle statistics           = true
  end


* The ``grad div`` parameter allows the enabling of the grad div stabilization for the Volume Averaged Navier Stokes equations `[1] <https://doi.org/10.1016/j.softx.2020.100579>`_. This allows for a much better mass conservation of the system.
* The ``void fraction time derivative`` parameter allows us to choose whether or not we want to account for the time derivative of the void fraction or take it equal to zero.
* The ``interpolated void fraction`` parameter allows us choose whether the void fraction used to calculate drag is the cell void fraction or the one interpolated at the position of the particle (Using the cell void fraction to calculate drag on each particle instead of the interpolated one is currently under investigation).
* The ``vans model`` parameter allows us to choose between the vans Model A or Model B. Details about the differences between the models are provided in Lethe's unresolved CFD-DEM theory guide :doc:`../../theory/multiphase/cfd_dem/unresolved_cfd-dem`.
* The ``drag force``, ``saffman lift force``, ``magnus lift force``, ``buoyancy force``, ``shear force``, and ``pressure force`` parameters allow us to enable or disable the respective forces in a cfd-dem simulation.

.. note::
    By setting ``set saffman lift force = true``, the applied Saffman lift force model is the often called Saffman-Mei model, developed by Mei (1992) `[2] <https://doi.org/10.1016/0301-9322(92)90012-6>`_ as an extension of the work by Saffman (1968) `[3] <https://doi.org/10.1017/S0022112065000824>`_. A complete description of the model is provided by Crowe et al. (2010) `[4] <https://doi.org/10.1201/b11103>`_.

.. note::
    By setting ``set magnus lift force = true``, the applied Magnus lift force model is detailed by Crowe et al. (2010) `[4] <https://doi.org/10.1201/b11103>`_, following the recommendation of using the correlation by Oesterlé & Dinh (1998) `[5] <https://doi.org/10.1007/s003480050203>`_ for :math:`1 < \Omega < 6` and :math:`10 < Re_p < 140`. :math:`\Omega` is calculated as:

    .. math::
        \Omega = \frac{d_p \omega}{2 \left | u - v \right |}

    where :math:`\omega` is the angular velocity of the particle.

 .. warning:: 
   We do not recommend using the Magnus lift force. The Magnus lift force model does not include any angular momentum dissipation mechanism in the solid-fluid coupling. Using the Magnus force may lead to unphysical results.

* The ``rotational viscous torque`` and ``vortical viscous torque`` parameter controls whether the fluid-particle contact generates torque on the particles due to viscosity.

.. note::

    When ``set rotational viscous torque = true`` and ``set vortical viscous torque = true``, the applied torque (:math:`\bf{M}_{viscous}`) is the one described by Derksen (2004) `[6] <https://doi.org/10.1002/aic.690491104>`_:

    .. math::
        \bf{M}_{viscous} = \pi d_p^3 \mu \left ( 0.5 \bf{\omega}_f - \bf{\omega}_p \right )

    where :math:`\bf{\omega}_f` is the fluid vorticity at particle's position and :math:`\bf{\omega}_p` is the particle's angular velocity. The rotational and vortical torques can be applied separately by setting one of them to `false`.

    In case ``set rotational viscous torque = false``, the particle's angular velocity :math:`\bf{\omega}_p` is removed from the equation.
    In case ``set vortical viscous torque = false``, :math:`0.5 \bf{\omega}_f` is removed from the equation.

.. warning::
    We do not recommend the use of ``vortical viscous torque`` with coarse meshes, especially when Q1 elements are used. In such case, the space resolution may not be enough to properly capture vorticity.
    Since the viscous torque model is not complete without the vortical component, ``rotational viscous torque`` should be used with caution.

* The ``drag model`` parameter allows one to choose the type of drag model to be implemented for the calculation of the drag force between the particles and the fluids. Given :math:`F_d = \beta (\bf{u} - \bf{v})`, the available drag models at the time are:

.. csv-table::
   :file: tables/drag_models_unresolved_cfd-dem.csv
   :header-rows: 1
   :align: center

* The ``particle statistics`` parameter, when enabled, outputs statistics about the particles' velocity, kinetic energy, and the amount of contact detection.
* The ``coupling frequency`` determines the number of DEM iterations per 1 CFD iteration.

.. note::
   The ``coupling frequency`` parameter is used to calculate the dem time step as it is not explicitly determined in the parameter file. It is calculated as: 

   .. math::
      \Delta t_{DEM} = \frac{\Delta t_{CFD}}{coupling frequency}

* The ``implicit stabilization`` parameter determines whether or not we calculate the :math:`\tau` for the SUPG/PSPG stabilization and the :math:`\gamma` for the grad-div stabilization using the current velocity (implicit stabilization) or the velocity at the previous time step (explicit stabilization). By default, this is set to true. If difficulties are encountered in the convergence of the non-linear solver, a good practice is to set this to false.
* The ``grad-div length scale`` parameter determines the value of the length scale constant :math:`c^*` in the calculation of :math:`\gamma = \nu + c^* \mathbf{u}`.

.. tip::
   Experience shows that simulations are more numerically stable when the ``grad-div length scale`` is of the same length as the characteristic length of the flow. For example, for a pipe, the recommended value for the ``grad-div length scale`` would be the pipe's diameter.

`[1] <https://doi.org/10.1016/j.softx.2020.100579>`_ B. Blais *et al.*, “Lethe: An open-source parallel high-order adaptative CFD solver for incompressible flows,” *SoftwareX*, vol. 12, p. 100579, Jul. 2020, doi: 10.1016/j.softx.2020.100579.

`[2] <https://doi.org/10.1016/0301-9322(92)90012-6>`_ R. Mei, “An approximate expression for the shear lift force on a spherical particle at finite reynolds number,” Int. J. *Multiph. Flow*, vol. 18, no. 1, pp. 145–147, Jan. 1992, doi: 10.1016/0301-9322(92)90012-6.

`[3] <https://doi.org/10.1017/S0022112065000824>`_ P. G. Saffman, “The lift on a small sphere in a slow shear flow,” *J. Fluid Mech.*, vol. 22, no. 2, pp. 385–400, Jun. 1965, doi: 10.1017/S0022112065000824.

`[4] <https://doi.org/10.1201/b11103>`_ 	C. T. C. Tsuji John D. Schwarzkopf, Martin Sommerfeld, Yutaka, *Multiphase Flows with Droplets and Particles*, 2nd ed. Boca Raton: CRC Press, 2011. doi: 10.1201/b11103.

`[5] <https://doi.org/10.1007/s003480050203>`_ B. Oesterlé and T. B. Dinh, “Experiments on the lift of a spinning sphere in a range of intermediate Reynolds numbers,” *Exp. Fluids*, vol. 25, no. 1, pp. 16–22, Jun. 1998, doi: 10.1007/s003480050203.

`[6] <https://doi.org/10.1002/aic.690491104>`_ J. J. Derksen, “Numerical simulation of solids suspension in a stirred tank,” *AIChE J.*, vol. 49, no. 11, pp. 2700–2714, 2003, doi: 10.1002/aic.690491104.

`[7] <https://doi.org/10.1016/0301-9322(94)90011-6>`_ R. Di Felice, “The voidage function for fluid-particle interaction systems,” *Int. J. Multiph. Flow*, vol. 20, no. 1, pp. 153–159, Feb. 1994, doi: 10.1016/0301-9322(94)90011-6.

`[8] <https://doi.org/10.1016/j.ces.2013.05.036>`_ L. W. Rong, K. J. Dong, and A. B. Yu, “Lattice-Boltzmann simulation of fluid flow through packed beds of uniform spheres: Effect of porosity,” *Chem. Eng. Sci.*, vol. 99, pp. 44–58, Aug. 2013, doi: 10.1016/j.ces.2013.05.036.

`[9] <https://doi.org/10.1080/07373937.2010.482714>`_ W. Sobieski, “Drag Coefficient in Solid–Fluid System Modeling with the Eulerian Multiphase Model,” *Dry. Technol.*, vol. 29, no. 1, pp. 111–125, Dec. 2010, doi: 10.1080/07373937.2010.482714.

`[10] <https://doi.org/10.1016/j.ces.2013.05.014>`_  D. Jajcevic, E. Siegmann, C. Radeke, and J. G. Khinast, “Large-scale CFD–DEM simulations of fluidized granular systems,” *Chem. Eng. Sci.*, vol. 98, pp. 298–310, Jul. 2013, doi: 10.1016/j.ces.2013.05.014.

`[11] <https://doi.org/10.1016/j.ijmultiphaseflow.2020.103425>`_ T. M. J. Nijssen, H. A. M. Kuipers, J. van der Stel, A. T. Adema, and K. A. Buist, “Complete liquid-solid momentum coupling for unresolved CFD-DEM simulations,” *Int. J. Multiph. Flow*, vol. 132, p. 103425, Nov. 2020, doi: 10.1016/j.ijmultiphaseflow.2020.103425.

`[12] <https://doi.org/10.1016/j.powtec.2019.10.058>`_ F. Marchelli, Q. Hou, B. Bosio, E. Arato, and A. Yu, “Comparison of different drag models in CFD-DEM simulations of spouted beds,” *Powder Technol.*, vol. 360, pp. 1253–1270, Jan. 2020, doi: 10.1016/j.powtec.2019.10.058.
