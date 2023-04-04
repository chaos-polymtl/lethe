***********************************************
CFD-DEM
***********************************************
This subsection includes parameters related to multiphase flow simulations using the both the gls_vans solver and the cfd-dem_coupling solver within Lethe.

.. code-block:: text

   subsection cfd-dem
      set grad div = true
      set void fraction time derivative = true
      set interpolated void fraction = true
      set vans model = modelA
      set drag force = true
      set drag model = difelice
      set saffman lift force = false
      set magnus lift force = false
      set buoyancy force = true
      set shear force = true
      set pressure force = true
      set coupling frequency = 100
      set implicit stabilization = true
      set grad-div length scale = 1
      set particle statistics = true
   end


* The ``grad div`` parameter allows the enabling of the grad div stabilization for the Volume Averaged Navier Stokes equations `[1] <https://doi.org/10.1016/j.softx.2020.100579>`_. This allows for a much better mass conservation of the system.
* The ``void fraction time derivative`` parameter allows us to choose whether or not we want to account for the time derivative of the void fraction or take it equal to zero.
* The ``interpolated void fraction`` parameter allows us choose whether the void fraction used to calculate drag is the cell void fraction or the one interpolated at the position of the particle (Using the cell void fraction to calculate drag on each particle instead of the interpolated one is currently under investigation).
* The ``vans model`` parameter allows us to choose between the vans Model A or Model B. Details about the differences between the models are provided in Lethe's unresolved CFD-DEM theory guide :doc:`../../theory/unresolved_cfd-dem/unresolved_cfd-dem`.
* The ``drag force``, ``saffman lift force``, ``magnus lift force``, ``buoyancy force``, ``shear force``, and ``pressure force`` parameters allow us to enable or disable the respective forces in a cfd-dem simulation.

.. note::
    By setting ``set saffman lift force = true``, the applied Saffman lift force model is the often called Saffman-Mei model, developed by Mei (1992) `[2] <https://doi.org/10.1016/0301-9322(92)90012-6>`_ as an extension of the work by Saffman (1968) `[3] <https://doi.org/10.1017/S0022112065000824>`_. A complete description of the model is provided by Crowe et al. (2010) `[4] <https://doi.org/10.1201/b11103>`_.

.. note::
    By setting ``set magnus lift force = true``, the applied Magnus lift force model is detailed by Crowe et al. (2010) `[4] <https://doi.org/10.1201/b11103>`_, following the recommendation of using the correlation by Oesterlé & Dinh (1998) `[5] <https://doi.org/10.1007/s003480050203>`_ for :math:`1 < \Omega < 6` and :math:`10 < Re_p < 140`. :math:`\Omega` is calculated as:

    .. math::
        \Omega = \frac{d_p \omega}{2 \left | u - v \right |}

    where :math:`\omega` is the angular velocity of the particle.

 .. warning:: 
   We do not recommend using the Magnus lift force. The current model does include any angular momentum dissipation mechanism in the solid-fluid coupling. Using the Magnus force may lead to unphysical results.

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

`[1] <https://doi.org/10.1016/j.softx.2020.100579>`_ B. Blais, L. Barbeau, V. Bibeau, S. Gauvin, T. E. Geitani, S. Golshan, R. Kamble, G. Mirakhori, J. Chaouki, Lethe: An open-source parallel high- order adaptative cfd solver for incompressible flows, SoftwareX 12 100579, 2020.

`[2] <https://doi.org/10.1016/0301-9322(92)90012-6>`_ R. Mei, An approximate expression for the shear lift force on a spherical particle at finite Reynolds number. International Journal of Multiphase Flow, v. 18, n. 1, p. 145-147, 1992.

`[3] <https://doi.org/10.1017/S0022112065000824>`_ P. G. Saffman, The lift on a small sphere in a slow shear flow. Journal of fluid mechanics, v. 22, n. 2, p. 385-400, 1965.

`[4] <https://doi.org/10.1201/b11103>`_ C.T. Crowe, J.D. Schwarzkopf, M. Sommerfeld, Y. Tsuji, . Multiphase Flows with Droplets and Particles (2nd ed.). CRC Pres, 2011.

`[5] <https://doi.org/10.1007/s003480050203>`_ B. Oesterlé, T. Dinh, Experiments on the lift of a spinning sphere in a range of intermediate Reynolds numbers. Experiments in Fluids 25, 16–22, 1998.

`[6] <https://doi.org/10.1016/0301-9322(94)90011-6>`_ R. Di Felice, The voidage function for fluid-particle interaction systems. International journal of multiphase flow 20 (1), 153–159, 1994.

`[7] <https://doi.org/10.1016/j.ces.2013.05.036>`_ L. Rong, K. Dong, A. Yu, Lattice-boltzmann simulation of fluid flow through packed beds of uniform spheres: Effect of porosity, Chemical engineering science 99, 44–58, 2013.

`[8] <https://doi.org/10.1080/07373937.2010.482714>`_ W. Sobieski. Drag Coefficient in Solid–Fluid System Modeling with the Eulerian Multiphase Model. Drying Technology, 29, 111-125, 2011.

`[9] <https://doi.org/10.1016/j.ces.2013.05.014>`_  D. Jajcevic, E. Siegmann, C. Radeke, J. G. Khinast, Large-scale cfd–dem simulations of fluidized granular systems. Chemical Engineering Science 98, 298–310, 2013.

`[10] <https://doi.org/10.1016/j.ijmultiphaseflow.2020.103425>`_ Tim M.J. Nijssen, Hans A.M. Kuipers, Jan van der Stel, Allert T. Adema, Kay A. Buist, Complete liquid-solid momentum coupling for unresolved CFD-DEM simulations, International Journal of Multiphase Flow, Volume 132, 2020.

`[11] <https://doi.org/10.1016/j.powtec.2019.10.058>`_ F. Marchelli, Q. Hou, B.Bosio, E. Arato, & A. Yu, Comparison of different drag models in CFD-DEM simulations of spouted beds. Powder Technology, 360, 1253-1270, 2020.
