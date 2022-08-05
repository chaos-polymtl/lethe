***********************************************
CFD-DEM
***********************************************
This subsection includes parameters related to multiphase flow simulations using the both the gls_vans solver and the cfd-dem_coupling solver within Lethe.

.. code-block:: text

   subsection cfd-dem
      set grad div = false
      set void fraction time derivative = false
      set interpolated void fraction = true
      set drag force = true
      set saffman lift force = false
      set buoyancy force = true
      set shear force = true
      set pressure force = true
      set drag model = dallavalle
      set post processing = true
      set coupling frequency = 200
      set inlet boundary id = 2
      set outlet boundary id = 3
      set implicit stabilization = true
      set grad-div length scale = 1
   end


* The ``grad div`` parameter allows the enabling of the grad div stabilization for the Volume Averaged Navier Stokes equations `[1] <https://doi.org/10.1016/j.softx.2020.100579>`_. This allows for a much better mass conservation of the system.
* The ``void fraction time derivative`` parameter allows us to choose whether or not we want to account for the time derivative of the void fraction or take it equal to zero.
* The ``interpolated void fraction`` parameter allows us choose whether the void fraction used to calculate drag is the cell void fraction or the one interpolated at the position of the particle (Using the cell void fraction to calculate drag on each particle instead of the interpolated one is currently under investigation).
* The ``drag force``, ``saffman lift force``, ``buoyancy force``, ``shear force``, and ``pressure force`` parameters allow us to enable or disable the respective forces in a cfd-dem simulation.

.. note::
    By setting ``set saffman lift force = true``, the applied Saffman lift force model is the often called Saffman-Mei model, developed by Mei (1992) `[2] <https://doi.org/10.1016/0301-9322(92)90012-6>`_ as an extension of the work by Saffman (1968) `[3] <https://doi.org/10.1017/S0022112065000824>`_. A complete description of the model is provided by Crowe et al. (2010) `[4] <https://doi.org/10.1201/b11103>`_.

* The ``drag model`` parameter allows one to choose the type of drag model to be implemented for the calculation of the drag force between the particles and the fluids. Available drag models at the time of writing are: Di Felice `[5] <https://doi.org/10.1016/0301-9322(94)90011-6>`_, Rong `[6] <https://doi.org/10.1016/j.ces.2013.05.036>`_, Dallavalle `[7] <https://doi.org/10.1080/07373937.2010.482714>`_, Koch and Hill `[8] <https://doi.org/10.1016/j.ces.2013.05.014>`_, Beetstra `[9] <https://doi.org/10.1002/aic.11065>`_, and Gidaspow `[10] <https://books.google.ca/books?id=fHecceQyaYkC&lpg=PP1&ots=uhExYvWrkv&lr&hl=pt-BR&pg=PP1#v=onepage&q&f=false>`_.
* The ``post processing`` parameter, when enabled, allows the calculation of the pressure drop, void fraction in the packed region, and the mass conservation in a packed bed at each time step.
* The ``coupling frequency`` parameter is only applicable for the cfd-dem solver and it determines the number of DEM iterations per 1 CFD iteration.

.. note::
   The ``coupling frequency`` parameter is used to calculate the dem time step as it is not explicitly determined in the parameter file. It is calculated as: 

   .. math::
      \Delta t_{DEM} = \frac{\Delta t_{CFD}}{coupling frequency}

* The ``inlet boundary`` and ``outlet boundary`` parameters allow us to specify the ID of the inlet and outlet for pressure drop calculations.
* The implicit stabilization parameter determines wehether or not we calculate the :math:`\tau` for the SUPG stabilization and the :math:`\gamma` for the grad-div stabilization using the current velocity (implicit stabilization) or the velocity at the previous time step (explicit stabilization). By default, this is set to true. If difficulties are encountered in the convergence of the non-linear solver, a good practice is to set this to false.
* The grad-div length scale determines the value of the length scale constant :math:`c^*` in the calculation of :math:`\gamma = \nu + c^* \mathbf{u}`. 

`[1] <https://doi.org/10.1016/j.softx.2020.100579>`_ B. Blais, L. Barbeau, V. Bibeau, S. Gauvin, T. E. Geitani, S. Golshan, R. Kamble, G. Mirakhori, J. Chaouki, Lethe: An open-source parallel high- order adaptative cfd solver for incompressible flows, SoftwareX 12 (2020) 100579.

`[2] <https://doi.org/10.1016/0301-9322(92)90012-6>`_ Mei, Renwei. An approximate expression for the shear lift force on a spherical particle at finite Reynolds number. International Journal of Multiphase Flow, v. 18, n. 1, p. 145-147, 1992.

`[3] <https://doi.org/10.1017/S0022112065000824>`_ SAFFMAN, Philip Geoffrey. The lift on a small sphere in a slow shear flow. Journal of fluid mechanics, v. 22, n. 2, p. 385-400, 1965.

`[4] <https://doi.org/10.1201/b11103>`_ Crowe, C.T., Schwarzkopf, J.D., Sommerfeld, M., & Tsuji, Y. (2011). Multiphase Flows with Droplets and Particles (2nd ed.). CRC Pres.

`[5] <https://doi.org/10.1016/0301-9322(94)90011-6>`_ R. Di Felice, The voidage function for fluid-particle interaction systems, International journal of multiphase flow 20 (1) (1994) 153–159.

`[6] <https://doi.org/10.1016/j.ces.2013.05.036>`_ L. Rong, K. Dong, A. Yu, Lattice-boltzmann simulation of fluid flow through packed beds of uniform spheres: Effect of porosity, Chemical engineering science 99 (2013) 44–58.

`[7] <https://doi.org/10.1080/07373937.2010.482714>`_ Sobieski, Wojciech. (2011). Drag Coefficient in Solid–Fluid System Modeling with the Eulerian Multiphase Model. Drying Technology. 29. 111-125. 10.1080/07373937.2010.482714.

`[8] <https://doi.org/10.1016/j.ces.2013.05.014>`_  D. Jajcevic, E. Siegmann, C. Radeke, J. G. Khinast, Large-scale cfd–dem simulations of fluidized granular systems, Chemical Engineering Science 98 (2013) 298–310.

`[9] <https://doi.org/10.1002/aic.11065>`_ R. Beetstra, M. A. van der Hoef, J. A. M. Kuipers. Drag Force of Intermediate Reynolds Number Flow Past Mono- and Bidisperse Arrays of Spheres. AIChE journal, v. 53, n. 2, p. 489-501, 2007.

`[10] <https://books.google.ca/books?id=fHecceQyaYkC&lpg=PP1&ots=uhExYvWrkv&lr&hl=pt-BR&pg=PP1#v=onepage&q&f=false>`_ D. Gidaspow. Multiphase flow and fluidization: continuum and kinetic theory descriptions. Academic press, 1994.