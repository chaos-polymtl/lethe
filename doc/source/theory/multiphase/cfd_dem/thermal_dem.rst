
--------------------------------
Thermal DEM in a Stagnant Gas
--------------------------------

The heat transfer mechanisms considered here are:
   - Conduction through the particles themselves
   -  Conduction between particles through the contact surface (microcontacts and macrocontacts)
   - Conduction through the interstitial fluid (microgap and macrogap between particles)

* The temperature is uniform within each particle.
* The temperature of each particle changes slowly enough that thermal disturbances do not propagate beyond immediate neighbors.
* Convection and radiation are neglected.
* The interstitial fluid is a stagnant gas (usually air), through which conduction occurs for contacting particles.


Using Batchelor and O'Brien's model [#Batchelor1977]_ for thermal DEM:

.. math::

   \frac{d T_i}{dt} = \frac{Q_i}{m_i c_i} 

with

.. math::

   Q_{i j} = H_{i j} (T_j - T_i)

where :math:`Q_i \approx \sum Q_{ij} + Q_{Si}` is the heat transfer rate, and :math:`H_{ij}` is the thermal conductance.

The thermal conductance, which is the inverse of the thermal resistance :math:`R_{ij}`, is calculated as follows by Beaulieu *et al*. [#Beaulieu2020]_:

.. math::

   \frac{1}{R_{ij}} = \frac{1}{R_L + \left( \frac{1}{R_s} + \frac{1}{R_g} \right)^{-1}} + \frac{1}{R_c + R_G}

Where:

* :math:`R_L` resistance of the contact surface [#Batchelor1977]_
* :math:`R_s` resistance of the microcontacts [#VanLew2016]_
* :math:`R_g` resistance of the interstitial gas microgap [#Bahrami2006]_
* :math:`R_c` resistance of the solid layers of the particles [#Beaulieu2020]_
* :math:`R_G` resistance of the interstitial gas macrogap [#Bahrami2006]_

.. figure:: Schematic_resistances.png
   :width: 100%
   :alt: Schematic of thermal resistances
   :align: center

   Schematic of the heat transfer between two particles in contact [#Beaulieu2020]_.

~~~~~~~~~~~~~~~~~~~~
Thermal Resistances
~~~~~~~~~~~~~~~~~~~~

The thermal resistances which model the heat transfer between particles are calculated as follows:

.. math::

   R_L &= \frac{1}{2 k_h r_c } \\
   R_s &= \left(\frac{H'}{P_0}\right)^{0.96} \frac{1.184}{\pi r_c^2 k_h}\left(\frac{\sigma}{\tau}\right) \\
   R_c &= R_{c,i} +R_{c,j}, \quad R_{c,i} = \frac{L_i}{k_i A_i} \\
   R_g &= \frac{2\sqrt{2}\sigma a_2}{\pi k_g r_c^2 \ln\left(1+\frac{a_2}{a_1+M/(2\sqrt{2}\sigma)}\right)} \\
   R_G &= \frac{2}{\pi k_g \left[S \ln\left(\frac{S-B}{S-A}\right) + B - A\right]} \\


The contact radius :math:`r_c` is calculated as follows:

.. math::

   r_c = \left( \frac{3F_n r^* }{4E^*}\right)^{1/3}

In [#Zhou2010]_, a factor c is introduced to correct the contact radius, which can be overestimated when the Young's modulus in the simulation is underestimated for computational efficiency.

.. math::

   r_c' = r_c c, \quad c = \left( \frac{E^*_{Sim}}{E^*_{Real}} \right)^{1/5}

.. note::
   For now, the parameter for the real young modulus of the particles is not implemented so the factor c is equal to 1.
   

Where:

* :math:`k_h = \frac{2k_ik_j}{k_i+k_j}` harmonic mean of the particles' thermal conductivities
* :math:`H'` harmonic mean of the particles' microhardnesses
* :math:`E^* = \left( \frac{(1-\nu_i^2)}{E_i} + \frac{(1-\nu_j^2)}{E_j}\right)^{-1}` effective Young's modulus
* :math:`r^* = \frac{r_ir_j}{r_i+r_j}` effective radius
* :math:`P_0 = \frac{2E^*\delta_n}{\pi r_c}` maximum pressure for hertzian contacts
* :math:`\sigma = \sqrt{\sigma_i^2 + \sigma_j^2}` equivalent surface roughness
* :math:`\tau = \sqrt{\tau_i^2 + \tau_j^2}` equivalent surface slope
* :math:`L_i = \frac{\pi r_i}{4}` characteristic length parallel to the heat flux
* :math:`A_i = \pi(r_i^2 - r_c^2)` characteristic area perpendicular to the heat flux
* :math:`a_1 = erfc^{-1}(2P_0/H'), \quad a_2 = erfc^{-1}(0.03P_0/H') - a_1`
* :math:`A = 2\sqrt{r_h^2 - r_c^2}, \quad B = 0 \quad` for simple cubic packing, :math:`S = 2\left(r_h - \frac{r_c^2}{2r_h}\right) + M.`
* :math:`M = \left( \frac{2-\alpha_{T_i}}{\alpha_{T_i}} + \frac{2-\alpha_{T_j}}{\alpha_{T_j}}\right)\left(\frac{2\gamma_g}{1+\gamma_g}\right)\frac{\Lambda}{Pr}` gas parameter
   with:

   * :math:`\alpha_{T_i}` thermal accomodation coefficient particle-gas
   * :math:`\gamma_g` specific heats ratio of the gas
   * :math:`Pr = \frac{\mu_g c_g}{k_g}` Prandlt number of the gas
   * :math:`\Lambda` molecular mean free path of the gas

   

-----------
References
-----------

.. [#Batchelor1977] \G. K. Batchelor and R. W. O’Brien, “Thermal or electrical conduction through a granular material,” Proc. R. Soc. Lond. A Math. Phys. Sci., vol. 355, no. 1682, pp. 313–333, Jul. 1977, doi: 10.1098/rspa.1977.0100 <https://doi.org/10.1098/rspa.1977.0100>_.

.. [#Beaulieu2020] \C. Beaulieu, Impact de la ségrégation granulaire sur le transfert de chaleur dans un lit rotatif, Ph.D. thesis, Polytechnique Montréal, 2020. Available: https://publications.polymtl.ca/4757/

.. [#VanLew2016] \J. T. Van Lew, On thermal characterization of breeder pebble beds with microscale numerical modeling of thermofluid and pebble-pebble interactions, Ph.D. dissertation, University of California, Los Angeles, 2016.

.. [#Bahrami2006] \M. Bahrami, M. M. Yovanovich, and J. R. Culham, “Effective thermal conductivity of rough spherical packed beds,” Int. J. Heat Mass Transf., vol. 49, no. 19–20, pp. 3691–3701, Sep. 2006, doi: 10.1016/j.ijheatmasstransfer.2006.02.021 <https://doi.org/10.1016/j.ijheatmasstransfer.2006.02.021>_.

.. [#Zhou2010] \Z. Y. Zhou, A. B. Yu, and P. Zulli, “A new computational method for studying heat transfer in fluid bed reactors,” Powder Technol., vol. 197, no. 1–2, pp. 102–110, Sep. 2010, doi: 10.1016/j.powtec.2009.09.002 <https://doi.org/10.1016/j.powtec.2009.09.002>_.