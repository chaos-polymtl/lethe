==================================================
Photon Count Calculation in a Cylindrical Vessel
==================================================

In this example, using a Monte Carlo technique, we perform the calculation of photon counts of a single radioactive particle that emits :math:`\gamma`-rays. The calculation is performed for a given set of positions inside a cylindrical vessel. The Monte Carlo method allows us to estimate the photon counts of a particle at a given position inside the vessel with respect to a given detector.

--------
Features
--------

- Solver: ``lethe-rpt-3d``
- Displays the use of the Monte Carlo method in the calculation of photon count


---------------------------
Files Used in This Example
---------------------------

All files mentioned below are located in the example's folder (``examples/rpt/count-calculation``).

- File containing detector positions: ``positions.detector``
- File containing particle positions for the first scenario:  ``positions-horizontalx.particle``
- File containing particle positions for the second scenario  ``positions-horizontaly.particle``
- File containing particle positions for the third scenario:  ``positions-vertical.particle``
- File containing particle positions for the fourth scenario:  ``positions-diagonal.particle``
- Parameter file: ``rpt-count-calculation.prm``
- Postprocessing Python script: ``rpt_count-calculation_plot.py``


-------------------------
Description of the Case
-------------------------
In this example, four different sets of particle positions are studied for a given detector position. The four different scenarios studied in this example are:

1. Horizontal translation of a particle along the x-axis
2. Horizontal translation of a particle along the y-axis
3. Vertical translation of a particle along the z-axis
4. Particle going across the vessel on a diagonal line


The illustration below depicts the geometry of the vessel, the detector, and the path traveled by the particle for each scenario:

.. image:: images/scenarios.png
    :alt: Scenarios
    :align: center
    :name: geometry_description

As a particle travels in the cylindrical vessel, its photon count (:math:`C`) measured by the detector varies according to the following relation:

.. math::
    C = \frac{T \nu R \phi \xi_i (\mathbf{X})}{1 + \tau \nu R \phi \xi_i (\mathbf{X})}
		
where

- :math:`T` is the sampling time (:math:`s`);
- :math:`\nu` is the number of :math:`\gamma`-rays emitted by each disintegration;
- :math:`R` is the activity of the tracer (:math:`Beq`);
- :math:`\phi` is the peak-to-total ratio;
- :math:`\tau` is the dead time of the detector (:math:`s`);
- :math:`\mathbf{X}` is the tracer particle's position, and
- :math:`\xi_i(\mathbf{X})` is the efficiency of the :math:`i_{th}` detector related to the position :math:`\mathbf{X}`.


The efficiency of the detector may be expressed by means of the following equation:

.. math::
	

    \xi_i (\mathbf{X}) = \int_{\Omega } \frac{\mathbf{r}\cdot d\overrightarrow{A}}{\left \| \mathbf{r} \right \|^{3}}f_{a}(\alpha ,\theta )f_{d}(\alpha ,\theta )
	

where

- :math:`\Omega` is the closed exposed surface of the detector;
- :math:`\mathbf{r}` is a vector going from the position of the tracer particle (:math:`\mathbf{X}`) to a variable point (:math:`\mathbf{P}`) on the exposed surface of the detector;
- :math:`d\mathbf{A}` is the external surface vector normal to the surface at the contact point on the detector's crystal;
- :math:`f_a(\alpha, \theta)` is the probability function of the non-interaction between the :math:`\gamma`-rays emitted within :math:`\Omega` and the material inside the vessel, and
- :math:`f_d(\alpha, \theta)` is the probability function of the interaction of the :math:`\gamma`-rays with the detector. 

The two last functions may be re-written the following way:

.. math::

    f_a(\alpha, \theta) = exp\{-\mu_r \ e(\alpha, \theta)\}

where :math:`\mu_r` is the reactor's attenuation coefficient and :math:`e(\alpha, \theta)` is the length of the path traveled by the photon inside the vessel.


And

.. math::

    f_d(\alpha, \theta) = 1 - exp\{ -\mu_d \ d(\alpha,\theta)\}

where :math:`\mu_d` is the detector's attenuation coefficient and :math:`d(\alpha,\theta)` is the length of the path traveled by the photon inside the detector.



Using the Monte Carlo algorithm, we approximate the previous closed surface integral by randomly selecting several thousands of photon path directions.

Thus, the efficiency of the :math:`i_{th}` detector is calculated as follows:

.. math::

    \xi_i (\mathbf{X}) = \frac{1}{N} \sum_{j=1}^{N} \omega(\alpha) \omega(\theta) f_a(\alpha_j, \theta_j) f_d(\alpha_j, \theta_j)


where

- :math:`N` is the number of randomly generated photons;
- :math:`\alpha_j` and :math:`\theta_j` are randomly generated angles that describe the direction of a ray emitted by a tracer particle;
- :math:`\omega(\alpha)` is the weighting factor associated with the angle :math:`\alpha`, and
- :math:`\omega(\theta)` is the weighting factor associated with the angle :math:`\theta`.


----------------
Parameter File
----------------

RPT Parameters
~~~~~~~~~~~~~~~

In the subsection ``rpt parameters``, we define the values of the set of parameters necessary for calculating the counts using the Monte Carlo method.  Among these parameters, we have the name of the file which contains a set of different positions of the particle inside the vessel (``particle position file``), the number of Monte Carlo iterations (``monte carlo iteration``), the seed that is used to generate a random number (``random number seed``) and other parameters that describe the studied :math:`\gamma`-ray model. We also define the name of the file in which the counts for each position will be exported in with the parameter ``counts file``. These common parameters used for the RPT simulation are described in the :doc:`../../../parameters/rpt/rpt_parameters` documentation page.

.. code-block:: text

    subsection rpt parameters
      set particle positions file          = positions-horizontalx.particle
      set verbosity                        = verbose
      set export counts                    = true
      set counts file                      = counts_horizontalx.csv
      set monte carlo iteration            = 100000
      set random number seed               = 0
      set reactor height                   = 0.3
      set reactor radius                   = 0.1
      set peak-to-total ratio              = 0.4
      set sampling time                    = 1
      set gamma-rays emitted               = 2
      set attenuation coefficient detector = 21.477
    end


Detector Parameters
~~~~~~~~~~~~~~~~~~~~

In the subsection ``detector parameters``, we specify the file that contains two positions located on the axis of symmetry of the detector. The first point is on the surface facing the vessel (face of the detector), and the second point can be any point located inside the detector. In the current example, the center position of the face is :math:`(0.200, 0, 0.075)`, and the second point on the axis is :math:`(0.238, 0, 0.075)`. We also specify the radius (``radius``) and the length (``length``) of the detector. A detailed description of these parameters can be found in the :doc:`../../../parameters/rpt/detector_parameters` documentation page.


.. code-block:: text

    subsection detector parameters
      set detector positions file         = positions.detector
      set radius                          = 0.0381
      set length                          = 0.0762
      set dead time                       = 1e-5
      set activity                        = 2e6
      set attenuation coefficient reactor = 10
    end

.. note::
    The parameters ``dead time``, ``activity`` and ``attenuation coefficient reactor`` are obtained using the blackbox optimization software `NOMAD <https://www.gerad.ca/en/software/nomad/>`_ . The second example :doc:`../tuning-parameters-with-nomad/tuning-parameters-with-nomad` explains how we can obtain the values of these parameters using NOMAD.


----------------------------------
Running the Simulation
----------------------------------

Launching the simulation is as simple as specifying the executable name and the parameter file. Assuming that the ``lethe-rpt-3d`` executable is within your path, the simulation can be launched by typing:

.. code-block:: text
  :class: copy-button

  lethe-rpt-3d rpt-count-calculation.prm
  
Lethe will generate a ``.csv`` file with the name specified next to the ``counts file`` parameter in ``rpt-count-calculation.prm``. The generated ``.csv`` file will contain the :math:`(x,y,z)` coordinates of the particle with its respective photon count measured by a given detector. Each detector is identified by its id number (``detector_id``). In this example, as we have only one detector, all values in the ``detector_id`` column should be :math:`0`.


.. warning::
    When running the code with different particle position files, don't forget to change the name of the exporting ``counts file`` in ``rpt-count-calculation.prm`` so that the previous ``.csv`` file isn't overwritten.


-----------------------
Results and Discussion
-----------------------

To visualize the data and obtain the figures shown below, a Python script (``rpt_count-calculation_plot.py``) is provided. When running the script, the name of the ``.csv`` file you wish to open and read must be specified as an argument.

.. tip::
    You may use the ``rpt_count-calculation_plot.py`` script to plot any other set of data saved in a ``.csv`` file format.


Scenario 1: Horizontal Translation of a Particle along the X-Axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: images/counts-along-x-axis.png
    :alt: Results for the horizontal translation of a particle along the x-axis (Scenario 1)
    :align: center
    :name: Results for the horizontal translation of a particle along the x-axis 


In the figure shown above, as one would expect, as the particle approaches the detector, the photon count grows. Such evolution may be explained by the efficiency of the detector getting greater as the particle advances toward the detector's exposed surface. Since the photon's path length in the vessel decreases, :math:`f_a(\alpha, \theta)` increases, and therefore the efficiency gets greater. In addition to that, as the particle approaches the detector, the solid angle gets greater, the product :math:`\omega(\alpha) \omega(\theta)` increases, making the efficiency increase also.

Scenario 2: Horizontal Translation of a Particle along the Y-Axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: images/counts-along-y-axis-case1.png
    :alt: Results for the horizontal translation of a particle along the y-axis results when reactor attenuation coefficient is set at 10 and detector attenuation coefficient is set at 21.477 (Scenario 2)
    :align: center
    :name: Results for the horizontal translation of a particle along the y-axis (Case I)

    Case I: :math:`\mu_r = 10, \ \mu_d = 21.477`


The figure shown above illustrates the photon count of the particle as it travels from the back to the front of the vessel along the y-axis. The Case I figure shows the evolution of the photon count for the system we are currently studying (:math:`\mu_r = 10, \ \mu_d = 21.477`). Let's analyze the resulting plot.

First, a symmetry of photon counts from the center axis of the detector can be seen. Such symmetry should be expected since the detector is symmetrical from its center axis.

Secondly, we can notice that the variation in photon count as the particle travel is quite small. The difference between its maximal and minimal values is approximately :math:`147`, which is one order of magnitude smaller than the other scenarios. This may mainly be explained by the small variations in the distance between the particle and the detector's exposed surface. In other words, the lengths of the paths traveled by the photon in the vessel and in the detector vary less than in the other scenarios.

Lastly, as the particle travels across the vessel, we notice fluctuations in the photon count. Starting from the back of the vessel, the photon count decreases rapidly until a local minimal value at approximately :math:`y = -6` cm and then increases until a local maximum at :math:`y = 0` cm (center of the detector's face). Then, from the center to the front of the vessel, a mirrored image of the photon count's evolution can be seen. To understand the fluctuations, let's look at three other figures (Case II, Case III, and Case IV) while focusing on the first half of the studied domain (:math:`y \in ]-10, 0]` cm) since the evolution of the count is symmetrical from :math:`y = 0` cm.

+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------+
|  .. figure:: images/counts-along-y-axis-case2.png                                                       |   .. figure:: images/counts-along-y-axis-case3.png                                                      |
|    :alt: Results for the horizontal translation of a particle along the y-axis when the efficiency of   |     :alt: Results for the horizontal translation of a particle along the y-axis when                    |
|       the detector is the product of the weighting factors; fa and fd are constant and tend to 1        |         reactor attenuation coefficient is set at 0; fa is fixed to 1  (case III)                       |
|       (case II)                                                                                         |     :align: center                                                                                      |
|    :align: center                                                                                       |     :name: Results for the horizontal translation of a particle along the y-axis (case III)             |
|    :name: Results for the horizontal translation of a particle along the y-axis (case II)               |                                                                                                         |
|                                                                                                         |     Case III: :math:`\mu_r = 0, \ \mu_d = 21.477`                                                       |
|    Case II: :math:`\mu_r = 0, \ \mu_d = 1e9`                                                            |                                                                                                         |
|                                                                                                         |                                                                                                         |
+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------+
|  .. figure:: images/counts-along-y-axis-case4.png                                                       | .. figure:: images/reactor-path-lengths.png                                                             |
|    :alt: Results for the horizontal translation of a particle along the y-axis when detector attenuation|     :alt: Reactor path lengths for the horizontal translation of a particle along the y-axis            |
|        coefficient is set at 1e9; fd tends to 1 (case IV)                                               |     :align: center                                                                                      |
|    :align: center                                                                                       |     :name: Reactor path lengths for the horizontal translation of a particle along the y-axis           |
|    :name: Results for the horizontal translation of a particle along the y-axis (case IV)               |                                                                                                         |
|                                                                                                         |     :math:`e(\alpha, \theta)` function of :math:`y`                                                     |
|    Case IV: :math:`\mu_r = 10, \ \mu_d = 1e9`                                                           |                                                                                                         |
|                                                                                                         |                                                                                                         |
+---------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------+

The Case II figure shows the evolution of the photon count in absence of attenuation due to the medium found inside the vessel and the vessel's wall, and in the absence of variation of the interaction between the emitted :math:`\gamma`-ray and the detector. By setting :math:`\mu_r = 0`, we set :math:`f_a(\alpha_j, \theta_j) = 1`. As a consequence, the count becomes independent of the path of the photon inside the vessel. In a similar manner, by setting :math:`\mu_d = 1e9`, we make :math:`f_d(\alpha_j, \theta_j)` tend to :math:`1`. Consequently, the path traveled by the photon in the detector doesn't affect the efficiency anymore. Only the weighting factors :math:`\omega(\alpha)` and :math:`\omega(\theta)` have an influence on the calculated efficiency and photon count :math:`(\xi_i \approx \omega(\alpha) \omega(\theta))`. Therefore, the Case II figure gives us an idea of how the photon count evolves according to the particle's position respective to the detector's position disregarding the interactions between the emitted ray and the medium inside the vessel and its walls, and disregarding the interactions between the ray and the detector. We can use this case as a base to understand the interactions that occur in other cases.

The Case III figure depicts the evolution of the photon count in absence of the attenuation due to the medium found inside the vessel and the vessel's wall. Since we use the same set of positions in all cases, :math:`\omega(\alpha)` and :math:`\omega(\theta)` remain the same for each given position of the tracer particle. The length of the path traveled by the photon inside the detector should also be the same since the same seed number is used. As seen on the Case III figure, when the particle is aligned with the axis of symmetry of the detector, the photon count reaches a maximum. At that position, the evolution of the product :math:`\omega(\alpha) \cdot \omega(\theta)` seen on the Case II figure also reaches a maximum. And the distance :math:`d(\alpha,\theta)` reaches a local maximum at that position. On the case III figure, we notice that the inflection points at :math:`y \approx -5.5` cm and at :math:`y \approx -3.7` cm (not too far from the edge of the detector's face), seen on the Case II figure, are not present anymore. This means that when :math:`y \in ]-10, -3.8[` cm, when the particle sees both the face and the lateral sides of the detector and as the particle approaches the detector's face, the distance :math:`d(\alpha,\theta)` increases making the count increase. And when :math:`y \in ]-3.8, -1.5[` cm the distance :math:`d(\alpha,\theta)` decreases in such way that it counters the rapid increase in weighting factors giving the evolution of the photon count a more parabolic shape. Finally, between :math:`y \in ]-1.5, 0]` cm, :math:`d(\alpha,\theta)` increases until reaching a local maximum.

The last case studied (Case IV) shows the evolution of the photon count when :math:`\mu_d` is so great that :math:`f_d(\alpha_j, \theta_j)` tends to :math:`1 \ \forall y \in ]-10, 10[` cm. By doing so, we can see the evolution of the count when the efficiency is independent of the interaction between the emitted :math:`\gamma`-ray and the detector. With this case, we isolate the effect of the evolution of :math:`f_a(\alpha, \theta)` on the count. More specifically, we're looking at the evolution of :math:`e(\alpha,\theta)` as the particle travels in the vessel, since :math:`\mu_r` remains constant in the studied domain. We notice that we have a local minimum at :math:`y \approx -4.6` where we saw the convex section on the Case II figure. Considering the Case II results, we can interpret the Case IV figure as follows. Starting from the back of the vessel, where :math:`f_a(\alpha, \theta)` is at its maximal value, :math:`f_a(\alpha, \theta)` decreases at a decreasing rate until reaching :math:`y \approx -4.6` cm. The maximal value of :math:`f_a(\alpha, \theta)` (minimal value of :math:`e(\alpha,\theta)`) being when the particle is the furthest away from the detector may be explained by the curvature of the vessel's wall. Since the wall of the vessel is curved to form a circle, the distance traveled by the photon inside the vessel on the average probable path isn't necessarily larger than the radius of the reactor. We know that at :math:`y = 0`, :math:`e(\alpha,\theta) = 10` cm. In other words, :math:`e(\alpha,\theta)` is equivalent to the radius of the reactor. On the :math:`e(\alpha,\theta)` *function of* :math:`y` figure, we can read :math:`e(\alpha,\theta) \approx 10.04` cm when :math:`y = 10` cm. We also know that an increasing distance :math:`e(\alpha,\theta)` leads to a decreasing efficiency, which means a decreasing count. Therefore, we may assume that :math:`e(\alpha,\theta)` is minimal when :math:`y \approx -10` cm or when :math:`y \approx 10` cm. And, it slowly increases until reaching :math:`y \approx -4.6` cm. When the particle reaches the :math:`y \approx -4.6` cm position (local minimum), the variation of :math:`f_a(\alpha, \theta)` is so little that :math:`f_a(\alpha, \theta)` behaves as a constant. This explains why we see the same pattern of evolution of the photon count as in Case II when :math:`y \in ]-4.6, -3.8[` cm. Similarly, when the particle sees only the face of the detector, the pattern of the counts evolution follows the same trend as the one seen on the Case II when :math:`y \in ]-3.8, 0]` cm. This also indicates very little fluctuations of :math:`e(\alpha,\theta)` as we may see on the :math:`e(\alpha,\theta)` *function of* :math:`y` figure. Therefore, the photon count is highly dependant of the weighting factors when :math:`y \in ]-3.8, 0]` cm.

Coming back to the Case I figure, we can see that photon count follows a pattern similar to the one seen in Case IV. We may interpret from it that :math:`f_d(\alpha, \theta)` varies very little as opposed to :math:`f_a(\alpha, \theta)` that fluctuates greatly. The local minimal values, in this case, are at :math:`y \approx -6` cm and :math:`y \approx 6` cm, as opposed to :math:`y \approx -4.6` cm and :math:`y \approx -4.6` cm for the fourth case. This is due to the change in the value of :math:`\mu_d`. :math:`f_d(\alpha,\theta)` function of :math:`y` increases at a slower rate, making the minimums further way from the center. To summarize, the fluctuations seen in the Case I figure is the result of the combined influence of the values of the attenuation coefficients, the variation of the path lengths of the photon in the vessel and the detector, and the evolution of the weighting factors.


Scenario 3: Vertical Translation of a Particle along the Z-Axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: images/counts-along-z-axis.png
    :alt: Results for the vertical translation of a particle along the z-axis (Scenario 3)
    :align: center
    :name: Results for the vertical translation of a particle along the z-axis
	

Similar to the first scenario, as the particle approaches the detector, we notice an increase in photon count. The photon count reaches its maximal value at around :math:`z = 7.1` cm, which is close to the center of the detector's face.


Scenario 4: Particle Going across the Vessel on a Diagonal Line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: images/counts-across-vessel-on-a-diagonal-line.png
    :alt: Results for the particle going across the vessel on a diagonal line (Scenario 4)
    :align: center
    :name: Results for the particle going across the vessel on a diagonal line
	

After analyzing the past three scenarios, we get much-expected results for this scenario. As seen in the first scenario, the photon count varies greatly with the :math:`x` coordinate of the position vector of the particle. That is because the path of the photon inside the vessel gets longer when :math:`x` gets smaller. In other words, the ray is more attenuated by the material inside the vessel before getting to the detector, therefore the photon count gets smaller. Consequently, even though the particle is further away from the detector if the :math:`x` coordinate of the tracer's position is closer to the detector's exposed surface, the photon count could get greater and that's what we see on the figure above for high :math:`z` values.

Sensitivity Analysis of the Monte Carlo Method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Looking back at the second scenario's results (Case A), we notice that the counts are a little scattered. This is caused by the stochastic nature of the Monte Carlo method. Increasing the number of Monte Carlo iterations (:math:`N`), generates much smoother results as seen in the Case C figure where we have multiplied :math:`N` by a factor of :math:`10`. By increasing :math:`N`, we're covering more of the solid angle, making the simulation more representative of the physical system. Therefore, we see a better continuity in the photon counts. In the Case B figure, :math:`N` was divided by a factor of :math:`10`. As expected, in this figure, we see much more scattering.

+---------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/counts-along-y-axis-case1.png                                                                                                             |
|    :alt: Results for the horizontal translation of a particle along the y-axis results when reactor attenuation coefficient is set at 10 and detector         |
|       attenuation coefficient is set at 21.477 (Scenario 2)                                                                                                   |
|    :align: center                                                                                                                                             |
|    :name: Sensitivity analysis when N = 100000 (Case A)                                                                                                       |
|    :scale: 60%                                                                                                                                                |
|                                                                                                                                                               |
|    Case A: :math:`N = 1e5`                                                                                                                                    |
|                                                                                                                                                               |
+-----------------------------------------------------------------------------+---------------------------------------------------------------------------------+
|  .. figure:: images/sensitivity-analysis-caseB.png                          | .. figure:: images/sensitivity-analysis-caseC.png                               |
|    :alt: Scenario 2 results when reactor N = 10000                          |     :alt: Scenario 2 results when N = 1000000                                   |
|    :align: center                                                           |     :align: center                                                              |
|    :name: Sensitivity analysis when N = 10000 (Case B)                      |     :name: Sensitivity analysis when N = 1000000 (Case C)                       |
|                                                                             |                                                                                 |
|    Case B: :math:`N = 1e4`                                                  |     Case C: :math:`N = 1e6`                                                     |
|                                                                             |                                                                                 |
+-----------------------------------------------------------------------------+---------------------------------------------------------------------------------+



-----------
References
-----------

`[1] <https://doi.org/10.1016/0029-554X(78)90081-2>`_ G. B. Beam, L. Wielopolski, R. P. Gardner, and K. Verghese, “Monte Carlo calculation of efficiencies of right-circular cylindrical NaI detectors for arbitrarily located point sources,” *Nucl. Instrum. Methods*, vol. 154, no. 3, pp. 501–508, Sep. 1978, doi: 10.1016/0029-554X(78)90081-2.

`[2] <https://doi.org/10.1016/0168-9002(94)91343-9>`_ F. Larachi, G. Kennedy, and J. Chaouki, “A γ-ray detection system for 3-D particle tracking in multiphase reactors,” *Nucl. Instrum. Methods Phys. Res. Sect. Accel. Spectrometers Detect. Assoc. Equip.*, vol. 338, no. 2, pp. 568–576, Jan. 1994, doi: 10.1016/0168-9002(94)91343-9.
