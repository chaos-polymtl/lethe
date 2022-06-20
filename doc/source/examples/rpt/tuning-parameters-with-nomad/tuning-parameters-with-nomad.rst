==================================================
Tuning Parameters with NOMAD
==================================================



Features
----------------------------------
- Solver: ``rpt_3d``
- 

Locations of files used in the example
---------------------------------------
- Parameter file: ``examples/rpt/count_calculation/rpt_count_calculation.prm``
- 

Description of the case
-------------------------
In this example,


The illustration below depicts the geometry of the vessel, the detector and the particle positioning for each case:

<insert image>
 <!-- ..image:: images/rpt_cases.png
	:alt: The geometry
	:align: center
	:name: geometry_description -->
	
As a particle travels in the cylindrical vessel, its photon count (:math:`C`) mesured by the detector varies according to the following relation:

.. math::
	\begin{align}
	C = \frac{T \nu R \phi \xi_i (\vec{X})}{1 + \tau \nu R \phi \xi_i (\vec{X})}
	\end{align}
		
In the previous expression, 

- :math:`T` is the sampling time (:math:`s`);
- :math:`\nu` is the number of :math:`\gamma`-rays emmited by each disintegration;
- :math:`R` is the activity of the tracer (:math:`Beq`);
- :math:`\phi` is the peak-to-total ratio;
- :math:`\tau` is the dead time of the detector (:math:`s`) and 
- :math:`\xi_i(\vec{X})` is the efficiency of the :math:`i_{th}` detector related to the position :math:`\vec{X}`.

The efficiency of the detector

Parameter file
----------------

RPT Parameters
~~~~~~~~~~~~~~~

In the subsection *rpt parameters*, we define the values of the set of parameters that is necessary for the calculation of the counts using the Monte Carlo method.  Amoung these parmeters, we have, the name of the file in which is found a set of different positions of the particle inside the vessel (``particle position file``), the number of Monte Carlo iterations (``monte carlo iteration``), the seed that is used to generate a random number (``random number seed``) and other parameters that describe the studied gamma-ray model. We also define the name of the file in which the counts for each position will be exported with the parameter ``counts file``. These common parameters used for the RPT simulutation are discribed in the **RPT Parameters**  subsection in the `RPT parameters <../../../parameters/rpt/rpt.html>`_ documentation page.

.. code-block:: text

	subsection rpt parameters
		set particle positions file           = positions_diagonal.particle
		set verbosity                         = verbose
		set export counts                     = true
		set counts file                       = counts_diagonal.csv
		set monte carlo iteration             = 100000
		set random number seed                = 0
		set reactor radius                    = 0.1
		set peak-to-total ratio               = 0.4
		set sampling time                     = 1
		set gamma-rays emitted                = 2
		set attenuation coefficient detector  = 21.477
	end


Detector Parameters
~~~~~~~~~~~~~~~~~~~~

In the subsection *detector parameters*, we specify the file that contains two positions located on the axis of symmetry of the detector. The first point is on the surface facing the vessel (face of the detector) and the second point can be any point located inside the detector. In the current example, the center position of the face is [0.15, 0, 0.08] and the second point located on the axis is [0.2, 0, 0.08]. We also specify the radius (``radius``) and the length (``length``) of the detector. Detailed description of these parameters can be found in the **Detector Parameters** subsection in the `RPT parameters <../../../parameters/rpt/rpt.html>`_ documentation page.

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
	The parameters ``dead time``, ``activity`` and ``attenuation coefficient reactor`` are obtained using the blackbox optimization software `NOMAD <https://www.gerad.ca/en/software/nomad/>`_ . The second example `Tuning Parameters with NOMAD <insert link>`_ explains how we can obtain the values of those parameters using NOMAD.





Results
--------


Senario 1 : Horizontal translation of a particle on the x-axis 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Senario 2 : Horizontal translation of a particle on the y-axis 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Senario 3 : Vertical translation of a particle 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Senario 4 : Particle going across the vessel on a diagonal line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


References
-----------
