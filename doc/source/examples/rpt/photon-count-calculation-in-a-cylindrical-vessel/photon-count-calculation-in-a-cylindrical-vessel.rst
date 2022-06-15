==================================================
Photon Count Calculation in a Cylindrical Vessel
==================================================

In this example, using a Monte-Carlo technique, we perform the calculation of photon counts of a single radioctive particle that emits gamma-ray. The calculation is performed for a given set of positions inside a cylindrical vessel. The Monte-Carlo method allows us to estimate the photon counts of a particle at a given position inside the vessel with respect to a given dectector. 

.. image:: images/rpt3.png
	:alt: The geometry
	:align: center
	:name: geometry_description 

Features
----------------------------------
- Solver: ``rpt_3d``
- 
- 

Locations of files used in the example
---------------------------------------
- Parameter file: ``examples/rpt/count_calculation/rpt_count_calculation.prm``
- File containing detector positions: ``examples/rpt/count_calculation/positions.detector``
- File containing particle positions on a horizontal line:  ``examples/rpt/count_calculation/positions_horizontal.particle``
- File containing particle positions on a vertical line:  ``examples/rpt/count_calculation/positions_vertical.particle``
- File containing particle positions on a diagonal line:  ``examples/rpt/count_calculation/positions_diagonal.particle``

Description of the case
-------------------------
In this example, three different sets of particule positions are studied.

The illustration below depicts the geometry of the vessel and the detector:


Parameter file
----------------

RPT Parameters
~~~~~~~~~~~~~~~

In the subsection ``rpt parameters``, we define the values of the set of parameters that is necessary for the calculation of the counts using the Monte Carlo method.  Amoung these parmeters, we have, the name of the file in which is found a set of different positions of the particle inside the vessel (``particle position file``), the number of Monte-Carlo iterations (``monte carlo iteration``), the seed that is used to generate a random number (``random number seed``) and other parameters that describe the studied gamma-ray model. We also define the name of the file in which the counts for each position will be exported with the parameter ``counts file``. These common parameters used for the RPT simulutation are discribed in the `RPT parameters <https://github.com/lethe-cfd/lethe/wiki/RPT-PARAMETERS#RPT-parameters>`_  subsection in the `RPT Parameters documentation <https://github.com/lethe-cfd/lethe/wiki/RPT-PARAMETERS>` page.

.. code-block:: text

	# Listing of Parameters
	# ---------------------
	# --------------------------------------------------
	# RPT Monte Carlo technique
	#---------------------------------------------------
	subsection rpt parameters
		set particle positions file          = positions.particle
		set verbosity                        = verbose
		set export counts                    = true
		set counts file                      = counts.csv
		set monte carlo iteration            = 100000
		set random number seed               = 0
		set reactor radius       			 = 0.1
		set peak-to-total ratio  			 = 0.4
		set sampling time        			 = 1
		set gamma-rays emitted        		 = 2
		set attenuation coefficient detector = 21.477
	end


Detector Parameters
~~~~~~~~~~~~~~~~~~~~

In the subsection ``detector parameters``, we specify the file that contains two positions located on the axis of symmetry of the detector. The first point is on the surface facing the vessel (face of the detector) and the second point can be any point located inside the detector. In the current example, the center position of the face is [0.15, 0, 0.08] and the second point located on the axis is [0.17, 0, 0.08]. We also specify the radius (``radius``) and the length (``length``) of the detector. Detailed description of these parameters can be found in the `Detector Parameters <https://github.com/lethe-cfd/lethe/wiki/RPT-PARAMETERS#Detector-parameters>`_ subsection in the `RPT Parameters domentation <https://github.com/lethe-cfd/lethe/wiki/RPT-PARAMETERS>` page.

.. code-block:: text

	#---------------------------------------------------
	# Detector parameters
	#---------------------------------------------------
	subsection detector parameters
		set detector positions file          = positions.detector
		set radius       			    	 = 0.0381 
		set length					         = 0.0762
		set dead time       				 = 1e-5
		set activity  						 = 2e6
		set attenuation coefficient reactor  = 10
	end

.. note::
	The parameters ``dead time``, ``activity`` and ``attenuation coefficient reactor`` are obtained using the black-box optimization software `NOMAD <https://www.gerad.ca/en/software/nomad/>`_ . The second example `"insert name here" <>` explains how we obtain the values of those parameters using NOMAD.






Results
--------



References
-----------
