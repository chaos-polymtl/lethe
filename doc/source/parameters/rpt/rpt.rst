***********************************************
Radioactive Particle Tracking (RPT)
***********************************************

Radioactive particle tracking (RPT) is a non-intrusive velocimetry method that is used to study the hydrodynamics of single and multiphase systems. Launching a RPT simulation for photon count calculation in Lethe requires a solver which is ``lethe-rpt-3d`, a parameter file, and two files including detector and particle positions inside the vessel. The ``particle positions file`` includes either the experimental calibration positions or a set of generated points inside the vessel by the user. The ``detector positions file`` contains the position of the detector face's center and the position of a point inside the detector on its axis. In Lethe-RPT code the middle of the bottom face of the cylindrical vessel is considered as the origin. An example of these two files can be found `here <https://github.com/lethe-cfd/lethe/tree/master/examples/rpt/count_calculation>`_. This section aims at describing the various parameters available within Lethe-RPT.

The parameter file is composed of different subsections. The principal subsections of a RPT parameter template file are :

.. toctree::
   :maxdepth: 1
   
   detector_parameters
   fem_reconstruction
   parameter_tuning
   rpt_parameters
