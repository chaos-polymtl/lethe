***********************************************
Radioactive Particle Tracking (RPT)
***********************************************

Radioactive particle tracking (RPT) is a non-intrusive velocimetry method that is used to study the hydrodynamics of single and multiphase systems. Launching a RPT simulation for photon count calculation in Lethe requires a solver which is ``rpt_3d``, a parameter file, and two files including detector and particle positions inside the vessel. The ``particle positions file`` includes either the experimental calibration positions or a set of generated points inside the vessel by the user. The ``detector positions file`` contains the position of the detector face's center and the position of a point inside the detector on its axis. In Lethe-RPT code the middle of the bottom face of the cylindrical vessel is considered as the origin. An example of these two files can be found `here <https://github.com/lethe-cfd/lethe/tree/master/examples/rpt/count_calculation>`_. This section aims at describing the various parameters available within Lethe-RPT.

Parameters file (.prm)
-----------------------

In the parameter file format, the parameters are established one by one using the following syntax, for instance: ``set reactor radius = 0.1`` would fix the reactor radius to 0.1 m. The arguments can be either doubles, integers, or a choice between a predefined condition. In the parameter files, comments are preceded by the sharp symbol (e.g., *#comment*).

The parameter file is composed of different subsections. The principal subsections of a RPT parameter template file are :

.. toctree::
   :maxdepth: 1
   
   rpt_parameters
   detector_parameters
   parameter_tuning
