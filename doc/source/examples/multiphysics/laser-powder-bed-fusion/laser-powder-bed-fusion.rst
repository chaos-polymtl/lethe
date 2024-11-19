===================================
Static Irradiation of a Bare Plate
===================================

This example simulates the static irradiation of a Ti6Al4V bare plate. It is based on the experimental work of Cunningham *et al.* and includes the relevant phenomena involved in the laser powder bed fusion manufacturing process. 

****

--------
Features
--------

- Solver: ``lethe-fluid`` 
- Volume of fluid (VOF) and Heat Transfer (HT)
- Unsteady problem with phase change handled by an adaptive BDF2 time-stepping scheme
- Python scripts for data postprocessing

****

---------------------------
Files Used in This Example
---------------------------

All files mentioned below are located in the example's folder (``examples/multiphysics/static-irradiation``).

- Parameter file: ``static-irradiation.prm``
- Mesh file: ``mesh/bare-plate.msh``
- Postprocessing Python script to extract quatities of interest: ``static-irradiation.py``

****

-----------------------
Description of the Case
-----------------------

Laser powder bed fusion is a manufacturing process using a laser to selectively melt and consolidate, layer-by-layer, a metal powder. Simply, it corresponds to 3D printing with metal powder. The main laser-material interaction takes place at the melt pool scale. The melt pool dynamics involve multiple driving forces:

- phase change due to laser heating
- surface tension due to the small scale of the melt pool
- evaporative cooling and recoil pressure due to temperature reaching the boiling point

In this example, we are considering the laser irradiation of Ti6Al4V bare plate (without power) to study the melt pool dynamics. The following figure shows the case setup, which is based on the experimental work of Cunningham *et al.*.

+-------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/new-benchmark.png                                                                             |
|     :align: center                                                                                                |
|     :width: 1000                                                                                                  |
|     :name: Case setup                                                                                             |
|                                                                                                                   |
+-------------------------------------------------------------------------------------------------------------------+