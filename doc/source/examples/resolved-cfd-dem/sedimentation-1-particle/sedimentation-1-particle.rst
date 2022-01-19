==============================================================================
Sedimentation of one particle.
==============================================================================

This example aims to reproduce numerically the results obtain by Ten Cate `.et al` `[1] <https://doi.org/10.1063/1.1512918>`_ for the E4 experience. This experience mesurethe velocity of the sedimentation of a 1.5 cm particle in a container filled with a viscous fluid. The container is sufficiently small to impact the particle sedimentation.


.. warning:: 
    This case is heavier then most exemple. It can take several hours to run.
    

Features
----------------------------------
- Solvers: ``gls_sharp_navier_stokes_3d`` (with Q1-Q1)
- Transient problem
- Displays the apability of the resolved cfd-dem solver for the flow around one particle

Description of the case
-----------------------
The E4 experiment consiste in the release of a 1.5 cm particle of Nylon (:math:`\rho_p=0.001120 \frac{kg}{cm^3}`) with it center 12.75 centimeter above the bottom of a 10x10x16 cm container. The viscosity of the fluid is :math:`\mu_f=0.00058 \frac{kg}{s cm}` which is equivalent to :math:`\mu_f=0.058 \frac{N s}{m^2}`. The density of the fluid is :math:`\rho_f=0.000960 \frac{kg}{cm^3}`. The gravity is :math:`g= -981 \frac{cm}{s^2}`


.. note:: 
    You will note where that we have transform every length units in centimeter. The reason is that the particle is very close to 1 cm. Representing the problem in this way helps the linear solver as it avoid extremly small value in the matrix due to the volume of the cell behing expressed in m^3 instead of cm^3. 
    


