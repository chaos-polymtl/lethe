***********************************************
Unresolved CFD-DEM
***********************************************
The CFD-DEM solver templates are available in both 2D ("cfd_dem_coupling_2d") and 3D ("cfd_dem>_coupling_3d"). The following parameter subsections are defined similarly to the CFD parameter template: "Simulation Control", "FEM", "Physical Properties" which represent the properties of the fluid, "Mesh", "Initial Conditions", and "Boundary Conditions". The following parameter subsections are defined similarly to the DEM parameter template: "Model Parameters" and "Lagrangian Physical Properties" which represents the physical properticles of the solid phase (particles). Two new subsections are introduced for the Volume Averaged Navier Stokes (VANS) solver as well as for the coupled unresolved CFD-DEM simulations. These are "Void Fraction" and "CFD-DEM" subsections. Parameters that are only applicable to the CFD-DEM solver are mentioned in their definitions.

.. toctree::
   :maxdepth: 1

   void_fraction
   cfd_dem

