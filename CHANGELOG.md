
# Change Log
All notable changes to the Lethe project will be documented in this file.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

### [Master] - 2025-07-25

### Added

- MINOR A feature that allows to log the particle-wall collision statistics in a .csv or a .dat file is added. The time, velocity and angular velocity at the start and at the end of a collision are logged for each collision between a particle and the boundary id(s) specified in the .prm file. [#1586](https://github.com/chaos-polymtl/lethe/pull/1586)

### Added

- MINOR The DG tracer implementation was not very good for pure advection equations since the high-order implementation did not have a limiter. This PR implements a Moe limiter that works for all scalar equations with DG. This is especially useful for the DG Tracer implementation and will be tested for the DG VOF implementation in the future. [#1595](https://github.com/chaos-polymtl/lethe/pull/1595)

### [Master] - 2025-07-23

### Fixed

- MINOR The default insertion frequency was causing a segmentation fault when using the DEM solver (division by zero when set unchanged). This PR adds an "if" statement that checks if this insertion frequency is equal to 0. [#1593](https://github.com/chaos-polymtl/lethe/pull/1593)

### [Master] - 2025-07-21

### Changed

- MINOR The QCM calculation method for the void fraction had multiple chunks of code which were duplicated. This PR refactors the QCM to make an additional use of functions and greatly simplifies the code. Nothing is changed in the way the calculations are made. [#1591](https://github.com/chaos-polymtl/lethe/pull/1591)

### [Master] - 2025-07-20

### Fixed

- MAJOR Deal.II-9.8pre deprecates some functions such as get_communicator and the parallel::distributed::SolutionTransfer class. This PR uses the modern version of these functions or classes following the new 9.8 release standard before they are deprecated. This PR officially deprecates the compatibility with deal.II 9.6.0.  [#1589](https://github.com/chaos-polymtl/lethe/pull/1589)

### Changed

- MINOR At multiple locations and for specific features, Lethe was still checking if deal.II was compiled with a version number larger than 9.4. Since we do not even compile with anything before deal.II 9.6, this PR removes all of these checks which are not relevant anymore. This increases readability and prevents confusion in terms of compatibility. [#1588](https://github.com/chaos-polymtl/lethe/pull/1588)

### [Master] - 2025-07-17

### Added

- MINOR Allow quad/hex meshes to use bubble enrichment function FE_Q_Bubbles. [#1584](https://github.com/chaos-polymtl/lethe/pull/1584)

### [Master] - 2025-07-16

### Changed

- MINOR Lethe could not be compiled when deal.II was using 64 bit indices due to some locations where unsigned int were used instead of types::global_dof_index. This PR fixes this by using the correct type. It compiles with deal.II with 64bit indices. [#1585](https://github.com/chaos-polymtl/lethe/pull/1585)

### [Master] - 2025-07-14

### Changed

- MAJOR The introduction of the deal.II-9.8-pre version introduced numerous warnings for functions that will be deprecated in the near future. This PR addresses all of these warnings. Furthermore, I reran all of the prototypes and realized there was an issue with the matrix_free_non_linear prototype in which the multigrid solution did not have their ghost value correctly calculated. I have addressed this issue at the same time. Finally, on my machine I encountered a small segfault when running the kelly error estimator with multiple variables. This is because the solver assumed the number of coarsen and refine flag would be the same in the for loop, while there is no reason for this assumption. I have thus fixed this problem within this PR also. [#1582](https://github.com/chaos-polymtl/lethe/pull/1582)

### [Master] - 2025-07-13

### Changed

- MAJOR The void fraction linear tolerance can now be specified directly within the .prm file. The previous default value which was 1e-15 has now been changed to the default value used for all linear solvers. We advise users to explicitely specify the minimum tolerance for the void fraction within their parameter files. [#1580](https://github.com/chaos-polymtl/lethe/pull/1580)

### [Master] - 2025-07-11

### Fixed

- MINOR The rotor mesh rotation was not in the right place in the code, so that at the start of the simulation the accessed mapping information was still incorrect. This PR fixes this, and adds an application test for a rotor rotated mesh. [#1579](https://github.com/chaos-polymtl/lethe/pull/1579)

### Added

- MINOR Removes the old gas-solid-spouted-bed example which has now been replaced with the much better gas-solid-spouted-rectangular-bed. The newer example features validation and other components which make it significnatly more relevant than the previous one which just displayed the capabilities of the solver. [#1578](https://github.com/chaos-polymtl/lethe/pull/1578)

### Added

- MINOR Dynamic particle insertion in CFD-DEM functionality is added using the insertion methods defined for DEM. [#1576](https://github.com/chaos-polymtl/lethe/pull/1576)

### [Master] - 2025-07-10

### Added

- MAJOR A new example consisting in the 3d simulation of the turbulent flow around a cylinder at Re=3900 has been added. The example includes several comparison with the literature including the drag coefficient, the Strouhal number and the pressure profile around the cylinder for several meshes and orders of convergence. The example will be extended in the near future to include additional validation metric such as the time-averaged velocity and reynolds stresses at serveral locations behind the cylinder. [#1570](https://github.com/chaos-polymtl/lethe/pull/1570)

# [Master] - 2025-07-08

### Fixed

- MINOR The center of rotation of the rotor was not being accounted for in the MortarManagerBase, which has been fixed. Two application tests for a MMS with mortar have also been added. [#1575](https://github.com/chaos-polymtl/lethe/pull/1575)

## [Master] - 2025-07-03

### Added

- MINOR Add tailored capability for validation script to run on Lucille (a compute node) and re-adapt rising bubble case. [#1574](https://github.com/chaos-polymtl/lethe/pull/1574)

## [Master] - 2025-07-07

### Changed

- MINOR The default parameter for the linear solver of the matrix-free solvers now uses eigenvalue estimation instead of not. This approach is generally more robust and should be preferred. [#1572](https://github.com/chaos-polymtl/lethe/pull/1572)

## [Master] - 2025-07-03

### Added

- MINOR Add functionality that allows to delete previous vtu and pvd files by using a flag '-R'. [#1569](https://github.com/chaos-polymtl/lethe/pull/1569)

## [Master] - 2025-06-30

### Added

- MINOR The surface area computation of the interface was added in this PR for VOF simulations when mass conservation monitoring is enable. For that, the NonMatching FEValues features of deal.II are used. Additionally, the geometric volume computation was changed to use the same features. [#1562](https://github.com/chaos-polymtl/lethe/pull/1562)

## [Master] - 2025-06-26

### Added

- MINOR This PR adds an Multiphysic DEM example that simulates the heating of a packed bed. It is based on the experimental and simulation results of Beaulieu with stainless steel particles. [#1561](https://github.com/chaos-polymtl/lethe/pull/1561)

## [Master] - 2025-06-25

### Added

- MAJOR This PR implements the coupling of the mortar method with the Navier-Stokes matrix-based solver. An application test of a Taylor-Couette flow has been added, where a mortar interface binds two hyper_shell geometries. [#1563](https://github.com/chaos-polymtl/lethe/pull/1563)

## [Master] - 2025-06-19

### Added

- MINOR This PR extends the rising bubble example with the addition of a second test case and the comparison of the regularization methods implemented in Lethe. It also documents the geometric reinitialization method in the parameters section and the theory guide. [#1558](https://github.com/chaos-polymtl/lethe/pull/1558)

## [Master] - 2025-06-16

### Fixed

- MINOR The function defining zero_constraints considered the types noslip, function, and none in the same else case. This has been fixed so that none is indeed a "do-nothing" boundary condition. [#1557](https://github.com/chaos-polymtl/lethe/pull/1557)

### Added

- MAJOR The geometric redistanciation method had only a tanh transformation to go from the signed distance to the phase fraction. With it, it is difficult (technically impossible) to clamp the value of the phase fraction indicator to 0 and 1 away from the interface. To get close, the max redistanciation distance must be increased leading to higher computational cost. Hence, a 4th degree piecewise polynomial transformation is added in this PR. It clamps the phase fraction to 0 or 1 at +/- the max redistanciation distance and it has a smooth change from 0 to 1. This PR also improves the treatment of the level-set field and the iso-level in the VOF and SignedDistanceSolver implementation. [#1546](https://github.com/chaos-polymtl/lethe/pull/1546)

## [Master] - 2025-06-15

### Fixed

- MINOR Add capability to calculate the void fraction in the matrix-free vans solver using information drawn from a particle handler instead of just a function. This is a minor extension that re-uses the feature of the regular VANS implementation. [#1556](https://github.com/chaos-polymtl/lethe/pull/1556)

## [Master] - 2025-06-13

### Fixed

- MINOR Fixed bug with constrained indices and edge constrained indices in the functions required by the operator when using local smoothing geometric multigrid in the matrix-free application. Added a cylinder test with kelly error estimator to avoid this in the future. [#1555](https://github.com/chaos-polymtl/lethe/pull/1555)

## [Master] - 2025-06-12

### Fixed

- MINOR Fixed restarts in CFD-DEM simulations with Lagrangian post-processing enabled. The checkpointing file for the pvd file was not named correctly. This has now been corrected to ensure consistency. [#1552](https://github.com/chaos-polymtl/lethe/pull/1552)

### Added

- MINOR Added a simple first multiphysic DEM example that consists of lined up particles that are heated by others or by a wall. It verifies the models of particle-particle and particle-wall heat transfer using two simple test cases. This PR also adds dem-mp in the documentation and examples directories. [#1548](https://github.com/chaos-polymtl/lethe/pull/1548)

- MINOR Added initial support for the dynamic flow control in lethe-fluid-sharp. Beta force particle is still missing from the calculation. [#1553] (https://github.com/chaos-polymtl/lethe/pull/1553)

## [Master] - 2025-06-11

### Fixed

- MINOR Added support for the box refinement in the matrix-free solver. [#1551](https://github.com/chaos-polymtl/lethe/pull/1551)

## [Master] - 2025-06-05

### Added

- MAJOR Added the calculation of particle-wall heat transfer and the updating of the temperature of solid objects when needed.
A parameter real_youngs_modulus is added, to be used for Multiphysic DEM when the Young's modulus is underestimated in the simulation.
A parameter disable_position_integration is also added to be able to freeze the position of particles. [#1542](https://github.com/chaos-polymtl/lethe/pull/1542)

## Release of Lethe v1.0.1 - 2025-05-28
The Lethe v1.0.1 release introduces some key new features including:

A new DEM rolling friction model
- Two new VOF redistanciation/sharpening algorithm based on geometric and algebraic redistantiation
- Improvement to the matrix-free algorithm enabling their coupling with heat transfer, but also non-Newtonian flows
- Minor bug fixes and improvement

The Lethe v1.0.1 release comes with a structured and stable syntax. All of the examples and test provided with lethe were re-ran with this v1.0.1 release to ensure that all of the results provided online are up to date. This version is production-ready and has been tested and validated extensively.

## [Master] - 2025-05-27

### Fixed

- MINOR The output of the mass conservation table in VOF was resulting in a seg fault in 3D. It was permanently fixed by adding the missing column name in 3D. [#1541](https://github.com/chaos-polymtl/lethe/pull/1541)

## [Master] - 2025-05-26

### Changed

- MINOR The rolling friction used in the hopper example as been change from 0.09 to 0.1776 to be equivalent to the one used in the reference article. This change will slightly affect the results from the validation script. [#1539](https://github.com/chaos-polymtl/lethe/pull/1539)

- MINOR The particle_particle_heat_transfer files are changed to calculate heat transfer for both particle-particle and particle-wall contacts with the same functions, as the calculation is very similar. [#1534](https://github.com/chaos-polymtl/lethe/pull/1534)

## [Master] - 2025-05-21

### Fixed

- MINOR The rotor mesh rotation in the mortar method context was being done by rotating the triangulation. This PR changes this so that the mapping is rotated instead. [#1536](https://github.com/chaos-polymtl/lethe/pull/1536)

## [Master] - 2025-05-20

### Fixed

- MINOR The functions used to define constraints in navier_stokes_base were calling reinit() without passing the locally_owned_dofs. This caused the constraints object to not have the updated locally owned DoFs in each process. This has now been fixed. [#1533](https://github.com/chaos-polymtl/lethe/pull/1533)

### Changed 

- MAJOR The architecture of the code for particle-wall contacts was changed to reproduce the one for particle-particle contacts and remove code duplication. 
The particle_wall_contact_info has been changed to a struct instead of a class and normal_overlap, tangential_relative_velocity and normal_relative_velocity_value are calculated on the fly and not stored in contact_info anymore. [#1520](https://github.com/chaos-polymtl/lethe/pull/1520)

## [Master] - 2025-05-09

### Added

- MAJOR Added the buoyancy term to the matrix free operators to be able to support two way coupling between the heat transfer solver and the matrix-free fluid solver. In the GMG preconditioner the temperature solution is also transfered to all the MG levels. [#1524](https://github.com/chaos-polymtl/lethe/pull/1524)

### Added

- MINOR Mortar: the matrix-based solver is now able to read 'gmsh' files, and a function to compute the rotor radius and the number of subdivisions at the rotor-stator interface has been added. [#1526] (https://github.com/chaos-polymtl/lethe/pull/1526)

### Added

- MAJOR Added the buoyancy term to the matrix free operators to be able to support two way coupling between the heat transfer solver and the matrix-free fluid solver. In the GMG preconditioner the temperature solution is also transfered to all the MG levels. [#1524](https://github.com/chaos-polymtl/lethe/pull/1524)

### Added

- MINOR Some of the functions related to the mortar coupling have been refactored, and new functionalities have been added. For instance, the matrix-based solver is able to read 'gmsh' files, and the matrix-free solver is able to read data from the mortar subsection. A new coupling operator has been added. [#] (https://github.com/chaos-polymtl/lethe/pull/)

## [Master] - 2024-05-07

### Fixed

- MINOR The solution preceding the algebraic interface reinitialization time-step is now also regularized using the same process to ensure consistency in the time integration scheme. This is done by replacing the relevant previous solution with the reinitialized previous solution. A new application-test (`vof-algebraic-interface-reinitialization-advected-circle-bdf2-frequency-check`) has been added to test this implementation. Also, a new validity map has been added to the VOFSubequationInterface. This ensures that subequations are solved using the correct VOF solution and in the right order. Subequations dependency checks are also made before solving each subequation. [#1507](https://github.com/chaos-polymtl/lethe/pull/1507)

## [Master] - 2025-04-30

### Added

- MAJOR Multiphysic DEM simulations can now be launched with solver type dem_mp to follow the evolution of the temperature of the particles. The temperature is also added in the properties printed by application tests for the solver type dem_mp. [#1514](https://github.com/chaos-polymtl/lethe/pull/1514)

## [Master] - 2025-04-25

### Changed 

- MINOR The find_cell_neighbors is now templated with a boolean so that neighboring cell vectors can be reciprocal, which means that vector i and j will contain cell j and i, respectively. [#1512](https://github.com/chaos-polymtl/lethe/pull/1512)

### Fixed

- MINOR The Jacobian of the non-Newtonian GLS assembler has an error in the advection term. This has been fixed. [#1513](https://github.com/chaos-polymtl/lethe/pull/1513)

## [Master] - 2025-04-24

### Added

- MINOR This PR Adds the find_line_sphere_intersection in lethe_grid_tools for use in the case of ray-particle intersection for the profilometry hackathon project. This function simply evaluates all intersection points between a line and a sphere. [#1511](https://github.com/chaos-polymtl/lethe/pull/1511)

### Added

- MAJOR Adds an experimental prototype that uses the DG method to solve the VOF equation. At the present stage, this prototype is compatible with surface tension, but does not support any of the redistanciation mechanism nor the smoothing of the initial condition. Furthermore, it does not yield satisfactory results when BDF2 is used as a time integration scheme. This feature is thus very experimental. [#1510](https://github.com/chaos-polymtl/lethe/pull/1510) 

## [Master] - 2025-04-24

### Changed

- MAJOR The parameters of the function which calculates the contact force, torque and heat transfer rate between particles were changed. Parameters torque, force and heat_transfer_rate were replaced by a class object ParticleInteractionOutcomes where they are stored. A function to resize these containers is also implemented in the class. [#1504](https://github.com/chaos-polymtl/lethe/pull/1504)

## [Master] - 2025-04-22

### Added 

- MINOR This PR completes the signed distance solver by filling the (almost) empty architecture of the SignedDistanceSolver implemented in [#1451](https://github.com/chaos-polymtl/lethe/pull/1451). Its coupling in the VOF solver (geometric reinitialization method) is also completed. The implemented routines follow the one of the already existing vof_advection prototype. [#1497](https://github.com/chaos-polymtl/lethe/pull/1497)

## [Master] - 2025-04-17

### Added

- MAJOR The template parameter PropertiesIndex was added to the insertion classes to initialize the temperature only in multiphysic DEM. The template parameter dim was added to InsertionInfo (Parameters Lagrangian) to initialize the temperature with a parsed function in the parameter file for plane and volume insertion, which were modified accordingly. [#1491](https://github.com/chaos-polymtl/lethe/pull/1491) 

## [Master] - 2025-04-16

### Added

- MAJOR The mortar feature from [#1462](https://github.com/chaos-polymtl/lethe/pull/1462) started to be integrated into the lethe-fluid application. Together with [#1481](https://github.com/chaos-polymtl/lethe/pull/1481), a new mortar subsection is now available the parameters file. Rotor and stator meshes generated using dealii grids are merged into a unique triangulation; coupling at the interface is not yet considered. [#1490](https://github.com/chaos-polymtl/lethe/pull/1490) 

## [Master] - 2025-04-14

### Changed

- MAJOR The "tangential overlap" terminology is a misused of language since an overlap means that two things are occupying the same physical space which is not the case this vector. For this reason, this terminology was replaced by "tangential displacement" every where in the code and in the documentation. [#1492](https://github.com/chaos-polymtl/lethe/pull/1492)

## [Master] - 2025-04-13

### Fixed

- MINOR In the sharp immersed boundaries solver, added a condition that clears the combined_shape cache if a particle has moved. This fixes a bug that would lead to an incorrect level set function when particles would move accross subdomains. This only affected the visualization of the particles. [#1479](https://github.com/chaos-polymtl/lethe/pull/1495)

- MINOR The Rayleigh-Benard convection example would sometimes not work correctly on some machines and instead generate a static solution. This has been fixed by adding an initial condition in which a linear temperature profile is imposed in the direction orthogonal to the gravity (along the x axis). This also increases the robustness of the example. [#1496](https://github.com/chaos-polymtl/lethe/pull/1496)

## [Master] - 2025-04-10

### Added

- MINOR Added an example demonstrating the use of the Method of Manufactured Solutions (MMS) with the incompressible solver, using quadrilateral and simplex meshes at varying resolutions. [#1461](https://github.com/chaos-polymtl/lethe/pull/1461)

## [Master] - 2025-04-07

### Added

- MINOR The GMG preconditioner in the matrix-free VANS solver now uses the void fraction on every level instead of only the finest one. [#1458](https://github.com/chaos-polymtl/lethe/pull/1458)

## [Master] - 2025-04-02

### Fixed

- MINOR In the sharp immersed boundaries solver, added a condition that checks if immersed solids have moved since the last time step before executing mesh refinement. [#1479](https://github.com/chaos-polymtl/lethe/pull/1479)

## [Master] - 2025-04-01

### Added

- MAJOR Added library for mortar feature, including tests using the Poisson operator and the Stokes operator. [#1462](https://github.com/chaos-polymtl/lethe/pull/1462)

## [Master] - 2025-03-31

### Added

- MAJOR Added a new feature in the DEM solver. The base weight of a cell is now defined by a parsing function in the parameter file. This allows for cells to have a personalized weight depending on the position of their barycenter in space and time. This feature will be useful when dealing with problems requiring few load balancing steps, like the granuheap example, where the tangential overlap is important and needs to be fully tracked. [#1446](https://github.com/chaos-polymtl/lethe/pull/1446)

- MINOR After a restart, the first time step is now automatically a load-balancing step, except if the load-balancing method used in "none".

## [Master] - 2025-03-28

### Added

- MINOR Added a functionality that allows to print the deal.II and Lethe versions when running an application via a '-V' flag. Also, a parameter that allows to print content of the input file to output is added. [#1475](https://github.com/chaos-polymtl/lethe/pull/1475)

### Fixed

- MINOR Phase fraction gradient projection is computed with the filtered phase fraction gradient. However, before the starting the VOF algebraic interface reinitialization process, the filter was not applied leading to wrong values of filtered phase fraction gradient and therefore curvature. This is now fixed. Also, in the calculate_momentum function of VOF postprocessing quantities, the FEValues of the fluid dynamics was initialized with the FE degree of the fluid dynamics physics instead of the VOF one leading to a vector size error when getting the function values. The quadrature formula has now been changed to the VOF one. [#1474](https://github.com/chaos-polymtl/lethe/pull/1474)

## [Master] - 2025-03-26

### Changed

- MAJOR Refactor the VOF interface sharpening and reinitialization parameters and merge them in one subsection called interface regularization method. [#1467](https://github.com/chaos-polymtl/lethe/pull/1467)

## [Master] - 2025-03-26

### Added

- MINOR Addition of the explicit template instantiations for the PropertiesIndex of multiphysic DEM simulations. They are added to all files where there were already DEMProperties template instantiations. [#1469](https://github.com/chaos-polymtl/lethe/pull/1469)
## [Master] - 2025-03-25

### Added

- MINOR Changed the simulation control to ensure that for transient simulation the last time-step is always written to a vtu file. Used this opportunity to improve the results for the rising bubble since the contour we were using to compare to the reference solution was not always the last output. Finally, changed the rising bubble example to use the BDF2 time integration instead of BDF1 [#1471](https://github.com/chaos-polymtl/lethe/pull/1471).

## [Master] - 2025-03-21

### Added

- MINOR Fix the simulation restart. When instantiating a SolutionTransfer object, the average_values flag can be set to true to average the contributions to the same DoF coming from different cells. This boolean must be set to true during a mesh adaptation, but not during a simulation restart. The average_values flag was therefore removed from the simulation restart. [#1466](https://github.com/chaos-polymtl/lethe/pull/1466)

## [Master] - 2025-03-20

### Added

- MINOR Addition of the periodic boundary conditions for the Cahn-Hilliard, heat transfer and tracer physics solvers. Parameters for periodic boundary conditions for all auxiliary physics (including VOF) now follow the same syntax as for the fluid dynamics solvers. [#1456](https://github.com/chaos-polymtl/lethe/pull/1456)

## [Master] - 2025-03-19

### Added

- MINOR Verify if non-Newtonian fluid parameters are used in the matrix-free solver. Since the MF solver does not support non-Newtonian fluids at the moment, an exception is now shown if this happens. [#1463](https://github.com/chaos-polymtl/lethe/pull/1463)

## [Master] - 2025-03-13

### Fixed

- MAJOR This PR switches the CI from the deal.ii-focal to the deal.ii-jammy image. This image uses Mold, which meant that some small refactor was necessary to prevent some compilation issues. Furthermore, this image uses P4est 2.2, which means that the restart files for the tests that require restart files had to be regenerated. [#1464](https://github.com/chaos-polymtl/lethe/pull/1464)

## [Master] - 2025-03-13

### Added

- MINOR This PR adds the (almost) empty architecture of the SignedDistanceSolver and its coupling in the VOF solver (geometric reinitialization method), following the structure of the already existing vof_advection prototype. The parameters associated to the geometric reinitialization method for the VOF solver are also added. [#1451](https://github.com/chaos-polymtl/lethe/pull/1451)

## [Master] - 2025-03-11

### Added

- MINOR Remove muParser and custom parsed function class (reversal of [#1143](https://github.com/chaos-polymtl/lethe/pull/1143), [#1153](https://github.com/chaos-polymtl/lethe/pull/1153), [#1160](https://github.com/chaos-polymtl/lethe/pull/1160)). [#1452](https://github.com/chaos-polymtl/lethe/pull/1452)

## [Master] - 2025-03-06

### Added

- MAJOR A first version of the lethe-fluid-vans-matrix-free solver is added. This version supports all of the preconditioners of the regular matrix-free solver. However, it is limited right now to void fraction fields which are defined using analytical functions and not from particles. Furthermore, it does not support any particle-fluid coupling. [#1448](https://github.com/chaos-polymtl/lethe/pull/1448)

## [Master] - 2025-03-03

### Added

- MINOR New functionalities for interface problems are added. The PR includes the computation of the volume enclosed by a given iso-level of a level-set field and the surface reconstruction of the interface described by the level 0. [#1429](https://github.com/chaos-polymtl/lethe/pull/1429)

## [Master] - 2025-02-28

### Fixed

- MINOR Simulation with dynamic mesh adaptation and the matrix-free solver could sometimes crash during the mesh adaptation stage due to the SolutionTransfer. This has been fixed by enabling averaging of the SolutionTransfer to reconcile incoherence between multiple cells. [#1436](https://github.com/chaos-polymtl/lethe/pull/1436)

- MINOR Simulation with dynamic mesh adaptation and time-averaging of the velocity field carried out using the lethe-fluid-matrix-free solver would crash with a thrown error. This was caused by the fact that solution transfer with deal.II vector requires ghosted vectors whereas solution transfer with Trilinos vectors requires locally_owned vectors. This has been fixed by adding a const expression to adapt the behaviour depending on the vector type. [#1436](https://github.com/chaos-polymtl/lethe/pull/1436)

## [Master] - 2025-02-27

### Added

- MINOR Add two parameters to the unresolved CFD-DEM simulations using the QCM void fraction scheme, allowing to choose the quadrature rule and the number of quadrature points used in the void fraction calculation. This implementation includes the Gauss quadrature rule (default) with 2 quadrature points as default and Gauss-Lobatto quadrature rule with 3 quadrature points as default. [#1432](https://github.com/chaos-polymtl/lethe/pull/1432)

## [Master] - 2025-02-24

### Fixed and added

- MINOR Temperature and heat flux time averaged quantities were computed at every output iteration instead of at every time iteration. A mechanism was also added to reset the averaged quantities to zero after a restart if the initial averaging time exceeds the simulation time. [#1431](https://github.com/chaos-polymtl/lethe/pull/1431)

## [Master] - 2025-02-24

### Fixed

- MINOR The function get_properties_name() did not contain the pair corresponding to "volumetric_contribution". As consequence, the CFD-DEM vtu files would not have this property and would have an empty property instead. This PR adds this modification to the changelog as the addition of the property was pushed straight to master. The hash or the commit that pushes it to master is: 3dd817096cd2f5056902ace464e1aef9358c4ef8. [#1430](https://github.com/chaos-polymtl/lethe/pull/1430)

## [Master] - 2025-02-22

### Added

- MINOR In the tracer physics, a new feature is added that allows to take into account single species reaction of arbitrary reaction constant and order. This feature works in CG or DG, with or without immersed solids. [#1411](https://github.com/chaos-polymtl/lethe/pull/1411)

## [Master] - 2025-02-20

### Added

- MINOR This PR adds the possibility to rotate the laser around one axis when using the gaussian_heat_flux_vof_interface model. [#1421](https://github.com/chaos-polymtl/lethe/pull/1421)

## [Master] - 2025-02-19

### Fixed

- MINOR The vtu output for 2D DEM and CFD-DEM simulations were not able to display the velocity and the angular velocity correctly. This is because the data is always 3D in the simulations even for 2D simulations. This has been fixed by always writing velocity and angular velocity as 3D vectors even for 2D simulations.[#1427](https://github.com/chaos-polymtl/lethe/pull/1427)

## [Master] - 2025-02-18

### Added

- MINOR A new prototype for the geometric redistanciation is added. It includes an advection solver and the geometric redistanciation itself.[#1418](https://github.com/chaos-polymtl/lethe/pull/1418)

### Fixed

- MAJOR Since [#1417](https://github.com/chaos-polymtl/lethe/pull/1417), the restart_moving_receptable test had a segfault at the end of the simulation in semi-reproducible ways. After using valgrind, it was found that the boost signal that was connecting the Insertion class and the triangulation was causing this segfault when the triangulation was destructed. This was because the Insertion class was destructed prior, but the boost signals were still connected. To solve this issue, the Insertion class now has a destructor that disconnect the boost signal before destructing the whole object. [#1425](https://github.com/chaos-polymtl/lethe/pull/1425)

## [Master] - 2025-02-17

### Added

- MINOR The neighboring threshold of some application test for the lethe-particles executable were using a value of 20 when it should be in the range of 1.0 to 1.5 . Concerned were updated to a value of 1.3 .  [#1427](https://github.com/chaos-polymtl/lethe/pull/1427)

### Added

- MINOR A new example showcasing heat transfer in a cooling fin.[#1426](https://github.com/chaos-polymtl/lethe/pull/1426)

## [Master] - 2025-02-15

### Added

- MAJOR A new rolling resistance model called the elastic-plastic spring-dashpot model has been added to the DEM. This model will help when dealing with static DEM problem, like in the realization of a heap of particle.[#1417](https://github.com/chaos-polymtl/lethe/pull/1417)

## [Master] - 2025-02-14

### Fixed

- MINOR When outputting vtu files from the algebraic reinitialization process, the PVDHandler was accumulating vtu files of previous reinitialization processes. This is now fixed by clearing the `times_and_names` vector of the PVDHandler object when setting the initial condition. [#1423](https://github.com/chaos-polymtl/lethe/pull/1423)

## [Master] - 2025-02-13

### Fixed

- MINOR The rotating drum example was updated to include a comparison with experimental results. The parameters of the example had regressed with time and we now obtain the same (or highly similar) results we obtained with the original Lethe-DEM version with the current master version. The main issue we had was related to the fact that the front and the back walls of the cylinders were not rotating in the rotating cylinders examples.[#1413](https://github.com/chaos-polymtl/lethe/pull/1422)

## [Master] - 2025-02-10

### Added

- MINOR To improve VOF interface definition through simulations, a new algebraic interface reinitialization process was added. To reinitialize the interface, a non-linear PDE is solved. Since all previous added subequations were linear, a new class for non-linear subequations was added, namely PhysicsNonlinearSubequationsSolver. This class inherits from PhysicsSolver and follows the same solving algorithm as other non-linear equations solved in Lethe. The algebraic reinitialization PDE is integrated to the VOFSubequationsInterface.[#1416](https://github.com/chaos-polymtl/lethe/pull/1416)

## [Master] - 2025-02-09

### Fixed

- MINOR Updated the installation instructions under Linux to specify installation since it was not previously explained. Took the opportunity to correct the graphviz for multiphysics problem. [#1419](https://github.com/chaos-polymtl/lethe/pull/1419)

## [Master] - 2025-01-24

### Added

- MINOR A new validation case is added. This new case used the full scale rotating drum example. All the validation cases can be launch using the contrib/validation/validate_lethe.sh script.[#1413](https://github.com/chaos-polymtl/lethe/pull/1413)

## [Master] - 2025-01-23

### Fixed

- MINOR There was a missing include (<deal.II/grid/cell_data.h>) for two tests. The include has been added. It has only been added for deal.II master since deal.II 9.6 does not have this include. [#1412](https://github.com/chaos-polymtl/lethe/pull/1412)

## [Master] - 2025-01-23

### Changed

- MINOR The clang-tidy CI now uses deal.II 9.6.0 instead of deal.II 9.5.1 [#1410](https://github.com/chaos-polymtl/lethe/pull/1410)

## [Master] - 2025-01-23

### Changed

- MINOR Update the installation instructions under WSL and Linux to use the deal.II 9.6.0 version instead of the 9.5.1 [#1409](https://github.com/chaos-polymtl/lethe/pull/1409)

## [Master] - 2025-01-20

### Changed

- MAJOR This change deprecates the "l2 smoothing factor" void fraction parameter and introduces the "l2 smoothing factor length". The previous parameters, named "l2 smoothing factor", meant the square of the desired smoothing length. The new parameter, named "l2 smoothing length", is the actual length of the smoothing region. This was done by adding an attribute to the void fraction class (l2_smoothing_factor) calculated as the parameter's square. This way, the meaning of the parameter are more understandable and more easily manipulated. [#1408](https://github.com/chaos-polymtl/lethe/pull/1408)

### Fixed

- MAJOR The change in #1406 introduced a major bug that would prevent the velocity, diameter, angular velocity and other particle properties from being adequately displayed in paraview. This PR fixes this bug by adequately positioning the output data when building the patches. [#1407](https://github.com/chaos-polymtl/lethe/pull/1407)

## [Master] - 2025-01-17

### Changed

- MAJOR The order of the particle properties used for DEM and CFD-DEM has been changed. The mass is now the third property (index 2) instead of being the last. This change, which is purely esthetic, allowed us to fix a leftover bug in the visualisation of DEM results where the array allocated for the particle properties in the VTU output was not the right size (it was one double too small). Surprisingly, this was not crashing, but it was wrong in any cases. Furthermore, we figured out that some of the restart CFD-DEM application tests (the particle sedimentation and the gas fluidized bed) had wrong restart files. This PR addresses all of those changes and refreshes the restart files due to the change in the order of the particle properties. [#1406](https://github.com/chaos-polymtl/lethe/pull/1406)

## [Master] - 2025-01-10

### Changed

- MAJOR The index of properties in the DEM and CFD-DEM simulations are now indexed using a template (PropertiesIndex) instead of using an hardcoded enum. This enables the different solvers to have different number of properties and different properties index. For large simulations, this gives minor performance gain for DEM simulations (e.g. around 5-10%), but this makes the solver significantly more flexible. This is a major change since it breaks DEM and CFD-DEM restart files from previous versions. [#1399](https://github.com/chaos-polymtl/lethe/pull/1399)

## [Master] - 2025-01-02

### Changed

- MINOR All lethe-particles tests that have a restart file now start with "restart_". This nomenclature change is made to ensure that restart tests are easier to identify. Furthermore, all restart file generators have been tested and missing generators have beem added. [#1402](https://github.com/chaos-polymtl/lethe/pull/1402)

## [Master] - 2025-01-01

### Changed

- MINOR All lethe-fluid-particles tests that have a restart file now start with "restart_". This nomenclature change is made to ensure that restart tests are easier to identify. Furthermore, all application_tests files now are seperated by underscores instead of a blend of hyphens and underscores. [#1400](https://github.com/chaos-polymtl/lethe/pull/1400)

## [Master] - 2024-12-17

### Added

- MAJOR A new mechanism is added that allows to launch a series of examples specified in contrib/validation/validation_cases.txt using the contrib/validation/validate_lethe.sh script. This script automatically launches the simulations that are specified in the validation cases and keeps the logs, the simulation results used to generate the plots, the plots and generates a pdf report with all of the main results of the validation cases. This will be used to monitor the stability of Lethe on more complicated test cases than those that are tested within the application_tests. [#1396](https://github.com/chaos-polymtl/lethe/pull/1396)

## [Master] - 2024-12-16

### Added

- MINOR A new post-processing code is now available for the rotating drum example. This code output the average velocity profile of particles perpendicular to the free surface. [#1394](https://github.com/chaos-polymtl/lethe/pull/1394)

- MINOR Time averaged temperature can now be attached to the vtus using a new parameter in the post-processing subsection. A new object AverageScalar has been created. It can be used to calculate the time average temperature from a certain point in time. [#1395](https://github.com/chaos-polymtl/lethe/pull/1395)

## [Master] - 2024-12-03

### Removed

- MAJOR The ability to bound the void fraction from below and above (using l2 lower bound and l2 upper bound parameters) has been removed. This bounding of the void fraction was highly problematic, since it could create discontinuities in the time derivative of the void fraction and contaminate the solution. In reality, this feature was never used. Importantly, this PR refactors the calculation of the void fraction outside of the VANS and CFD-DEM solvers of Lethe so that it exists only independently in a seperate subequation solver. This change is in preperation for the matrix free implementation of the VANS equations. [#1392](https://github.com/chaos-polymtl/lethe/pull/1392)

## Release of Lethe v1.0 - 2024-11-30

The lethe v1.0 release marks the transition of lethe to a numbered released format. Following numbered releases will be generated more often and we aim to release a numbered version 4 to 6 times per year or when sufficient changes have been realized.

The lethe v1.0 release comes with a structured and stable syntax. All of the examples and test provided with lethe were re-ran with this v1.0 release to ensure that all of the results provided online are up to date. This version is production-ready and has been tested and validated extensively.

## [Master] - 2024-11-28

### Changed

- MINOR The squared term that was added in the correction direction vector of the VOF DCDD stabilization [#1103](https://github.com/chaos-polymtl/lethe/pull/1103) has been removed after observing regression in computed results of VOF examples. [#1390](https://github.com/chaos-polymtl/lethe/pull/1390)

## [Master] - 2024-11-21

### Changed

- MINOR The default value for the particle weight for the load balancing was changed from 10K to 2K. [#1347](https://github.com/chaos-polymtl/lethe/pull/1347)

## [Master] - 2024-11-18

### Fixed

- MINOR Time step can now be changed after a restart if adaptive time stepping is disabled. [#1343](https://github.com/chaos-polymtl/lethe/pull/1343)

- MINOR Application of the immersed solid tanh diffusivity model is now faster thanks to using optimizations from the Shape class. [#1343](https://github.com/chaos-polymtl/lethe/pull/1343)

- MINOR The specific heat model in BDF2 was not protected against a division by zero (dH/dT). A tolerance is now added to the denominator (dH/(dT+tol)) to avoid this problem. [#1367](https://github.com/chaos-polymtl/lethe/pull/1367)
 
## [Master] - 2024-11-11

### Changed

- MAJOR The time step in the simulation control class is no longer modified by default to be exactly the end time of the simulation. Moreover, the time step is no longer modified to output Paraview files at certain times, therefore, the time output for transient simulations was refactored. [#1336](https://github.com/chaos-polymtl/lethe/pull/1341)
-
### Added

- MINOR Added bubble detachment in liquid shear flow example to the documentation [#1334](https://github.com/chaos-polymtl/lethe/pull/1334)

## [Master] - 2024-11-04

### Changed

- MINOR Curvature L2 projection within the VOF auxiliary physics now go through the VOF subequations interface. [#1336](https://github.com/chaos-polymtl/lethe/pull/1336)

### Added

- MINOR Added capillary rise example to the documentation [#1328](https://github.com/chaos-polymtl/lethe/pull/1328)

## [Master] - 2024-11-03

### Fix

- MINOR Force chains between local-ghost particle were being written multiple time when running in parallel. An arbitrary rule was added so that only one of the process is writing the force chain.[#1342](https://github.com/chaos-polymtl/lethe/pull/1342)


### Fix

- MINOR It was not possible to differentiate cohesive and repulsive forces with the force chains since the code was using the ".norm()" function. Now, the force chain calculation uses the scalar product between the normal force and the normal unit vector. [#1339](https://github.com/chaos-polymtl/lethe/pull/1339)

## [Master] - 2024-11-02

### Fix

- MINOR The checkpoint files generated from the lethe-particles solver were now incompatible with the CFD-DEM solvers (lethe-fluid-particles and lethe-fluid-vans) since the prefix are now ending with a "_0" or "_1". The "read_dem" function now reads the ".checkpoint_controller" file at first. [#1338](https://github.com/chaos-polymtl/lethe/pull/1338)

## [Master] - 2024-10-31

### Changed

- MAJOR The way Lethe manages boundary conditions was changed dramatically. As of this change, every boundary condition defined in the triangulation (the mesh), but have a corresponding boundary condition defined in the parameter file for every physics that is enabled. This will be tested by the solver at run time. For example, if a mesh containing four boundary conditions (id 0, 1, 2, 3) is used, then a boundary condition must be defined for all of these ids. In the past, a default boundary condition was automatically applied to these boundaries. As of this change, there is no default boundary condition. Although this change requires more work and slightly longer parameter files, it has enabled us to extensively refactor and simplify the way boundary conditions are managed. It also greatly enhances the sanity checks we can do with boundary conditions.

## [Master] - 2024-10-28

### Added

- MINOR The DEM solver generates two versions of checkpoint files alternately, so that at least one version is not corrupted when the simulation stops during a checkpoint procedure. [#1327](https://github.com/chaos-polymtl/lethe/pull/1327)

## [Master] - 2024-10-22

### Changed

- MAJOR Secondary equations (subequations) solved within the VOF auxiliary physics now go through a subequations interface similarly to how auxiliary physics go through the multiphysics interface. At the moment, only the implementation for the L2 projection of the phase fraction gradient has been refactored. Furthermore, only a linear equation solver has been implemented so far. [#1318](https://github.com/chaos-polymtl/lethe/pull/1318)

- MAJOR All scratch data objects now inherit from a base class, namely PhysicsScratchDataBase. In a similar manner, all base assemblers are now specialized types of the PhysicsAssemblerBase. [#1318](https://github.com/chaos-polymtl/lethe/pull/1318)
## [Master] - 2024-10-18

### Added

- MINOR Create compiler flag to use float precision for the geometric multigrid preconditioner in the matrix-free application. [#1319](https://github.com/chaos-polymtl/lethe/pull/1319)

## [Master] - 2024-10-16

### Changed

- MINOR Change default coarse grid solver to a direct solver for the geometric multigrid preconditioner in the matrix application. [#1322](https://github.com/chaos-polymtl/lethe/pull/1322)

### Fixed

- MINOR Added missing flag to update quadrature points for the FEFaceValues of the tracer physics. [#1323](https://github.com/chaos-polymtl/lethe/pull/1323)

### Added

- MAJOR Tracer physics can now be solved by using a discontinuous Galerkin method instead of a continuous Galerkin method. [#1320](https://github.com/chaos-polymtl/lethe/pull/1320)

## [Master] - 2024-10-15

### Added

- MINOR Five examples of the Cahn-Hilliard-Navier-Stokes solver were added to keep a trace of the work made with this solver. [#1307](https://github.com/chaos-polymtl/lethe/pull/1307)

## [Master] - 2024-10-11

### Fixed

- MAJOR Unexpected segmentation faults were occurring when using large simulation with a periodic boundary condition. The cause of this bug was identified and fixed. It was related to the resizing of the force, torque and displacement vectors that were not considering the ghost particles when resized. [#1316](https://github.com/chaos-polymtl/lethe/pull/1316)

## [Master] - 2024-10-09

### Changed

- MAJOR The boundary conditions for the fluid dynamics applications: matrix-based, matrix-free and block applications are now implemented only in the NavierStokesBase class. [#1313](https://github.com/chaos-polymtl/lethe/pull/1313)

## [Master] - 2024-10-04

### Changed

- MINOR Values outputted on the terminal have been uniformized using the log precision parameter. [#1310](https://github.com/chaos-polymtl/lethe/pull/1310)

## [Master] - 2024-10-03

### Fixed

- MINOR Simulations with time-dependent boundary conditions would not restart correctly in the matrix-free solver because the initial guess of the solution was reset to zero. This is fixed now, and the initial guess of the solution is identical. [#1302](https://github.com/chaos-polymtl/lethe/pull/1302)

### Changed

- MINOR The parameter "Lagrangian post-processing" has been renamed "lagrangian post-processing" to ensure that our parameter nomenclature (all lowercase) is homogenous. [#1303](https://github.com/chaos-polymtl/lethe/pull/1303)

## [Master] - 2024-09-30

### Fixed

- MINOR Simulations with time-averaging of the velocity fields were unable to restart when the domain was very large due to the fact that the restart vectors were read into the wrong vectors (locally_owned instead of locally_relevant). This PR fixes this. This also ensures that Lethe is able to restart with a different core-count than what was used to generate the restart file. [#1300](https://github.com/chaos-polymtl/lethe/pull/1300)

## [Master] - 2024-09-26

### Changed

- MINOR Made discontinuity-capturing directional dissipation (DCDD) stabilization optional for the VOF auxiliary physics. [#1296](https://github.com/chaos-polymtl/lethe/pull/1296)

## [Master] - 2024-09-25

## Changed

- MINOR The solutions for the mulitphysics are updated right after the non linear solve of the fluid dynamics matrix free solver. This ensures that the physics solved after the fluid dynamics have the most recent solution. [#1294](https://github.com/chaos-polymtl/lethe/pull/1294)

### Changed

- MINOR The tracer gradient used for DCDD stabilization calculations is now from the previous solution, which reduces the number of non-linear iterations. [#1293](https://github.com/chaos-polymtl/lethe/pull/1293)

## [Master] - 2024-09-25

### Changed

- MINOR Remove the GIFs from the repository to reduce its size. [#1289](https://github.com/chaos-polymtl/lethe/pull/1289)

### Added

- MINOR Added verbosity for the interface thickness in the CHNS solver. An application test was added to verify the feature resists to future changes made in the code. The documentation was updated to display the new feature in the Cahn-Hilliard section of the documentation. [#1291](https://github.com/chaos-polymtl/lethe/pull/1291)

## [Master] - 2024-09-24

### Added

- MINOR The outlet boundary condition was implemented in the matrix-free operators. A test in 2D with a cylinder close to the outlet was added to test this boundary condition for both the matrix-based and the matrix-free application. [#1287](https://github.com/chaos-polymtl/lethe/pull/1287)

## [Master] - 2024-09-23

### Fixed

- MINOR In the [PR #994](https://github.com/chaos-polymtl/lethe/pull/994), BDF extrapolation of the velocity in the VOF auxiliary physics was moved to the scratch data. However, the implementation was bypassed by an if condition. This bug is fixed with this PR, and an application test was added to ensure that the feature remains intact with future implementations. [#1286](https://github.com/chaos-polymtl/lethe/pull/1286)

## [Master] - 2024-09-20

### Added

- MINOR Added dimensionality for the mobility and interface thickness in the CHNS solver, and surface tension related parameters in the VOF solver. [#1274](https://github.com/chaos-polymtl/lethe/pull/1274)

## [Master] - 2024-09-12

### Fixed

- MINOR The previous implementation to calculate cell diameter used the `measure` function which assumed a linear mapping. A function that calculates the volume by summing JxW values returned by the FEValues object now replaces the `measure` function. [#1279](https://github.com/chaos-polymtl/lethe/pull/1279)

## [Master] - 2024-09-10

### Fixed

- MINOR A bug in the tracer physics where non-zero constraints were not reinitialized with the locally relevant DOFs made it that processes ran out of memory in large scale simulations. This is fixed. [#1278](https://github.com/chaos-polymtl/lethe/pull/1278)

## [Master] - 2024-09-09

### Added

- MINOR Added missing documentation entries regarding mesh parameters for the DEM solver and geometrical manipulations at initialization. [#1277](https://github.com/chaos-polymtl/lethe/pull/1277)

## [Master] - 2024-08-29

### Changed

- MINOR Rendered discontinuity-capturing directional dissipation optional for the heat transfer. [#1268](https://github.com/chaos-polymtl/lethe/pull/1268)

## [Master] - 2024-08-27

### Changed

- MINOR Refactored cell diameter computation as an inline function in the utilities. [#1265](https://github.com/chaos-polymtl/lethe/pull/1265)

## [Master] - 2024-08-23

### Added

- MINOR Added the capacity to print the parameter from the parameter file that were changed compared to the default parameters as well as the lethe and deal.II commit hash. This is achieved by adding the runtime argument --print-parameters to the command line arguments. [#1255](https://github.com/chaos-polymtl/lethe/pull/1255) and [#1257](https://github.com/chaos-polymtl/lethe/pull/1257)

## [Master] - 2024-08-09

### Changed

- MINOR Renamed the main file related to the `lethe-fluid` application: `gls_navier_stokes` to `fluid_dynamics_matrix_based`. The name of the `GLSNavierStokesSolver` class was changed to `FluidDynamicsMatrixBased`. [#1236](https://github.com/chaos-polymtl/lethe/pull/1236)

## [Master] - 2024-08-07

### Changed

- MINOR Renamed the main file related to the `lethe-fluid-sharp` application: `gls_sharp_navier_stokes` to `fluid_dynamics_sharp`. The name of the `GLSSharpNavierStokesSolver` class was changed to `FluidDynamicsSharp`. [#1231](https://github.com/chaos-polymtl/lethe/pull/1231)

### Changed

- MINOR Renamed the main file related to the `lethe-fluid-nitsche` application: `gls_nitsche_navier_stokes` to `fluid_dynamics_nitsche`. The name of the `GLSNitscheNavierStokesSolver` class was changed to `FluidDynamicsNitsche`. [#1228](https://github.com/chaos-polymtl/lethe/pull/1228)

## [Master] - 2024-08-06

### Changed

- MINOR Renamed the main file related to the `lethe-fluid-block` application: `gd_navier_stokes` to `fluid_dynamics_block`. The name of the `GDNavierStokesAssembler...` classes were changed to `BlockNavierStokesAssembler...`, and the main class `GDNavierStokesSolver` was also renamed to `FluidDynamicsBlock`. [#1226](https://github.com/chaos-polymtl/lethe/pull/1226)

## [Master] - 2024-08-05

### Changed

- MINOR Renamed the main file related to the `lethe-fluid-vans` application: `gls_vans` to `fluid_dynamics_vans`. The name of the `GLSVansAssembler...` classes were changed to `VANSAssembler...`, and the main class `GLSVANSSolver` was also renamed to `FluidDynamicsVANS`. [#1225](https://github.com/chaos-polymtl/lethe/pull/1225)

### Changed

- MINOR Renamed all the files related to the `lethe-fluid-matrix-free` application: `mf_navier_stokes` to `fluid_dynamics_matrix_free` and `mf_navier_stokes_operators` to `fluid_dynamics_matrix_free_operators`. The main class `MFNavierStokesSolver` was also renamed to `FluidDynamicsMatrixFree`. [#1222](https://github.com/chaos-polymtl/lethe/pull/1222)

### Added

- MINOR A post-processing feature for the tracer flow rate through boundaries is added. This allows to produce data useful for residence time distribution analyses. [#1197](https://github.com/chaos-polymtl/lethe/pull/1197)

### Fixed

- MINOR Conditions for whether a boundary condition should be updated were made less restrictive. Now, each physics is responsible for updating or not their boundary conditions. [#1197](https://github.com/chaos-polymtl/lethe/pull/1197)

### Removed

- MINOR Removed the gear3 DEM integrator. It was never used, it did not possess any unit test and it was not even clear if it still worked. [#1221](https://github.com/chaos-polymtl/lethe/pull/1211)

## [Master] - 2024-08-01

### Added

- MINOR The possibility to specify an intermediate level as coarse grid solver for the gcmg preconditioner was added as a parameter. It allows to perform several v-cycles at the level chosen. [#1211](https://github.com/chaos-polymtl/lethe/pull/1211)

## [Master] - 2024-07-31

### Fixed

- MINOR The Coulomb's criterion was wrong in the particle-particle contact force for Hertz-Mindlin with limit force, Hertz and Linear in DEM. The normal force norm was explicitly positive when doing normal_force.norm(), making the Coulomb's criterion always positive even if the particles are in repulsion, so in slidling. This has been fixed using the same method as for Hertz-Mindlin with limit overlap is calculated. [#1216](https://github.com/chaos-polymtl/lethe/pull/1216)

### Added

- MINOR P-multigrid was added to the gcmg preconditioner of the lethe-fluid-matrix-free application. It supports three different coarsening strategies to define the degree p of the different levels. It also allows to use hybrid hp- and ph-multigrid strategies. [#1209](https://github.com/chaos-polymtl/lethe/pull/1209)

## [Master] - 2024-07-24

### Changed

- MINOR Forces a contact search at the last DEM iteration of a CFD iteration for more robustness related to the update of the reference location of the particles prior the void fraction calculation [#1205](https://github.com/chaos-polymtl/lethe/pull/1205)

## [Master] - 2024-07-23

### Changed

- MINOR PBC with QCM in CFD-DEM are working without forcing the particle displacement at the last DEM iteration of the CFD iteration. It is removed. [#1204](https://github.com/chaos-polymtl/lethe/pull/1204)

## [Master] - 2024-07-21

### Changed

- MINOR The load balancing functions of the DEM and coupling CFD-DEM solvers are encapsulated in a new class, removing the duplicated code. [#1199](https://github.com/chaos-polymtl/lethe/pull/1199)

### Fixed

- MINOR The load balancing is fixed and working for the coupling CFD-DEM solver. It was not working after a small refactoring. An application test is also added. [#1199](https://github.com/chaos-polymtl/lethe/pull/1199)

## [Master] - 2024-07-20

### Fixed

- MINOR The ratio of the critical Rayleigh time step was wrong in CFD-DEM and was modified as done in DEM. [#1203](https://github.com/chaos-polymtl/lethe/pull/1203)

## [Master] - 2024-07-14

### Added

- MINOR Following the implementation of temperature-dependent solid domain constraint for both single-phase and two-phase flows, a need for an additional condition was identified. Indeed, in some simulations, splatters may happen, and with the previous implementation, splatter could freeze "mid-air" which is non-physical and therefore undesired. A plane to divide the domain in two parts and only the one in the opposite direction to the normal vector is considered for stasis constraints to be applied. The plane is defined through a point and an outward pointing normal vector. [#1193](https://github.com/chaos-polymtl/lethe/pull/1193)

### Fixed

- MINOR Previously, the pressure DOFs were also constrained with the "constrain stasis" feature, but this resulted in an ill-posed problem when sources terms in the momentum equation were pressure-dependent. Therefore, the constraints on pressure DOFs are now removed. [#1193](https://github.com/chaos-polymtl/lethe/pull/1193)

## [Master] - 2024-07-12

### Added

- MINOR A tracer diffusivity model is added, to be used with the `lethe-fluid-sharp` solver. It makes it possible to assign distance-based (depth) diffusivity to the immersed solids. [#1185](https://github.com/chaos-polymtl/lethe/pull/1185)

## [Master] - 2024-07-09

### Added

- MINOR The multigrid output now also prints the workload imbalance and vertical communication efficiency of the multigrid hierarchy being used. In addition, the mulrigrid timers now print the min max and average times correctly with the appropriate rank. [#1194](https://github.com/chaos-polymtl/lethe/pull/1194)

### Fixed

- MINOR The lethe-fluid-nitsche solver could not checkpoint once it had been restarted. This has been fixed by ensuring that the particles are properly checkpointed and read from the restart files. [#1192](https://github.com/chaos-polymtl/lethe/pull/1192)

## [Master] - 2024-07-08

### Changed

- MINOR  The timer output for the lethe-fluid and lethe-fluid-matrix-free applications contains now more detailed elements. Including post-processing capabilities that were not timed, the time spent in reading mesh and manifolds, and setting the initial conditions. [#1187](https://github.com/chaos-polymtl/lethe/pull/1187) 

## [Master] - 2024-07-05

### Changed

- MINOR  The NS-VOF assembler was not optimized as the NS assembler. The if conditions on the dof components were added to not compute multiplications that are a priori 0 in the matrix and rhs assembly. [#1188](https://github.com/chaos-polymtl/lethe/pull/1188)

## [Master] - 2024-07-04

### Changed

- MINOR  The `entry_string_to_tensor` functions were duplicated; they were in `utilities.h` and `parameters.h`. In one case, it was templated with `spacedim` and the other not. This PR merges both functions in `utilities.h` as `value_string_to_tensor`. Since these are small functions, they are now inline functions. [#1186](https://github.com/chaos-polymtl/lethe/pull/1186)

## [Master] - 2024-07-03

### Fixed

- MINOR The term associated with the viscosity jump in the strong residual for the NS-VOF assembler added in [#1049] (https://github.com/chaos-polymtl/lethe/pull/1149) led to an ill-posed formulation in solid-liquid phase change cases where the solid is represented by a highly viscous fluid. Hence, this term was removed from the current formulation. [#1180](https://github.com/chaos-polymtl/lethe/pull/1180)

### Changed

- MINOR LPBF benchmark example prm files are up to date with the current working version. [#1180](https://github.com/chaos-polymtl/lethe/pull/1180)

### Added

- MINOR There was no appropriate application test for lpbf/phase change with vof cases that could have allowed us to identify this ill-posed formulation prior to the merge of [#1049] (https://github.com/chaos-polymtl/lethe/pull/1149). A new application test is now included to avoid future mistake in lpbf/phase change with vof cases. [#1180](https://github.com/chaos-polymtl/lethe/pull/1180)

## [Master] - 2024-07-02

### Fixed

- MINOR Box refinement was applied when restarting simulations, changing the restarting triangulation and leading to undesirable behavior. A boolean stating if the simulation restarts is now checked before applying the box refinement. [#1184](https://github.com/chaos-polymtl/lethe/pull/1184) 

## [Master] - 2024-06-16

### Changed

- MINOR Some application-test restart files have been updated using p4est 2.3.6. Some test results have changed for lethe-fluid-particles and lethe-fluid-vans, since the DEM solver has slightly changed since the previous restart files generation, and it is now impossible to regenerate the exact same initial condition. [#1181](https://github.com/chaos-polymtl/lethe/pull/1181)


## [Master] - 2024-06-16

### Changed

- MINOR A new subsection ``solid surfaces`` needs to be used when defining solid objects in the parameter file. For more information, see the solid object documentation. [#1169](https://github.com/chaos-polymtl/lethe/pull/1169)
- MINOR The center of rotation of a solid object is no longer being defined using a subsection and three parameters. It is now defined with one parameter and a list of doubles. For more information, see the solid object documentation. [#1169](https://github.com/chaos-polymtl/lethe/pull/1169)

## [Master] - 2024-06-13

### Added

- MINOR The DEM solver supports deprecated parameters when 3 individual component parameters are changed to a list of values parameter. [#1171](https://github.com/chaos-polymtl/lethe/pull/1171)

### Changed

- MINOR The parameters <gx>, <gy> and <gz> of the DEM solver are changed to the parameter <g>. [#1171](https://github.com/chaos-polymtl/lethe/pull/1171)

## [Master] - 2024-06-05

### Added

- MINOR A cylindrical manifold type has been added in the manifold parameter file option. Input arguments arg1, arg2, arg3, etc. used to describe the manifold geometry were modified to point coordinates and direction vector. A manifold section has been added to the documentation [#1167](https://github.com/chaos-polymtl/lethe/pull/1167)

## [Master] - 2024-05-27

### Changed

- MINOR The "insertion file name" parameter has been renamed to "list of input files". It now supports a list of input files that will be used in order when inserting particles. [#1135](https://github.com/chaos-polymtl/lethe/pull/1135)

## [Master] - 2024-05-26

### Fixed

- MINOR Simulation restart would not work adequately when a boundary condition was time-dependent and this would lead unstable simulation restarts. [#1158](https://github.com/chaos-polymtl/lethe/pull/1158) 

## [Master] - 2024-05-23

### Added

- MINOR dependency on muParser. This requires to load an additional module "muparser/2.3.2" when compiling deal.II and Lethe on clusters. [#1143](https://github.com/chaos-polymtl/lethe/pull/1143) 

## [Master] - 2024-05-20

### Removed

- MINOR The initial conditions executable which could be used to test some initial conditions has now been removed completely. This executable was never used and was unstable to the multiphysics interface. [#1148](https://github.com/chaos-polymtl/lethe/pull/1148)

## [Master] - 2024-05-13

### Changed

- MAJOR The source term specifications for the Navier-Stokes equations were using the label "xyz" which made no sense. This has been changed to "fluid dynamics" [#1130](https://github.com/chaos-polymtl/lethe/pull/1130)

## [Master] - 2024-05-09

### Added

- MAJOR Added advection of particles in the Adaptive Sparse Contacts (ASC) feature in DEM and CFD-DEM (Disabling Contacts is renamed Adaptive Sparse Contacts). [#1113](https://github.com/chaos-polymtl/lethe/pull/1113)

## [Master] - 2024-05-09

### Added

- MINOR Added a condition to the generation of .vtu and .pvd files in addition to the output generation frequency. If specified in the parameter file, only the results within the specified simulation time interval will be generated. [#1120] (https://github.com/chaos-polymtl/lethe/pull/1120) 


## [Master] - 2024-05-02

### Fixed

- MAJOR Periodic boundary conditions for the DEM were not working if they weren't the last boundary conditions being declare in the parameter file. Now, every boundary condition work in which even order they are being declared. [#1110](https://github.com/chaos-polymtl/lethe/pull/1110)

## [Master] - 2024-04-30

### Fixed

- MINOR The lethe-fluid-nitsche solver was unable to restart when the immersed triangulation was made of simplices. This has been fixed, however, mesh adaptation crashes when it is done after the restart process. I (problembär) have an idea why (the previous particle_handler of the previous checkpoint is still registered somehow in the triangulation), but I will need more time to come up with an adequate solution. [#1106](https://github.com/chaos-polymtl/lethe/pull/1106)

## [Master] - 2024-04-24

### Added

- MAJOR A new parameter to reuse the preconditioner between Newton iterations was added. Therefore, the logic of the Newton solver was also modified as it required a different call for the assembly of the matrix and the set up of the preconditioner. These calls were separated for all relevant physics solvers. [#1102](https://github.com/chaos-polymtl/lethe/pull/1102)

## [Master] - 2024-04-24

### Changed

- MAJOR For the multigrid preconditioners of the matrix-free application a new class was added, where the constructor sets the operators, constraints and transfers, and the initialize function sets the smoother, coarse-grid solver and final multigrid object. This allows to reuse the operators, constraints and transfers whenever possible reducing computational cost. [#1102](https://github.com/chaos-polymtl/lethe/pull/1102)

## [Master] - 2024-04-18

### Fixed

- MINOR Checkpoints were not established correctly with the lethe-fluid-particles solver. The void fraction was delayed in time and the restart file would not yield the same result because of this, resulting in crashes when the time derivative of the void fraction was used in the calculations. [#1096](https://github.com/chaos-polymtl/lethe/pull/1096)

## [Master] - 2024-04-17

### Added

- MINOR Added Cahn-Hilliard equations energy computation and output. The documentation was updated accordingly. [#1095](https://github.com/chaos-polymtl/lethe/pull/1095)

## [Master] - 2024-04-16

### Fixed

- MINOR Simulations with adaptive time-stepping would checkpoint correctly, but would not restart with the correct time-step due to the way the time-step was read from the checkpointing files. This is now fixed and tested. [#1093](https://github.com/chaos-polymtl/lethe/pull/1093)

## [Master] - 2024-04-15

### Fixed

- MINOR The previous flow rate post processing tool would add a line to the .dat file for each defined boundary (in the boundary condition subsection of the .prm file) for every time iteration. This is not the expected output as one would expect a single line to be added for each time step with the flow rates of every boundary. This is now fixed. [#1092](https://github.com/chaos-polymtl/lethe/pull/1092)

## [Master] - 2024-04-07

### Added

- MINOR The lethe-fluid-matrix-free solver would not work adequately on mesh that were not refined homogenously if all the boundary conditions of the domain were periodic or dirichlet boundary condition. This is fixed by adding a new parameter to fix the pressure constant in the coarse grid solver. This prevents the coarse grid linear solver from failing. [#1089](https://github.com/chaos-polymtl/lethe/pull/1089)

## [Master] - 2024-04-05

### Changed

- MINOR The output at every iteration in the terminal of the non-linear solver now displays the norms of the correction for both the velocity and pressure contributions for the Navier-Stokes solver when `verbosity` in the `non-linear solver` subsection is set to `verbose`. [#1076](https://github.com/chaos-polymtl/lethe/pull/1076)

## [Master] - 2024-04-04

### Added

- MINOR Added the turbulent Taylor-Couette incompressible flow example. [#1083](https://github.com/chaos-polymtl/lethe/pull/1083)

## [Master] - 2024-03-28

### Changed

- MINOR The ``velocity x/y/z`` parameters have been replaced by one parameter named ``initial velocity``. The ``omega x/y/z`` parameters have been replaced by one parameter named ``initial angular velocity``. [#1078](https://github.com/chaos-polymtl/lethe/pull/1078)
 [#1078](https://github.com/chaos-polymtl/lethe/pull/1078)


## [Master] - 2024-03-20

### Changed

- MINOR The ``insertion box minimum/maximum x/y/z`` parameters have been replaced by one parameter named ``insertion box points coordinates``. The ``insertion first/second/third direction `` parameters have been replaced by one parameter named ``insertion direction sequence``. [#1074](https://github.com/chaos-polymtl/lethe/pull/1074)

## [Master] - 2024-03-19

### Added

- MINOR The mass conservation feature of the VOF model now also monitors the momentum of both phases. This should be useful for some cases and improve the sharpening mechanism. [#1073](https://github.com/chaos-polymtl/lethe/pull/1073)

## [Master] - 2024-03-16

### Fixed

- MINOR The recently implemented VOF stasis constraint ([#1048](https://github.com/chaos-polymtl/lethe/pull/1048)) was not compatible with adaptive mesh refinements and checkpointing system. When refining the mesh or when restarting a simulation, the filtered phase fraction was reinitialized, but not filled with the appropriate values. This issue is now solved by adding `apply_phase_filter()` to `post_mesh_adaptation()` and `read_checkpoint()` VOF functions. [#1072](https://github.com/chaos-polymtl/lethe/pull/1072)

## [Master] - 2024-03-13

### Changed

- MINOR In a similar way to "VOF Barycenter", the output of the "VOF Mass Conservation" is now displayed at every time iteration on the terminal when `calculate mass conservation` in the `post-processing` is enabled (default) and the `verbosity` is set to `verbose`. [#1062](https://github.com/chaos-polymtl/lethe/pull/1062)

- MINOR A new boolean argument in `make_table` functions has been added, namely `display_scientific_notation` to decide if the display of table values is displayed at a fixed decimal or in scientific notation. [#1062](https://github.com/chaos-polymtl/lethe/pull/1062)

### Added

- MINOR Assertion on fields required for physical property models were added to the code. This is done to ensure that if a physical property model necessitating a certain field is employed and that field is not defined, the exception message thrown is comprehensible. [#1065](https://github.com/chaos-polymtl/lethe/pull/1065)

- MINOR A new `make_table_scalars_vectors` function has been added in `utilities.h`(`.cc`). [#1062](https://github.com/chaos-polymtl/lethe/pull/1062)

## [Master] - 2024-03-11

### Changed

- MINOR The eigenvalue estimation for the multigrid preconditioners in the lethe-fluid-matrix-free solver is now done internally by the PreconditionRelaxation class instead of using a PreconditionChebyshev and an estimate omega function [#1064](https://github.com/chaos-polymtl/lethe/pull/1064).

## [Master] - 2024-03-11

### Added

- MINOR Temperature-dependent `stasis constraint` is now featured in the Melting Cavity heat transfer example. [#1061](https://github.com/chaos-polymtl/lethe/pull/1061)

## [Master] - 2024-03-05

### Added

- MINOR The new temperature-dependent solid domain constraints for two-phase Volume of Fluid (VOF) simulations has been implemented (`constrain_stasis_with_temperature_vof()`). To select cells on which the constraints are applied, we check if the filtered phase fraction at quadrature points in the cell are within a range of accepted values. This range of accepted values is defined with the `phase fraction tolerance` parameter. This parameter defines an absolute tolerance on filtered phase fraction. Related documentation was updated and a new application test was also added. [#1048](https://github.com/chaos-polymtl/lethe/pull/1048)

- MINOR Added the capability to locally refine the mesh in the vicinity of the boundary conditions as specified by a list of boundary conditions. [#1056](https://github.com/chaos-polymtl/lethe/pull/1056)

### Changed

- MINOR Previously implemented temperature-dependant solid domain constraints for one-fluid simulations (`constrain_solid_domain()`) was renamed to `constrain_stasis_with_temperature()` and its content changed a bit to avoid copied lines in `constrain_stasis_with_temperature_vof()`. [#1048](https://github.com/chaos-polymtl/lethe/pull/1048)

- MINOR The parameter subsection `constrain solid domain` was renamed to `constrain stasis`. Related documentation was also updated. [#1048](https://github.com/chaos-polymtl/lethe/pull/1048)

### Added - 2024-03-04

- MINOR The "file" insertion method has been added to the DEM and CFD-DEM solvers.[#1054](https://github.com/chaos-polymtl/lethe/pull/1054)

### Fixed

- MINOR In LPBF simulations, recoil pressure formulation $`p_\mathrm{rec} = 0.55p_\mathrm{sat} + [1/\rho]m_\mathrm{dot}^2`$ was accounting 2 times for the term $`[1/\rho]m_\mathrm{dot}^2`$ because the latter is included in the term $`0.55p_\mathrm{sat}`$. The formulation is now corrected and reads $`p_\mathrm{rec} = 0.55p_\mathrm{sat}`$.

## [Master] - 2024-03-02

### Fixed

- MINOR The previously implemented weak compressibility ([#815](https://github.com/chaos-polymtl/lethe/pull/815)) was not working with the heat transfer (HT) physic. The missing pressure field was added to scratch data. Two application tests were added to test the compressibility with HT for single-fluid and VOF simulations. [#1051](https://github.com/chaos-polymtl/lethe/pull/1051)

## [Master] - 2024-02-27

### Added

- MINOR Temperature-dependant solid domain constraints were added for one-fluid simulations. In order to improve calculation time and linear system conditioning, homogeneous constraints are imposed on velocity DOFs in a defined solid domain through a temperature range. Additionally, pressure DOFs in the same solid domain that are not connected to fluid cells are also imposed with homogeneous constraints. The temperature range is specified through a new parameter subsection, namely `constrain solid domain`. The use of a `dynamic_zero_constraints` AffineConstraints object, rather than reconstruction of the `zero_constraints` one at every time step was determined to be optimal after some profiling tests. Documentation and a new application test have been also added. [#1038](https://github.com/chaos-polymtl/lethe/pull/1038)

### Changed

- MINOR The `Variable` enum class was moved from `parameters.h` to `multiphysics.h`. Using variables from the `Variable` enum class, the constrained and constraining fields could be parametrized in the future according to user needs. [#1038](https://github.com/chaos-polymtl/lethe/pull/1038)

### Fixed

- MINOR In `establish_solid_domain()` the pressure was not constrained in solid cells that were not connected to a fluid due to a small logic mistake. This was corrected. [#1038](https://github.com/chaos-polymtl/lethe/pull/1038)

## [Master] - 2024-02-26

### Fixed

- MAJOR The first level of the fine search candidates for particle-particle, particle-wall and particle-floating walls was never cleared even when the value of the unordered_map became an empty unordered_map. In simulations with load balancing or particles that were moving significantly, this could potentially lead to a scenario where the size of the first level of the fine search candidate became equal to the number of total particles in the simulation, potentially leading to a crash.

## [Master] - 2024-02-09

### Fixed

- MAJOR Restarts using lethe-fluid-particles with "void fraction time derivative = true" would be incoherent because the void fraction was reinitialized when restarting, leading to a wrong time derivative of the void fraction. This has been patched.

## [Master] - 2024-01-31

### Changed

- MAJOR Robin boundary condition was renamed to incorporate the imposed heat flux (Neumann). The corresponding type parameter was ``convection-radiation`` and the current is ``convection-radiation-flux``. All application tests were updated accordingly.

## [Master] - 2024-01-28

### Fixed

- MAJOR All application_tests that use files have now mpirun=1 in their output to ensure that they are run from a mpirun=1 folder. Since https://github.com/dealii/dealii/pull/16551 deal.II 9.6 uses a serial folder for serial tests. However, deal.II 9.5 does not. To maintain portability between the two versions, we manually force a folder called mpirun=1 when running tests and adapt the copy of the files to this folder. This is a temporary fix that can be reverted once we drop support for deal.II 9.5

## [Master] - 2024-01-24

### Added

- MINOR Darcy-like penalization is now extended to two-fluid VOF simulations with the `PhaseChangeDarcyVOFAssembler` assembler class. Related documentation has been updated and a new application test was added. [#990](https://github.com/chaos-polymtl/lethe/pull/990)

- MINOR Extrapolation in time of temperature and temperature gradient values has been added. This is used for source terms in the Navier-Stokes equations requiring them (Marangoni and evaporation recoil force terms). [#994](https://github.com/chaos-polymtl/lethe/pull/994)

### Changed

- MINOR Velocity extrapolation has been moved to `include/solvers/vof_scratch_data.h`. Solution extrapolations are now located in `ScratchData` to avoid constraint on assemblers order. [#994](https://github.com/chaos-polymtl/lethe/pull/994)

## [Master] - 2024-01-24

### Fixed

- MINOR In the previous implementation of the heat flux postprocessor ([#953](https://github.com/chaos-polymtl/lethe/pull/953)), there was a confusion related to the `material_id` of a cell. As a cell's `material_id` is defined in the mesh, fluids moving from a cell to another did not have the right `material_id` when outputted. Furthermore, some deal.II meshes have an intrinsic way of numbering `material_ids` (e.g. colorized subdivided_hyper_rectangle) that were leading to inadequate comparison of the `material_ids`. As of now, in the heat flux postprocessor, the `material_id` is only used to indentify cell material when solids are present. When simulating with fluids and solids, cells of the mesh in which the fluids lie should have a `material_id` of `0`. For _fluids-only_ simulations, their heat fluxes are outputted on the entire domain, as the current `DataPostprocessor` can only hold 1 `dof_handler`. During postprocessing, the user could clip the domain of interest of each fluid using the `phase` or `filtered_phase` fields outputted. [#988](https://github.com/chaos-polymtl/lethe/pull/988)

## [Master] - 2023-12-27

### Changed

- MINOR  Phase change thermal expansion now considers the solid thermal expansion within the mushy zone (between Tsolidus and Tliquidus) instead of using a blending rule. This leads to more stable results in the melting cavity benchmark

## [Master] - 2023-12-22

### Fixed

- MINOR  PVD handles for solid DEM objects were written with the wrong text format ("_" separators instead of ".")

### [Master] - Changed - 2023-12-17

- MINOR The "insertion random number range" and "insertion random number seed" parameters got renamed to "insertion maximum offset" and "insertion prn seed" respectively. The old names didn't make sense, as they're not random because they're defined in the parameter file. [#970](https://github.com/chaos-polymtl/lethe/pull/896)

## [Master] - 2023-12-11

### Fixed

- MINOR The number of remaining particle to insert of each type is being checkpointed adequately. This means that no modification are required to the "number of particle" parameter after restarting a simulation. [#964](https://github.com/chaos-polymtl/lethe/pull/964)

### Fixed

- MINOR Solid objects can now be restarted adequately in DEM. They will resume at the position they had at the end of the simulation [#959](https://github.com/chaos-polymtl/lethe/pull/959)

## [Master] - 2023-11-27

### Added

- MINOR A few tables were omitted from being "checkpointed" causing a lost of valuable information when simulations were restarted (Issue: [#916](https://github.com/chaos-polymtl/lethe/issues/916)). The missing tables have been added to the `write_checkpoint()` and `read_checkpoint()` of each physic. [#938](https://github.com/chaos-polymtl/lethe/pull/938)

## [Master] - 2023-11-27

### Removed

- MINOR The average diameter of the uniform size distribution with the DEM module was specified with an "average diameter" parameter. It is now specified directly from the diameter parameter. This is correctly documented. [#940](https://github.com/chaos-polymtl/lethe/pull/940).

### Fixed

- MINOR The DEM time step verification was outputting the most permissive time step (the biggest) and not the most restrictive (the smallest). This bugfix doesn't affect the uniform particle size simulation. [#939](https://github.com/chaos-polymtl/lethe/pull/939).


## [Master] - 2023-11-23

### Fixed

- MINOR The plane insertion for the DEM was only supporting the uniform diameter distribution. Now it supports all types of distribution.

## [Master] - 2023-11-16

### Changed

- MINOR The maximum number of boundary conditions for all physics was fixed to 14 since the boundary conditions had to be declared before being parsed. A new mechanism is now in place which parses the "number" parameter for each physics and keeps the maximal value. Then, this maximal value is used to pre-declare the boundary conditions. This enables much more robust sanity checking of the input parameter. The major drawback of this (and this is a major one) is that if we ever have another parameter with the name "number" then this parameter would also be parsed and used to establish the maximum number of boundary conditions. In this case, the best approach would be to replace "number" with "number of boundary conditions" in the parameter file. I (B-saurus-rex) did not want to do this at the time of this change to not have a massive PR which breaks every parameter files.

- MAJOR The "number" parameter within "subsection lagrangian physical properties" and "particle type n" was changed to "number of particles" to prevent confusions with the "number" used for boundary conditions. The "number" for boundary conditions will be changed to "number of boundary conditions" in the near future.

## [Master] - 2023-11-12

### Deprecated

- MINOR The uniform insertion method had been removed. The non-uniform insertion method has been renamed to volume method to remain coherent with the plane method. If you want to use an insertion method equivalent to the uniform insertion method, use the volume method with an "insertion random number range " equal to zero. [#926](https://github.com/chaos-polymtl/lethe/pull/926)


## [Master] - 2023-10-01

### Fixed

- MINOR The calculation of the velocity on rotating walls was calculated inadequately [#896](https://github.com/chaos-polymtl/lethe/pull/896): The velocity of rotating boundary conditions (e.g. boundaries of the mesh) was calculated inadequately in the case where the normal vector of the wall was aligned with the rotation axis. The whole calculation procedure was slightly messed up and only worked for cylinders. This has been fixed and made general, paving the way for a full refactor of the calculation of the particle-wall contact force.

## [Master] - 2023-09-30

### Fixed

- MINOR Affine constraints used to bound the void fraction and used as boundary conditions within the heat transfer and the block Navier-stokes solver were corrected [#885](https://github.com/chaos-polymtl/lethe/pull/895): There was an error in the application of the affine constraints used to clip the void fraction in the lethe-fluid-vans and lethe-fluid-particles solvers. This led to assertions being thrown in debug. This has been corrected by reiniting the constraints with the appropriate size. For the heat transfer and the block Navier-stokes, the issue was that the constraints were never reinited with the correct size, so they contained ghosted elements. This was caught by a new assert introduced in 2023-09 within deal.II master.

## [Master] - 2023-09-20

### Deprecated

- MINOR The calculation of the source term was enabled using a parameter called "enable". This parameter was used in some physics, not in others and was poorly implemented. We deprecate the usage of this parameter and always enable source term, considering the fact that the default value of the source term is zero anyway. This prevents the false perception that source terms could be enabled or disabled, while the behavior was inconsistent across physics.


## [Master] - 2023-09-19

### Changed

- MAJOR All the applications were renamed [#882](https://github.com/chaos-polymtl/lethe/pull/882): `gls_navier_stokes` is now `lethe-fluid`, `gd_navier_stokes` is now `lethe-fluid-block`, `nitsche_navier_stokes` is now `lethe-fluid-nitsche`, `gls_sharp_navier_stokes` is now `lethe-fluid-sharp`, `gls_vans` is now `lethe-fluid-vans`, `mf_navier_stokes` is now `lethe-fluid-matrix-free`, `dem` is now `lethe-particles`, `cfd_dem_coupling` is now `lethe-fluid-particles`, `rpt3d` is now `lethe-rpt-3d`, `rpt_cell_reconstruction_3d` is now `lethe-rpt-cell-reconstruction-3d`, `rpt_fem_reconstruction_3d` is now `lethe-rpt-fem-reconstruction-3d`, and `rpt_l2_projection_3d` is now `lethe-rpt-l2-projection-3d`.


## [Master] - 2023-10-30

- MINOR The rotational vector for the rotational boundary condition in the lethe-particles solver is now define with one line in the parameters file. [#920](https://github.com/chaos-polymtl/lethe/pull/920)



## [Sample] - YYYY/MM/DD

### Added

- MAJOR/MINOR/PATCH Description (#PR).

### Changed

- MAJOR/MINOR/PATCH Description (#PR).

### Deprecated

- MAJOR/MINOR/PATCH Description (#PR).

### Removed

- MAJOR/MINOR/PATCH Description (#PR).

### Fixed

- MAJOR/MINOR/PATCH Description (#PR).
