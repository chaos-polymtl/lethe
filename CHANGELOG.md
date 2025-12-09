
# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

### [Master] - 2025-12-09

### Added

- MAJOR The user can now generate particle using a ``lognormal`` distribution. Two new parameter have been added, namely the `minimum diameter cutoff ` and the `maximum diameter cutoff` which are used to specify the minimal and maximal particle diameter that can be inserted in a DEM simulation. By default, for a ``normal`` the smallest and biggest diameter that can be generated are 2.5 standard deviation from the means, which wasn't the case previously. [#1837](https://github.com/chaos-polymtl/lethe/pull/1837)

### [Master] - 2025-12-08

### Fixed

- MINOR The SerialSolid class, which is in the core library, had become dependent on the DEM library which is not something that fits our architecture. This was because of a single functor used for comparisons. This has been fixed by moving these functors used for cell comparisons to utilities.h. In the future we might coalesce these functions into a separate header file. [#1842](https://github.com/chaos-polymtl/lethe/pull/1842)

### [Master] - 2025-12-04

### Fixed

- MINOR New changes to the portable matrix-free architecture in dealii (https://github.com/dealii/dealii/pull/19042) created a warning due to the change from apply_for_each_quad_point to for_each_quad_point. This PR ports the Kokkos infrastructure to this and fixes the warnings. [#1838](https://github.com/chaos-polymtl/lethe/pull/1838)

- MINOR Tracer reaction models failed for reaction order below 1 due to undefined NaN values arising from negative powers of concentration, specifically at the initial condition of tracer = 0. We added a smoothed minimum value (`tracer reaction threshold`) that prevents this issue and a MMS test to ensure that the proper solution is recovered even at a `tracer reaction order = 0.5`. [#1825](https://github.com/chaos-polymtl/lethe/issues/1825)

- MINOR The lethe-fluid-vans example for the flow through a particles was not working adequately since it was using a semi-implicit drag formulation. This is now fixed and the example now uses the fully implicit coupling. The results are much better now since they compare quite favorably with Ergun empirical equation (and the bonus is that the example actually runs correctly!). [#1840](https://github.com/chaos-polymtl/lethe/pull/1840)

### [Master] - 2025-12-03

### Added

- MINOR This PR adds the option to fix pressure in a node within the matrix-based solver. The process is the same as what is already done in the coarse multigrid level of the matrix-free solver when fix pressure constant is set to true. [#1828](https://github.com/chaos-polymtl/lethe/pull/1828)

### Fixed

- MINOR The velocity used to compute the stabilization parameter (tau) did not account for the ALE velocity component when the mortar feature is enabled. This led to inconsistent stabilization terms in the rotor and stator sides, and it has been fixed in this PR. [#1834](https://github.com/chaos-polymtl/lethe/pull/1834)

### [Master] - 2025-12-01

### Added

- MINOR This PR adds the option of using the filtered particle-fluid forces in the VANS equations for the matrix-based CFD-DEM solver, following the logic implemented starting from PR [#1618]. In this PR, only an explicit coupling is implemented for the fluid drag and one test case is performed, which is that of a single sedimenting particle. Subsequent PRs including the implicit and semi-implicit drag couplings, as well as more application tests are to follow. [#1813](https://github.com/chaos-polymtl/lethe/pull/1813)

### [Master] - 2025-11-26

### Fixed

- MINOR As pointed out in issue [#1820] (https://github.com/chaos-polymtl/lethe/issues/1820), the `set_initial_condition_fd` function in `fluid_dynamics_matrix_based`, when the L2_projection boundary condition is apply, does not work if the solver is not GMRES. This PR solves this issue by introducing a new dedicated function `solve_L2_system` that will always use GMRES regardless of the choice of solver by the user. This choice is justified because the L2 projection corresponds to a trivial linear problem involving only a mass matrix, for which GMRES converges rapidly. [#1827](https://github.com/chaos-polymtl/lethe/issues/1827)

### Changed
- MINOR As an effort to remove the number of vectors dangling around all over the place in Lethe, we are trying to move everything to use deal.II parallel vectors instead of a weird blend between Trilinos and deal.II parallel vectors. A first step towards this that is within reach is to refactor the ParticleProjector since that class needs to work with both the matrix-based and the matrix-free architecture. This change refactors the fields functionality of the ParticleProjector to only work with deal.II parallel vectors. This greatly reduces the number of vectors and simplifies the workflow of the code, yet does not change anything in the results. A small consequence of this is that the CG solver used is now the deal.II CG solver instead of the Trilinos one, but they behave similarly. A follow-up change will be to do this complete refactor with void fraction within the ParticleProjector class. [#1826](https://github.com/chaos-polymtl/lethe/pull/1826)

### [Master] - 2025-11-25

### Added

- MAJOR We now need the possibility to solve linear system of equations for future physics, consequently the `physics_solver` and the `non_linear_solver` classes have now an overloaded constructor to enable this functionality. To improve the semantics with this new architecture, the general name  "solver_strategy" is used instead of "non_linear" where appropriate. A new `linear_solver_strategy` class is also introduced with its associated unit test. [#1804] (https://github.com/chaos-polymtl/lethe/pull/1804)

### Fixed

- MAJOR As identified in issue [#1518](https://github.com/chaos-polymtl/lethe/issues/1518), in MultiphysicsInterface solutions and DoFHandler were stored and shared through raw pointers. To avoid memory leak issues, these are now handled with shared pointers (`std::shared_ptr`). [#1823](https://github.com/chaos-polymtl/lethe/pull/1823)

### Added

- MINOR Add the SDIRK time integration method to the lethe-fluid-matrix-free solver. The SDIRK scheme is now available for both MatrixBased and MatrixFree solvers. However, it does not support multiphysics components at the moment. This will be added in a near future. [#1616](https://github.com/chaos-polymtl/lethe/pull/1616/)

### Fixed

- MINOR In the calculate_particle_fluid_interactions function of the VANS assemblers for lift, pressure, and shear forces, the expressions for forces applied to the fluid were corrected. This includes adjustments to the signs and division by the fluid density where needed in the computation of explicit_particle_volumetric_acceleration_on_fluid. For the pressure and shear forces, only model B of the vans equation is affected since model A does not use the explicit pressure and shear forces. [#1815](https://github.com/chaos-polymtl/lethe/pull/1815)

### [Master] - 2025-11-24

### Added

- MINOR Added option to override the checkpointed time-step upon restart. [#1824](https://github.com/chaos-polymtl/lethe/pull/1824)

### [Master] - 2025-11-23

### Removed

- MINOR Reynolds stress solutions were to be integrated into MultiphysicsInterface. The integration work began in 2022 but was never completed. Since they are not used through MultiphysicsInterface and the implementation that had been started was flawed, the relevant code elements have been removed. This can be reworked when it becomes necessary to integrate Reynolds stress solutions into MultiphysicsInterface. [#1822](https://github.com/chaos-polymtl/lethe/pull/1822)

### [Master] - 2025-11-22

### Fixed

- MINOR Contacts between particles and solid surfaces had a problem when the contact was occurring with more than one triangle from the same solid object at the same time. Contacts were duplicated even if the solid surface represent a single surface. This resulting in a particle contacting the same solid object twice, meaning that the particle would see a surface that has twice the stiffness at this location. In this PR, duplicated contacts are identified and removed following the logic described in https://doi.org/10.1002/nme.4487. This ensures that there is only one contact between the particle and a solid object and keeps the contact regular as a particle moves over a solid object. [#1745](https://github.com/chaos-polymtl/lethe/pull/1745)

### [Master] - 2025-11-21

### Added

- MAJOR Adds a new solver (cfd_dem_coupling_matrix_free) which uses the matrix-free architecture for the VANS equations and couples it with the discrete element method. Thus far, this solver is a lot more robust (and faster), especially when using high-order element. It is, however, a very experimental solver and it has not been fully tested nor validated. It should thus be used with a grain of salt. To integrate this solver, some refactoring had to be done in the particle projector and some slight issues were identified. There is some code duplication between the new cfd_dem solver and the previous one, but to avoid having a large PR that combines a new feature and refactoring, we will live with this code duplication for now and it will be refactored step-by-step in a near future. [#1801](https://github.com/chaos-polymtl/lethe/pull/1801)

### Fixed

- MINOR In MultiphysicsInterface, there were two methods with the exact same function `set_block_previous_solutions` and `set_previous_block_solutions`. `set_previous_block_solutions` was removed to avoid redundancy and improve code clarity. [#1816](https://github.com/chaos-polymtl/lethe/pull/1816)

### [Master] - 2025-11-20

### Fixed

- MINOR As noted in issue [#1503], the SimulationControl class contained a duplicated function to set the time-step. This was overly confusing and unclear. There is now only a single function that sets the time-step and this is all it does. It does not append the time-step to the time list, it just sets the time-step. Some documentation was added also on some methods and members of the class to enhance their readability. [#1812](https://github.com/chaos-polymtl/lethe/pull/1812)

### [Master] - 2025-11-15

### Fixed

- MINOR Issue [#1803] pointed out that the parameter `initial_step` in solve_linear_system function is not used. This solves this issue and removes the parameter. In addition, this parameter was coming upstream from solve_non_linear_system so it was removed from there also. [#1813](https://github.com/chaos-polymtl/lethe/pull/1813)

### [Master] - 2025-11-15

### Fixed

- MINOR Since PR [#1722], the sedimentation-1-cuboid example crashes when the particle collides with the bottom wall. The previous case used to run, but the non-linear and linear solver struggled pretty hard during the collision. I suspect that the case itself was actually unstable and that PR [#1722] did not break anything in particular. Using adaptive time-stepping with a max CFL of 0.5 resolves the issue, leads to the same results, and reduces the computational time. We will keep on monitoring this test case to ensure that it remains stable. [#1805](https://github.com/chaos-polymtl/lethe/pull/1805)

### [Master] - 2025-11-14

### Fixed

- MINOR Remove the deprecated argument "renewed_matrix" of the solve_linear_system definition. 

### [Master] - 2025-11-14

### Fixed

- MAJOR A bug had been introduced in [#1752] that prevented paraview from opening the vtu, pvtu and pvd files for the particle results in the CFD-DEM solver. This was because the new property (momentum_transfer_coefficient) was not given a name and this prevented the output from being adequately named. This PR fixes it by giving an appropriate name to that property which fixes the output of the particles. [#1800](https://github.com/chaos-polymtl/lethe/pull/1800)

### [Master] - 2025-11-04

### Added

- MAJOR The residual displayed in the non-linear and linear solver are affected by the total volume of the triangulation (mesh). Consequently, if a simulation is carried out on a triangulation with a small volume, very small residuals are obtained at the initial solution step. This complicates the choice of the linear and non-linear tolerances since they need to be adjusted in a case-dependent fashion. This change introduces a new parameter ("rescale residual") to the linear solver subsection. When the parameter is true, all residuals (linear and non-linear) are rescaled by the squared root of the volume of the triangulation. This is very convenient because, when activated, tolerances are independent of the domain size. [#1728](https://github.com/chaos-polymtl/lethe/pull/1728)

### [Master] - 2025-11-08

### Fixed
- MINOR A segmentation fault was appearing when using uniform refinement with velocity average postprocessing. This was because the velocity average postprocessor prepare_for_mesh_adaptation() function was never being called for the uniform mesh adaptation. A call was added and now the segfault is fixed. [#1797](https://github.com/chaos-polymtl/lethe/pull/1797)

### Added
- MAJOR The particle-fluid coupling in the VANS matrix-free solver only allowed for explicit particle-fluid force coupling. This is because the Jacobian matrix of the particle-fluid coupling is difficult to establish. This PR changes this by adding the capability of filtering the particle-fluid force and the momentum transfer coefficient between the particles and the fluid. The momentum transfer coefficient is used to establish the implicit particle-fluid coupling whereas in the explicit coupling mode, the drag force is applied directly onto the fluid. The explicit mode is much faster since it does not require solving as many non-linear iterations, but it has a constraint on the time step. That constraint is, however, not a severe limitation in the case where the fluid is a liquid. With this PR, full particle-fluid coupling within the matrix-free vans solver may be established. Many test cases were added to test this coupling with gas and with a dense fluid (water). [#1757](https://github.com/chaos-polymtl/lethe/pull/1757)

### [Master] - 2025-11-03

### Added
- MAJOR The particle-fluid coupling the VANS and CFD-DEM solver currently rely on a semi-implicit coupling where the drag force is calculated explicitely and partially implicited on the fluid side. This is a bit problematic for some cases and I want to introduce capability to control what we are essentially doing. This change introduces 3 coupling mode. The semi-implicit (which is the one we actually use right now) as well as a fully implicit and a fully explicit coupling mode. In the fully implicit mode, the particle-fluid forces are calculated using the full velocity solution at time t+delta t. Surprisingly, this is a lot more stable than I thought and it leads to very good results. In the fully explicit formulation, the particle-fluid forces are calculated at time t and are applied to both the particle and the fluid. In this case, there is no Jacobian matrix and the solver is fully explicit. This means that these is a stability criterion. To make this work well, I had to move the CFD-DEM parameters to the core solver so that the library may be used within the scratch data (since I need to know the drag mode in the calculation of the velocity at the particle location). I also found some slight issues with the two-way coupling for the lift forces and there were some signs that were confusing (for example the lift force was applied to both the fluid and the particle with the same sign and the - sign came later), this was kinda confusing and I took the opportunity to clarify this here. As a consequence, some tests results have slightly changed. [#1752](https://github.com/chaos-polymtl/lethe/pull/1752)

### [Master] - 2025-11-02

### Fixed

- MINOR  Two of the files using the std::numbers namespace were actually missing the "numbers" include. Ths created a compilation issue on Apple CLANG. This has been fixed by adding the correct include in the sdirk and dem files. [#1750](https://github.com/chaos-polymtl/lethe/pull/1750)

- MINOR Fixes a bug introduced in [#1731]. Essentially, the wrong function was used to write the boundaries checkpoint pvd file (read instead of save) and this would crash the simulation when trying to write a new checkpoint file from zero. This PR fixes this by calling the correct function. [](https://github.com/chaos-polymtl/lethe/pull/)

### [Master] - 2025-11-01

### Fixed

- MINOR A warning related to Trillinos and Kokkos was occuring when compiling Lethe with gcc. The warning was occuring in the "get_face_transformation_jacobian" function. This PR fixes this warning by adding an "if" statement. [#1748](https://github.com/chaos-polymtl/lethe/pull/1748)   

### Fixed

- MINOR Fixes the output of the void fraction at iteration 0 when using the CFD-DEM solver. The void fraction was not calculated before the first post-processing iteratio, which means that a zero void fraction was always displayed as the initial condition. This is now fixed and the void fraction is correctly calculated before the first output. [#1730](https://github.com/chaos-polymtl/lethe/pull/1730)

### [Master] - 2025-10-29

### Changed

- MINOR This PR changes the parameter use_manifold_for_normal, used in slip BCs, to false when mortar is enabled. This way only the mapping is used to compute normal vectors.[#1743](https://github.com/chaos-polymtl/lethe/pull/1743)

- MINOR Some DEM applications-tests results changed due to a modification related to the smallest solid object mapping criterion. This criterion now uses the std::numbers::inv_sqrt3 function instead of a hardcoded value. The position of the particles does not change, only the number of contact detection iteration. [#1740](https://github.com/chaos-polymtl/lethe/pull/1740)

### Added

- MINOR This PR adds ALE terms for transient mortar problems in the matrix-free solver. It is similar to what has been done for the matrix-based solver in #1597 and #1638.[#1744](https://github.com/chaos-polymtl/lethe/pull/1744)

- MINOR This PR adds a post processing function to write pvd files for boundaries when output boundaries is set to true in the simulation control section. This function is used when the mortar feature is enabled, so that rotating boundaries can be printed at every time step. [#1731](https://github.com/chaos-polymtl/lethe/pull/1731)

### [Master] - 2025-10-28

### Changed

- MINOR During [PR#1618](https://github.com/chaos-polymtl/lethe/pull/1618), the `solve_linear_system_and_update_solution` method in `PhysicsLinearSubequationsSolver` was renamed to `solve_void_fraction_linear_system`, which was relevant to the specific use of that PR. However, it also accidentally changed the method's name in the `VOFLinearSubequationsSolver` class, where the "void fraction" naming is less appropriate. This PR reverts that change. [#1742](https://github.com/chaos-polymtl/lethe/pull/1742)
- MINOR The normal of the triangle was calculated in LetheGridTools::find_point_triangle_distance but not being used. This removes this element. [#1736](https://github.com/chaos-polymtl/lethe/pull/1736)

### Added

- MAJOR A new application called lethe-particle-ray-tracing has been added. This application allows to trace rays through a domain containing particles and to extract the intersection point between those rays and particles. This application can be used generate numerical profilometry measurements which can be use to validate DEM simulations. [#1629](https://github.com/chaos-polymtl/lethe/pull/1629)

### [Master] - 2025-10-25

### Changed

- MINOR The array used to calculate the void fraction at the particle location was not being resized to the number of particles, but it was kept at a constant size corresponding to the maximum number of particles. This could lead to a weird segmentation fault in some edge cases that I have not found an easy way to reproduce. This is fixed by resizing the array to be consistent with the number of particles. [#1732](https://github.com/chaos-polymtl/lethe/pull/1732)

### [Master] - 2025-10-23

### Changed

- MINOR The serialization/deserialization of tables in classes derived from AuxiliaryPhysics have been refactored to match the same code structure as that used for the tables in solvers derived from NavierStokesBase. Functions for serializing/deserializing multiple tables are now defined in include/core/utilities.h. Also, the output_struct.h file has been moved to the core folder. [#1725](https://github.com/chaos-polymtl/lethe/pull/1725)

### [Master] - 2025-10-16

### Changed

- MINOR The scratch data for all of the FEM-based physics used std::vector<type> for both the evaluation of the solution at the gauss point, but also to store the shape function (or their gradient, divergence, etc.) at the gauss points for the degrees of freedom. The issue with this is that it creates std::vector<std::vector<type>> data structures. These vectors of vectors are not contiguous in memory and may lead to a large number of cache misses as the degree of the polynomial of the shape function is increased. This PR fixes this by changing the data structure used to store the shape function (and their gradient, divergence, etc.) to use the deal.II Table<n_dim,type> data structure. These data structures are identical to vectors of vectors in the way they are used, but they store the information in a contiguous memory block (see https://dealii.org/current/doxygen/deal.II/classTableBase.html for more information). This change has minimal impact on the code, but should lead to much better scalability of the matrix and RHS assembly when the degree of the polynomial is increased and the matrix-based solvers are used. The improvement will be problem dependent, but it should never worsen performance. [#1726](https://github.com/chaos-polymtl/lethe/pull/1726)

### [Master] - 2025-10-14

### Changed

- MAJOR The DG Tracer formulation for the equation - D * Laplacian T +  u * gradT - f=0 is formulated in the weak sense. Consequently, we solve it by solving -D * Laplacian T + div(u*T) - f = 0 and taking the weak form of that expression. However, this assumes that the divergence of the velocity field is exactly zero, yet the divergence of the velocity field is not exactly zero if the velocity field arises from a CG formulation. This is especially true near corners or hanging nodes. To fix this, we explicitly account for the divergence of the velocity field and solve for -D * Laplacian T + div(u*T) -T*div(u) - f = 0. This essentially explicitly removes the term T div(u). It has proven to be more conservative at edges and corners. The dg_tracer_pipe application test was sensitive to this change, and will be used to track it. The use of the new formulation seems more sensitive to the CFL to keep the solution closer to its bounds.

### Fixed

- MINOR In our efforts to migrate from Clang Tidy v12 to Clang Tidy v20, this PR does the first batch of refactoring of the core library of lethe. The main things that is changed is that many loops over the dimensions used to use unsigned integer as indices and, when compared to the dimension (dim), this would raise a warning regarding comparison between signed and unsigned integer. This PR fixes that by using integers for these loops. Other small things are to use enum with a smaller data footprint and some usage of range-based interators. Some spelling mistakes there and there were also fixed. A flag was also added to the CMakeList (-std=c++20) to ensure that compilation still works on Apple Clang. Overall this removes half of the CLang V20 warnings. [#1722](https://github.com/chaos-polymtl/lethe/pull/1722)


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
