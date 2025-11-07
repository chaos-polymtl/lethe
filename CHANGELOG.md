
# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

### [Master] - 2025-11-04

### Added

- MAJOR The residual displayed in the non-linear and linear solver are affected by the total volume of the triangulation (mesh). Consequently, if a simulation is carried out on a triangulation with a small volume, very small residuals are obtained at the initial solution step. This complicates the choice of the linear and non-linear tolerances since they need to be adjusted in a case-dependent fashion. This change introduces a new parameter ("rescale residual") to the linear solver subsection. When the parameter is true, all residuals (linear and non-linear) are rescaled by the sqrt of the volume of the triangulation. This is very convenient because, when activated, tolerances are independent of the domain size. [1728](https://github.com/chaos-polymtl/lethe/pull/1728)

### [Master] - 2025-11-08

### Fixed
- MINOR A segmentation fault was appearing when using uniform refinement with velocity average postprocessing. This was because the velocity average postprocessor prepare_for_mesh_adaptation() function was never being called for the uniform mesh adaptation. A call was added and now the segfault is fixed. [#1797](https://github.com/chaos-polymtl/lethe/pull/1797)

### [Master] - 2025-11-03

### Added
- MAJOR The particle-fluid coupling the VANS and CFD-DEM solver currently rely on a semi-implicit coupling where the drag force is calculated explicitely and partially implicited on the fluid side. This is a bit problematic for some cases and I want to introduce capability to control what we are essentially doing. This change introduces 3 coupling mode. The semi-implicit (which is the one we actually use right now) as well as a fully implicit and a fully explicit coupling mode. In the fully implicit mode, the particle-fluid forces are calculated using the full velocity solution at time t+delta t. Surprisingly, this is a lot more stable than I thought and it leads to very good results. In the fully explicit formulation, the particle-fluid forces are calculated at time t and are applied to both the particle and the fluid. In this case, there is no Jacobian matrix and the solver is fully explicit. This means that these is a stability criterion. To make this work well, I had to move the CFD-DEM parameters to the core solver so that the library may be used within the scratch data (since I need to know the drag mode in the calculation of the velocity at the particle location). I also found some slight issues with the two-way coupling for the lift forces and there were some signs that were conmfusing (for example the lift force was applied to both the fluid and the particle with the same sign and the - sign came later), this was kinda confusing and I took the opportunity to clarify this here. As a consequence, some tests results have slightly changed. [#1752](https://github.com/chaos-polymtl/lethe/pull/1752)

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
