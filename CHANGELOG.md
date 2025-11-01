
# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

### [Master] - 2025-11-01

### Fixed

- MINOR A warning related to Trillinos and Kokkos was occuring when compiling Lethe with gcc. The warning was occuring in the "get_face_transformation_jacobian" function. This PR fixes this warning by adding an "if" statement. [#1748](https://github.com/chaos-polymtl/lethe/pull/1748)   


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
