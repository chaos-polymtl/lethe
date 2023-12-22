
# Change Log
All notable changes to the Lethe project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## [Master] - 2023-12-22

### Fixed

- MINOR  PVD handles for solid DEM objects were written with the wrong text format ("_" seperarators instead of ".")

### [Master] - Changed - 2023-12-17

- MINOR The "insertion random number range" and "insertion random number seed" parameters got renamed to "insertion maximum offset" and "insertion prn seed" respectively. The old names didn't make sense, as they're not random because they're defined in the parameter file. [#970](https://github.com/lethe-cfd/lethe/pull/896)


## [Master] - 2023-12-11

### Fixed

- MINOR  The number of remaining particle to insert of each type is being checkpointed adequately. This means that no modification are required to the "number of particle" parameter after restarting a simulation. [#964](https://github.com/lethe-cfd/lethe/pull/964) 

### Fixed

- MINOR Solid objects can now be restarted adequately in DEM. They will resume at the position they had at the end of the simulation [#959](https://github.com/lethe-cfd/lethe/pull/959) 

## [Master] - 2023-11-27

### Removed

- MINOR The average diameter of the uniform size distribution with the DEM module was specified with a "average diameter" paramter. It is now specified directly from the diameter parameter. This is correctly documented. [#940](https://github.com/lethe-cfd/lethe/pull/940).

### Fixed

- MINOR The DEM time step verification was outputting the the most permissive time step (the biggest) and not the most restrive (the smallest). This bugfix doesn`t affect the uniform particle size simulation. [#939](https://github.com/lethe-cfd/lethe/pull/939).


## [Master] - 2023-11-23

### Fixed

- MINOR The plane insertion for the DEM was only supporting the uniform diameter distribution. Now it supports all types of distribution. 

## [Master] - 2023-11-16
  
### Changed

- MINOR The maximum number of boundary conditions for all physics was fixed to 14 since the boundary conditions had the be declared before being parsed. A new mechanism is now in place which parses the "number" parameter for each physics and keeps the maximal value. Then, this maximal value is used to pre-declare the boundary conditions. This enables much more robust sanity checking of the input parameter. The major drawback of this (and this is a major one) is that if we ever have another parameter with the name "number" then this parameter would also be parsed and used to establish the maximum number of boundary conditions. In this case, the best approach would be to replace "number" with "number of boundary conditions" in the parameter file. I (B-saurus-rex) did not want to do this at the time of this change to not have a massive PR which breaks every parameter files.

- MAJOR The "number" parameter within "subsection lagrangian physical properties" and "particle type n" was changed to "number of particles" to prevent confusions with the "number" used for boundary conditions. The "number" for boundary conditions will be changed to "number of boundary conditions" in the near future.

## [Master] - 2023-11-12
  
### Deprecated

- MINOR The uniform insertion method had been removed. The non-uniform insertion method has been renamed to volume method to remain coherent with the plane method. If you want to use an insertion method equivalent to the uniform insertion method, use the volume method with a "insertion random number range " equal to zero. [#926](https://github.com/lethe-cfd/lethe/pull/926)


## [Master] - 2023-10-01
  
### Fixed

- MINOR The calculation of the velocity on rotating walls was calculated inadequately [#896](https://github.com/lethe-cfd/lethe/pull/896): The velocity of rotating boundary conditions (e.g. boundaries of the mesh) was calculated inadequately in the case where the normal vector of the wall was aligned with the rotation axis. The whole calculation procedure was slightly messed up and only worked for cylinders. This has been fixed and made general, paving the way for a full refactor of the calculation of the particle-wall contact force.

## [Master] - 2023-09-30
  
### Fixed

- MINOR Affine constraints used to bound the void fraction and used as boundary conditions within the heat transfer and the block Navier-stokes solver were corrected [#885](https://github.com/lethe-cfd/lethe/pull/895): There was an error in the application of the affine constraints used to clip the void fraction in the lethe-fluid-vans and lethe-fluid-particles solvers. This led to assertions being thrown in debug. This has been corrected by reiniting the constraints with the appropriate size. For the heat transfer and the block Navier-stokes, the issue was that the constraints were never reinited with the correct size, so they contained ghosted elements. This was caught by a new assert introduced in 2023-09 within deal.II master.

## [Master] - 2023-09-20
  
### Deprecated

- MINOR The calculation of the source term was enabled using a parameter called "enable". This parameter was used in some physics, not in others and was poorly implemented. We deprecate the usage of this parameter and always enable source term, considering the fact that the default value of the source term is zero anyway. This prevent the false perception that source terms could be enabled or disabled, while the behavior was inconsistent across physics.

 
## [Master] - 2023-09-19
  
### Changed

- MAJOR All the applications were renamed [#882](https://github.com/lethe-cfd/lethe/pull/882): `gls_navier_stokes` is now `lethe-fluid`, `gd_navier_stokes` is now `lethe-fluid-block`, `nitsche_navier_stokes` is now `lethe-fluid-nitsche`, `gls_sharp_navier_stokes` is now `lethe-fluid-sharp`, `gls_vans` is now `lethe-fluid-vans`, `mf_navier_stokes` is now `lethe-fluid-matrix-free`, `dem` is now `lethe-particles`, `cfd_dem_coupling` is now `lethe-fluid-particles`, `rpt3d` is now `lethe-rpt-3d`, `rpt_cell_reconstruction_3d` is now `lethe-rpt-cell-reconstruction-3d`, `rpt_fem_reconstruction_3d` is now `lethe-rpt-fem-reconstruction-3d`, and `rpt_l2_projection_3d` is now `lethe-rpt-l2-projection-3d`. 


## [Master] - 2023-10-30

- MINOR The rotational vector for the rotational boundary condition in the lethe-particles solver is now define with one line in the parameters file. [#920](https://github.com/lethe-cfd/lethe/pull/920)

 

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


