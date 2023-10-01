# Change Log
All notable changes to the Lethe project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## [Master] - 2023-09-30
  
### Fixed

- MINOR There was an error in the application of the affine constraints used to clip the void fraction in the lethe-fluid-vans and lethe-fluid-particles solvers. This led to assertions being thrown in debug. This has been corrected by reiniting the constraints with the appropriate size.

## [Master] - 2023-09-20
  
### Deprecated

- MINOR The calculation of the source term was enabled using a parameter called "enable". This parameter was used in some physics, not in others and was poorly implemented. We deprecate the usage of this parameter and always enable source term, considering the fact that the default value of the source term is zero anyway. This prevent the false perception that source terms could be enabled or disabled, while the behavior was inconsistent across physics.

 
## [Master] - 2023-09-19
  
### Changed

- MAJOR All the applications were renamed [#882](https://github.com/lethe-cfd/lethe/pull/882): `gls_navier_stokes` is now `lethe-fluid`, `gd_navier_stokes` is now `lethe-fluid-block`, `nitsche_navier_stokes` is now `lethe-fluid-nitsche`, `gls_sharp_navier_stokes` is now `lethe-fluid-sharp`, `gls_vans` is now `lethe-fluid-vans`, `mf_navier_stokes` is now `lethe-fluid-matrix-free`, `dem` is now `lethe-particles`, `cfd_dem_coupling` is now `lethe-fluid-particles`, `rpt3d` is now `lethe-rpt-3d`, `rpt_cell_reconstruction_3d` is now `lethe-rpt-cell-reconstruction-3d`, `rpt_fem_reconstruction_3d` is now `lethe-rpt-fem-reconstruction-3d`, and `rpt_l2_projection_3d` is now `lethe-rpt-l2-projection-3d`. 

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