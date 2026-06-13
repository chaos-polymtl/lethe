## [Master] - 2026/06/11

### Added

- MAJOR Adds the Carman-Kozeny non-linear phase change permeability model (``carman-kozeny phase change``). The model is added in both the single fluid and the NS-CLS assemblers. Both source terms are tested through application tests (phase_change_carman_kozeny and  cls_phase_change_carman_kozeny). The model is set through the new ``permeability model`` parameter that replaces the now deprecated ``Darcy type``  parameter. The model introduces two new parameters: ``Carman-Kozeny permeability area`` and ``Carman-Kozeny division tolerance`` in the ``velocity source`` subsection. The user documentation was updated with the new model and changes to the parameters. [#2013](https://github.com/chaos-polymtl/lethe/pull/2013)

### Fixed

- MAJOR Fixes the RHS of the Darcy source term in the ``PhaseChangeDarcyCLSAssembler``. The Darcy penalty was precalculated and appended to a vector, but the value added to the source term was wrongly indexed. [#2013](https://github.com/chaos-polymtl/lethe/pull/2013)
