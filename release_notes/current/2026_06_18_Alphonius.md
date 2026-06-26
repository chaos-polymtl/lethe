## [Master] - 2026/06/18

### Fixed

- MAJOR For box refinements, the transformations of translation and rotation were applied in the wrong order, resulting in a less intuitive transformation. We reorder them here so that the rotation comes first and then the translation. The existing application tests were updated to include a combined transformation of rotation, translation and scaling. The maximum number of box refinements is also increased to 20. [#2038](https://github.com/chaos-polymtl/lethe/pull/2038)