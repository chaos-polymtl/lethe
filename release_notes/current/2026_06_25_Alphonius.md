## [Master] - 2026/06/25

### Added

- MAJOR Implements the ability track ``temperature`` and ``phase`` (CLS phase indicator) isocontours in time. This can be useful to track for example the interface between two fluids or the melt pool depth for LPBF simulations. The ``isocontour bounding box`` subsection containing all relevant parameters is added to the ``post-processing`` subsection. The method to compute the bounding values is implemented within ``InterfaceTools`` and it is deployed in the ``HeatTransfer`` and ``ConservativeLevelSet`` classes. Three lethe-fluid application regression tests were modified to test this new feature: one with isocontours on the ``temperature`` only, one with isocontours on the ``phase`` only, and one with isocontours on both the ``temperature`` and the ``phase``. [#2028](https://github.com/chaos-polymtl/lethe/pull/2028)
