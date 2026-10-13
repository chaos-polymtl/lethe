## [Master] - 2026/10/07

### Fixed

- MINOR Multiple checkpointing functions were performed by every rank in parallel, which is not needed. This is problematic for many reason, one of them being that the checkpoint_controller file was updated by rank N while rank 0 did not finish writing all the other checkpoint files. As a result, the controller file was corrupted and had to be changed manually. [#2150](https://github.com/chaos-polymtl/lethe/pull/2150)
