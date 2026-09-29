## [Master] - 2026/09/29

### Added

- MAJOR This PR adds the post-processing of the force and torque exerted by the particles on the solid surfaces in ``lethe-particles``. The force and torque are accumulated from the contact forces already calculated for the particle-solid contacts, and the torque is calculated about the center of rotation of each solid surface. They are enabled with the new ``solid forces`` subsection of the ``post-processing`` subsection, and written in one file per solid surface at the requested output frequency. The bunny drill example now outputs the force and torque on the bunny. [#XXXX](https://github.com/chaos-polymtl/lethe/pull/XXXX)
