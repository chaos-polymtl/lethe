# Contact between two sedimenting spheres (drafting, kissing, tumbling)

Reproduction of section 4.3 and Fig. 12 of L. Barbeau, S. Golshan, J. Deng,
S. Étienne, C. Béguin and B. Blais, *High-order moving immersed boundary and its
application to a resolved CFD-DEM model*, Computers & Fluids 268 (2024) 106094.
The case is from Glowinski et al., J. Comput. Phys. 169 (2001) 363–426, §8.5.

This case is part of the Lucille validation suite
(`contrib/validation/validation_cases_lucille.txt`). It is not yet a documented
example: this file is its documentation.

## Running

`lethe-fluid-sharp` must be in the PATH. As a validation case, which runs Q1Q1
only (128 cores):

```bash
bash validate.sh -p <absolute output path> -np 128
```

One discretization, with its post-processing:

```bash
mpirun -np 128 lethe-fluid-sharp sedimentation-2-particles-q1q1-h128.prm
python3 post-process-sedimentation-2-particles.py -d q1q1-h128

mpirun -np 128 lethe-fluid-sharp sedimentation-2-particles-q2q1-h128.prm
python3 post-process-sedimentation-2-particles.py -d q2q1-h128
```

| File | Discretization | Cells at t = 0 | DoFs at t = 0 |
|---|---|---|---|
| `sedimentation-2-particles-q1q1-h128.prm` (validation) | Q1Q1, h = 1/128 | 67k | 0.31 M |
| `sedimentation-2-particles-q2q1-h128.prm` (reference) | Q2Q1, h = 1/128 | 67k | 1.8 M |

The two files differ only by the degree of the velocity, the stencil and the
output path. Q2Q1 needs much more memory than its number of DoFs suggests: more
than 35 GB before the end of the first linear solve on 14 ranks. Since the CFL
divides the cell size by the degree of the velocity, Q2Q1 also takes about
twice as many time steps. A Q2Q1 step costs roughly 20 times a Q1Q1 step, so
the whole Q2Q1 simulation is estimated at 4000 to 5000 core-hours, against
about 110 for Q1Q1 (7.7 h on 14 ranks of a Ryzen 9 3950X). Q2Q1 is therefore
not part of the validation for now. Add `q2q1-h128` to `discretizations` in
`validate.sh` to include it.

## Definition of the case (g, cm, s)

| Quantity | Value | Origin |
|---|---|---|
| Tank | 1 x 1 x 4 | paper |
| Spheres | diameter 1/6, density 1.14 | paper |
| p0 (upper) | (0.5 + 1/300, 0.5, 3.5) | paper |
| p1 (lower) | (0.5, 0.5, 3.16) | paper, see below |
| Fluid | density 1, kinematic viscosity 0.01 | paper (0.001 Pa s) |
| Gravity | (0, 0, -980) | paper |
| Contact | E = 1e12 dyn/cm² (100 GPa), Poisson 0.3, restitution 1, no friction | paper |
| Time integration | BDF2, initial time step 0.001, then adapted to a CFL of 0.8 | the paper uses a constant time step of 0.001 |
| Refinement around the spheres | inside factor 0.5, outside factor 1.5 | paper |
| Lubrication range max | 2 | paper (lambda = 2) |

**Initial height of p1.** Glowinski et al. give 3.1666... (3.16 with an overbar:
the centers are two diameters apart). The overbar was lost in our paper, and the
simulations did use 3.16: integrating the velocities of Fig. 12 closes a gap of
0.1732–0.1740 cm before contact, against 0.1733 for 3.16 and 0.1667 for 3.1666...
This case keeps 3.16 to be comparable with the paper.

**Diameter.** Text extracted from the PDF of the paper reads "16 cm". It is
1/6 cm, as in the LaTeX source and in Fig. 11.

## Choices that the paper does not give

| Choice | Value | Reason |
|---|---|---|
| Top of the tank | `function_weak`, zero velocity, beta = 1 | no-slip like Glowinski, but well-posed: the pressure level is fixed |
| Other walls | `noslip` | Glowinski |
| Base mesh | h = 1/16 | |
| Kelly refinement | velocity, every step, 5 % refined, 30 % coarsened | `sedimentation-1-particle` |
| Maximum number of cells | 750k | safety cap |
| Stencil | degree of the velocity, length ratio 2 (Q1) or 4 (Q2) | §2.2.1 of the paper and parameter documentation |
| Lubrication range min | 0.1 | default, also in November 2022 |
| DEM sub-steps | 10000 per time step (they scale with the adapted time step) | for a time step of 0.001, about 60–130 sub-steps per collision at 0.1–3 cm/s. The default of 1000 gives 8–14, and up to 0.4 % error on the restitution |
| Contact properties of the walls | those of the spheres | |
| End time | 0.7 | the paper stops when a sphere reaches the bottom (0.56–0.64 s) |
| Solver tolerances | Newton 1e-7, particles 1e-5, GMRES 1e-5 relative and 1e-9 minimum | |

Lethe multiplies the lubrication range by the diagonal of the smallest cell
divided by sqrt(2), which is 1.22 h in 3D: a range max of 2 is thus 2.45 h. The
code did the same when the results of the paper were produced.

## Differences to expect with the paper

The figure of the paper was created on 2022-11-09. Compared with master of
2022-11-08 (e7f5b651b):

- **Lubrication.** That code multiplied the force by the distance between the
  centers (about 1/6 cm at contact) and evaluated `3 / 2` as an integer. The
  force was about 1/9 of the correct one, which `sharp-edge-bugfix-3` computes.
  No parameter reproduces the old force, so kissing and tumbling should differ.
- **Contact.** That code used the full overlap in the Hertz model. The use of
  half the overlap came later and is fixed by `sharp-edge-bugfix-3`, so the
  contact stiffness is again that of the paper.
- **Chaos.** Tumbling depends on the discretization (see the paper), and can
  depend on round-off and on the number of ranks. Judge the drafting and the
  kissing quantitatively, the tumbling qualitatively. Always validate on the
  same number of cores.

## Post-processing

`validate.sh` produces, in `<output path>/sedimentation-2-particles/`, for
every discretization `<d>` that it runs (`q1q1-h128` for now):

- `log-<d>`: output of the simulation;
- `two-particles-velocity-<d>.pdf`: vertical velocities against the paper (same
  discretization dashed, the other one dotted) and the literature;
- `two-particles-trajectories-<d>.pdf`: gap, horizontal and vertical positions;
- `solution-particle-{0,1}-<d>.dat`: velocity and position every 0.01 s up to 0.7 s
  (NaN if not reached), so that the number of entries never changes;
- `solution-metrics-<d>.dat`: value of Lethe and value of the paper of the same
  discretization, with the same definitions, for

| Metric | Definition |
|---|---|
| `peak_velocity_p0` | minimum v_z of p0 for t ≤ 0.4 |
| `kissing_time` | first time the vertical gap is below 1/128 cm |
| `overtaking_time` | first time p0 is below p1 |
| `bottom_time` | first time a sphere is within 1/256 cm of the bottom. For the paper, end of the curves |
| `minimum_gap` | smallest distance between the surfaces, negative for an overlap |
| `rms_error_paper_p*` | RMS difference with the curve of the paper for t ≤ 0.28 (drafting) |
| `rms_error_<reference>_p*` | same with the literature. Second column: error of the paper |

## Reference data

`reference/` comes from the original vector file of Fig. 12
(`validation_contact.pdf`, in `article 2/im/` of the LaTeX sources of the thesis
of L. Barbeau). The coordinates are read from the PDF and calibrated on its
ticks.

- The Lethe curves are the plotted ones, as simplified by matplotlib. Only the
  two discretizations with h = 1/128 are kept; the figure also has Q1Q1 1/256.
- Nie & Lin (2020) and Apte et al. (2009) are **second-hand**: they are the
  markers digitized for the paper, not data from the original articles.
- In the figure, the p0 markers of both references end with a copy of the last
  p1 markers (hidden under them). They are removed: 7 markers at t = 0.569–0.611
  for Nie & Lin, 3 markers at t = 0.592–0.624 for Apte et al.

## Turning it into an example

1. Write the documentation page in
   `doc/source/examples/sharp-immersed-boundary/sedimentation-2-particles/`
   from this file, list it in the `toctree` of the category, and remove this file.
2. Add a release note.
