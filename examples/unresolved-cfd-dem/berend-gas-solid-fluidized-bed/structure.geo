// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

H=0.5; // height of bed
W=0.09; // width of bed
D=0.008; // depth of bed

Point(1) = {0, 0, 0};
Point(2) = {W, 0, 0};
Point(3) = {W, D, 0};
Point(4) = {0, D, 0};

Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 1};

Transfinite Line {1, 3} = 10 Using Progression 1;
Transfinite Line {2, 4} = 1 Using Progression 1;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface {1};

out[] = Extrude {0, 0, H} {
  Surface{1};
  Layers{50};
  Recombine;
};

Physical Surface(0) = {1};                        // Bottom
Physical Surface(1) = {out[0]};                   // Top
Physical Surface(2) = {out[2], out[3], out[4], out[5]}; // Sides
Physical Volume(3)  = {out[1]};
Transfinite Volume {out[1]};