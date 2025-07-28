// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

H=0.5; // height of bed
W=0.09; // width of bed
D=0.008; // depth of bed

Point(1) = {0, 0, -0.04};
Point(2) = {W, 0, -0.04};
Point(3) = {W, D, -0.04};
Point(4) = {0, D, -0.04};

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

out[] = Extrude {0, 0, H+0.04} {
  Surface{1};
  Layers{50};
  Recombine;
};

Physical Surface(0) = {1};                        // Bottom
Physical Surface(1) = {out[0]};                   // Top
// Sides : 2 and 4: front and back, 3 and 5: left and right
Physical Surface(2) = {out[2],out[4]};
Physical Surface(3) = {out[3],out[5]};
Physical Volume(3)  = {out[1]};
Transfinite Volume {out[1]};