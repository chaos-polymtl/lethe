// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Define variables

// Charging zone (nozzle)

wc = 0.0072 / Sqrt(2); // width/height of charging zone
sq = 0.75 * wc; // width/height of square inside circle
dc = 0.075; // depth in z of charging zone
tc = 0.00125; // thickness of nozzle tube
nhc = 5; // number of points in height/width of charging zone
nwc = 3; // number of points in the width/height of the square
ntc = 2; // number of points across the nozzle wall thickness (annulus)
wo = wc/2 + tc/Sqrt(2); // half-width of the nozzle outer wall (outer O-grid ring)
zc = 8; // number of cells in the z direction (depth) of charging zone
R = Sqrt(2) * wc; // radius of big circle
h = (R - wc) / 2; // height of arcs


// Projection zone (main chamber)

hp = 0.2; // height of projection zone
wp = 0.2; // width of projection zone
dp = 0.05; // depth in z of projection zone
dt = dc + dp; // total depth of chamber
hb = hp/2 - wc/2 - tc/Sqrt(2); // height of anterior segment of chamber
zp = 6; // number of cells in the z direction (depth) of projection zone
nhp = 15; // number of points in height of section of projection zone
nwp = 15; // number of points in width of section of projection zone
exp = 1.05;


Point(1) = {-wp/2, hp/2, 0};
Point(2) = {-wp/2, wo, 0};
Point(3) = {-wp/2, -wo, 0};
Point(4) = {-wp/2, -hp/2, 0};
Point(5) = {-wo, hp/2, 0};
Point(6) = {-wc/2, wc/2, 0};
Point(7) = {-wc/2, -wc/2, 0};
Point(8) = {-wo, -hp/2, 0};
Point(9) = {wo, hp/2, 0};
Point(10) = {wc/2, wc/2, 0};
Point(11) = {wc/2, -wc/2, 0};
Point(12) = {wo, -hp/2, 0};
Point(13) = {wp/2, hp/2, 0};
Point(14) = {wp/2, wo, 0};
Point(15) = {wp/2, -wo, 0};
Point(16) = {wp/2, -hp/2, 0};
Point(17) = {-wc/2, wc/2, -dc};
Point(18) = {-wc/2, -wc/2, -dc};
Point(19) = {wc/2, wc/2, -dc};
Point(20) = {wc/2, -wc/2, -dc};
Point(21) = {0, 0, -dc};
Point(22) = {0, 0, 0};
Point(23) = {-sq/2, sq/2, 0};
Point(24) = {-sq/2, -sq/2, 0};
Point(25) = {sq/2, -sq/2, 0};
Point(26) = {sq/2, sq/2, 0};
Point(27) = {-sq/2, sq/2, -dc};
Point(28) = {-sq/2, -sq/2, -dc};
Point(29) = {sq/2, -sq/2, -dc};
Point(30) = {sq/2, sq/2, -dc};

// Outer-wall O-grid corner points (z = 0), radially outside the inner-circle
// corners 6, 7, 11, 10. These define the outer circle of the nozzle wall.
Point(31) = {-wo, wo, 0};   // outer of 6 (top-left)
Point(32) = {-wo, -wo, 0};  // outer of 7 (bottom-left)
Point(33) = {wo, -wo, 0};   // outer of 11 (bottom-right)
Point(34) = {wo, wo, 0};    // outer of 10 (top-right)

Line(1) = {1, 2};
Line(2) = {2, 31};
Line(3) = {31, 5};
Line(4) = {5, 1};
Line(5) = {3, 2};
Circle(6) = {6, 22, 7};
Line(7) = {32, 3};
Line(8) = {3, 4};
Line(9) = {4, 8};
Line(10) = {8, 32};
Line(11) = {5, 9};
Line(12) = {9, 34};
Circle(13) = {10, 22, 6};
Circle(14) = {7, 22, 11};
Circle(15) = {11, 22, 10};
Line(16) = {33, 12};
Line(17) = {12, 8};
Line(18) = {34, 14};
Line(19) = {14, 13};
Line(20) = {13, 9};
Line(21) = {14, 15};
Line(22) = {15, 33};
Line(23) = {12, 16};
Line(24) = {16, 15};
Circle(25) = {17, 21, 18};
Circle(26) = {18, 21, 20};
Circle(27) = {20, 21, 19};
Circle(28) = {19, 21, 17};
Line(29) = {23, 24};
Line(30) = {24, 25};
Line(31) = {25, 26};
Line(32) = {26, 23};
Line(33) = {23, 6};
Line(34) = {10, 26};
Line(35) = {25, 11};
Line(36) = {7, 24};
Line(37) = {27, 28};
Line(38) = {28, 29};
Line(39) = {29, 30};
Line(40) = {30, 27};
Line(41) = {27, 17};
Line(42) = {19, 30};
Line(43) = {29, 20};
Line(44) = {18, 28};

// Outer circle of the nozzle wall (z = 0), concentric with the inner circle
Circle(100) = {31, 22, 32}; // left
Circle(101) = {32, 22, 33}; // bottom
Circle(102) = {33, 22, 34}; // right
Circle(103) = {34, 22, 31}; // top

// Radial connectors across the wall thickness (inner corner -> outer corner)
Line(104) = {6, 31};
Line(105) = {7, 32};
Line(106) = {11, 33};
Line(107) = {10, 34};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 2, 100, 7};
Line Loop(3) = {8, 9, 10, 7};
Line Loop(4) = {10, 101, 16, 17};
Line Loop(5) = {16, 23, 24, 22};
Line Loop(6) = {102, 18, 21, 22};
Line Loop(7) = {12, 18, 19, 20};
Line Loop(8) = {3, 11, 12, 103};

Line Loop(9) = {29, 30, 31, 32};

Line Loop(10) = {33, -13, 34, 32};
Line Loop(11) = {34, -31, 35, 15};
Line Loop(12) = {35, -14, 36, 30};
Line Loop(13) = {36, -29, 33, 6};

Line Loop(14) = {37, 38, 39, 40};

Line Loop(15) = {41, -28, 42, 40};
Line Loop(16) = {42, -39, 43, 27};
Line Loop(17) = {43, -26, 44, 38};
Line Loop(18) = {44, -37, 41, 25};

// Nozzle-wall annulus segments (between inner and outer circle, z = 0).
// These are extruded forward (fluid in the projection zone) but NOT backward,
// so their backward sweep is the un-meshed solid nozzle wall.
Line Loop(19) = {6, 105, -100, -104};   // left
Line Loop(20) = {14, 106, -101, -105};  // bottom
Line Loop(21) = {15, 107, -102, -106};  // right
Line Loop(22) = {13, 104, -103, -107};  // top

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};
Plane Surface(14) = {14};
Plane Surface(15) = {15};
Plane Surface(16) = {16};
Plane Surface(17) = {17};
Plane Surface(18) = {18};
Plane Surface(19) = {19};
Plane Surface(20) = {20};
Plane Surface(21) = {21};
Plane Surface(22) = {22};

Transfinite Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};
Recombine Surface{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};


// Each Extrude below groups all of its base surfaces into a single command so
// shared edges stay conformal. The returned list is captured: every base
// surface is a quad, so it occupies a 6-entry block [top, volume, lat0, lat1,
// lat2, lat3], where lat_k is the swept k-th curve of that surface's Line Loop.
// Block i therefore starts at index 6*i (in the order the surfaces are listed).

// Projection zone (forward, +dp): full cross-section = outer ring (1-8)
// + nozzle-wall annulus (19-22) + inner circle (9-13). Here the wall has
// ended, so the annulus region is open fluid.
// Order: 1..13, 19..22  (i = 0..16)
pA[] = Extrude {0, 0, dp} {
  Surface{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 19, 20, 21, 22};
  Layers{zp}; Recombine;
};

// Charging zone / nozzle interior (backward, -dc): inner circle only.
// Order: 14..18  (i = 0..4)
pB[] = Extrude {0, 0, dc} {
  Surface{14, 15, 16, 17, 18};
  Layers{zc}; Recombine;
};

// Surrounding fluid drafted around the nozzle (backward, -dc): outer ring
// only. The annulus (1-8 vs inner circle) is deliberately NOT extruded here,
// leaving the solid nozzle wall as a void between this volume and the nozzle
// interior. The cross-section at z = 0 matches surfaces 1-8 exactly, so this
// is conformal with the projection zone.
// Order: 1..8  (i = 0..7)
pC[] = Extrude {0, 0, -dc} {
  Surface{1, 2, 3, 4, 5, 6, 7, 8};
  Layers{zc}; Recombine;
};

Transfinite Line {5, 6, 15, 21, 25, 27, 11, 13, 14, 17, 28, 26, 29, 30, 31, 32, 37, 38, 40, 39, 100, 101, 102, 103} = Ceil(nhc) Using Progression 1;
Transfinite Line {41, 44, 43, 42, 33, 36, 34, 35} = Ceil(nwc) Using Progression 1;
Transfinite Line {104, 105, 106, 107} = Ceil(ntc) Using Progression 1;
Transfinite Line {4, -2, 7, -9, -20, 18, -22, 23} = Ceil(nwp) Using Progression exp;
Transfinite Line {-1, 3, -12, 19, 8, -10, 16, -24} = Ceil(nhp) Using Progression exp;


// Fluid volumes: 17 (projection) + 5 (nozzle interior) + 8 (surrounding) = 30,
// created sequentially by the three Extrudes. The nozzle wall is never
// extruded, so it is absent from every volume (un-meshed solid).
Physical Volume(0) = {1:30};

// --- Boundary groups, addressed through the captured Extrude lists so the IDs
//     survive any geometry renumbering. Index = 6*block + slot, where slot 0 is
//     the top face and slots 2-5 are the swept Line Loop curves (in order).

// (0) Inlet: nozzle interior back face (z = -dc), stable base surfaces.
Physical Surface(0) = {14, 15, 16, 17, 18};

// (1) Anterior face of the chamber (z = -dc, the face the nozzle exits through):
//     top faces of the backward outer-ring extrusion pC, cells 1-8.
Physical Surface(1) = {pC[0], pC[6], pC[12], pC[18], pC[24], pC[30], pC[36], pC[42]};

// (2) Downstream wall (z = +dp): every top face of the forward extrusion pA.
Physical Surface(2) = {pA[0], pA[6], pA[12], pA[18], pA[24], pA[30], pA[36], pA[42],
                       pA[48], pA[54], pA[60], pA[66], pA[72], pA[78], pA[84],
                       pA[90], pA[96]};

// (3) Nozzle wall (no-slip): inner cylinder + outer cylinder + rim.
//   inner cylinder = swept inner arcs 28,27,26,25  (pB cells 15,16,17,18)
//   outer cylinder = swept outer arcs 100,101,102,103  (pC cells 2,4,6,8)
//   rim            = annulus end caps at z = 0  (stable base surfaces 19-22)
Physical Surface(3) = {pB[9], pB[17], pB[21], pB[29],
                       pC[10], pC[21], pC[32], pC[47],
                       19, 20, 21, 22};

// (4) Lateral chamber walls (the four long box sides, parallel to the axis):
//     swept rectangle-edge curves of cells 1-8, from BOTH extrusions pA and pC.
Physical Surface(4) = {pA[2], pA[5], pA[8], pA[14], pA[15], pA[23], pA[27], pA[28],
                       pA[34], pA[40], pA[41], pA[45],
                       pC[2], pC[5], pC[8], pC[14], pC[15], pC[23], pC[27], pC[28],
                       pC[34], pC[40], pC[41], pC[45]};
