// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Define variables

// Projection zone

hp = 0.18; // height of projection zone
wp = 0.18; // width of projection zone
dp = 0.018; // depth in z of projection zone
zp = 20; // number of cells in the z direction (depth) of projection zone
nhp = 31; // number of points in height of section of projection zone
nwp = 31; // number of points in width of section of projection zone
exp = 1.05;

// Charging zone

hc = 0.0054; // height of charging zone
wc = 0.0054; // width of charging zone
dc = 0.036; // depth in z of charging zone
nhc = 10; // number of points in height of charging zone
nwc = 10; // number of points in the width of the charging zone
zc = 40; // number of cells in the z direction (depth) of charging zone

Point(1) = {-wp/2, hp/2, 0};
Point(2) = {-wp/2, hc/2, 0};
Point(3) = {-wp/2, -hc/2, 0};
Point(4) = {-wp/2, -hp/2, 0};
Point(5) = {-wc/2, hp/2, 0};
Point(6) = {-wc/2, hc/2, 0};
Point(7) = {-wc/2, -hc/2, 0};
Point(8) = {-wc/2, -hp/2, 0};
Point(9) = {wc/2, hp/2, 0};
Point(10) = {wc/2, hc/2, 0};
Point(11) = {wc/2, -hc/2, 0};
Point(12) = {wc/2, -hp/2, 0};
Point(13) = {wp/2, hp/2, 0};
Point(14) = {wp/2, hc/2, 0};
Point(15) = {wp/2, -hc/2, 0};
Point(16) = {wp/2, -hp/2, 0};
Point(17) = {-wc/2, hc/2, -dc};
Point(18) = {-wc/2, -hc/2, -dc};
Point(19) = {wc/2, hc/2, -dc};
Point(20) = {wc/2, -hc/2, -dc};

Line(1) = {1, 2};
Line(2) = {2, 6};
Line(3) = {6, 5};
Line(4) = {5, 1};
Line(5) = {3, 2};
Line(6) = {6, 7};
Line(7) = {7, 3};
Line(8) = {3, 4};
Line(9) = {4, 8};
Line(10) = {8, 7};
Line(11) = {5, 9};
Line(12) = {9, 10};
Line(13) = {10, 6};
Line(14) = {7, 11};
Line(15) = {11, 10};
Line(16) = {11, 12};
Line(17) = {12, 8};
Line(18) = {10, 14};
Line(19) = {14, 13};
Line(20) = {13, 9};
Line(21) = {14, 15};
Line(22) = {15, 11};
Line(23) = {12, 16};
Line(24) = {16, 15};
Line(25) = {17, 18};
Line(26) = {18, 20};
Line(27) = {20, 19};
Line(28) = {19, 17};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 2, 6, 7};
Line Loop(3) = {8, 9, 10, 7};
Line Loop(4) = {10, 14, 16, 17};
Line Loop(5) = {16, 23, 24, 22};
Line Loop(6) = {15, 18, 21, 22};
Line Loop(7) = {12, 18, 19, 20};
Line Loop(8) = {3, 11, 12, 13};
Line Loop(9) = {6, 14, 15, 13};
Line Loop(10) = {25, 26, 27, 28};

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

Transfinite Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Recombine Surface{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

Extrude {0, 0, dp} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{9}; Layers{zp}; Recombine;
}

Extrude {0, 0, dc} {
  Surface{10}; Layers{zc}; Recombine;
}

Transfinite Line {5, 6, 15, 21, 25, 27} = Ceil(nhc) Using Progression 1;
Transfinite Line {11, 13, 14, 17, 28, 26} = Ceil(nwc) Using Progression 1;
Transfinite Line {4, -2, 7, -9, -20, 18, -22, 23} = Ceil(nwp) Using Progression exp;
Transfinite Line {-1, 3, -12, 19, 8, -10, 16, -24} = Ceil(nhp) Using Progression exp;


Physical Volume(0) = {1:10};

Physical Surface(0) = {10};
Physical Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Surface(2) = {182,160,138,94, 72, 50, 204, 116, 226};
Physical Surface(3) = {243, 239, 235, 247};
Physical Surface(4) = {81, 59, 37, 49, 195, 181, 177, 155, 133, 129, 115, 85};