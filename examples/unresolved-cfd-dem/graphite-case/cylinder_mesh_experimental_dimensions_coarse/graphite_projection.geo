// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Define variables

// Projection zone

hp = 0.25; // height of projection zone
wp = 0.25; // width of projection zone
dp = 0.02; // depth in z of projection zone
zp = 10; // number of cells in the z direction (depth) of projection zone
nhp = 27; // number of points in height of section of projection zone
nwp = 27; // number of points in width of section of projection zone
exp = 1.05;

// Charging zone

wc = 0.007; // width/height of charging zone
sq = 0.0054; // width/height of square inside circle
dc = 0.036; // depth in z of charging zone
nhc = 5; // number of points in height/width of charging zone
nwc = 3; // number of points in the width/height of the square
zc = 20; // number of cells in the z direction (depth) of charging zone
R = Sqrt(2) * wc; // radius of big circle
h = (R - wc) / 2; // height of arcs

Point(1) = {-wp/2, hp/2, 0};
Point(2) = {-wp/2, wc/2, 0};
Point(3) = {-wp/2, -wc/2, 0};
Point(4) = {-wp/2, -hp/2, 0};
Point(5) = {-wc/2, hp/2, 0};
Point(6) = {-wc/2, wc/2, 0};
Point(7) = {-wc/2, -wc/2, 0};
Point(8) = {-wc/2, -hp/2, 0};
Point(9) = {wc/2, hp/2, 0};
Point(10) = {wc/2, wc/2, 0};
Point(11) = {wc/2, -wc/2, 0};
Point(12) = {wc/2, -hp/2, 0};
Point(13) = {wp/2, hp/2, 0};
Point(14) = {wp/2, wc/2, 0};
Point(15) = {wp/2, -wc/2, 0};
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

Line(1) = {1, 2};
Line(2) = {2, 6};
Line(3) = {6, 5};
Line(4) = {5, 1};
Line(5) = {3, 2};
Circle(6) = {6, 22, 7};
Line(7) = {7, 3};
Line(8) = {3, 4};
Line(9) = {4, 8};
Line(10) = {8, 7};
Line(11) = {5, 9};
Line(12) = {9, 10};
Circle(13) = {10, 22, 6};
Circle(14) = {7, 22, 11};
Circle(15) = {11, 22, 10};
Line(16) = {11, 12};
Line(17) = {12, 8};
Line(18) = {10, 14};
Line(19) = {14, 13};
Line(20) = {13, 9};
Line(21) = {14, 15};
Line(22) = {15, 11};
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

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 2, 6, 7};
Line Loop(3) = {8, 9, 10, 7};
Line Loop(4) = {10, 14, 16, 17};
Line Loop(5) = {16, 23, 24, 22};
Line Loop(6) = {15, 18, 21, 22};
Line Loop(7) = {12, 18, 19, 20};
Line Loop(8) = {3, 11, 12, 13};

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

Transfinite Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
Recombine Surface{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};


Extrude {0, 0, dp} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{9}; Surface{10}; Surface{11}; Surface{12}; Surface{13}; Layers{zp}; Recombine;
}

Extrude {0, 0, dc} {
  Surface{14}; Surface{15}; Surface{16}; Surface{17}; Surface{18}; Layers{zc}; Recombine;
}

Transfinite Line {5, 6, 15, 21, 25, 27, 11, 13, 14, 17, 28, 26, 29, 30, 31, 32, 37, 38, 40, 39} = Ceil(nhc) Using Progression 1;
Transfinite Line {41, 44, 43, 42, 33, 36, 34, 35} = Ceil(nwc) Using Progression 1;
Transfinite Line {4, -2, 7, -9, -20, 18, -22, 23} = Ceil(nwp) Using Progression exp;
Transfinite Line {-1, 3, -12, 19, 8, -10, 16, -24} = Ceil(nhp) Using Progression exp;


Physical Volume(0) = {1:18};

Physical Surface(0) = {14, 15, 16, 17, 18};
Physical Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Surface(2) = {110, 88, 66, 220, 198, 176, 154, 132, 264, 330, 242, 286, 308};
Physical Surface(3) = {402, 389, 364, 427};
Physical Surface(4) = {193, 171, 149, 197, 211, 65, 97, 75, 53, 101, 131, 145};