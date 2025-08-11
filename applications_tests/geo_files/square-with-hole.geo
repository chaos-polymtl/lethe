// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Variables
xc=1; // x center
yc=1; // y center
r=1; // radius
// Auxiliary distances
srq=r;
Ls=3*r;
ss=Ls;

// Mesh parameters
nr=3; // refinement radial
nt=3; // refinement tangential
lc = 1; // cell size

// Points
// Outside
Point(1) = {xc-ss, yc-ss, 0, lc};
Point(2) = {xc+ss, yc-ss, 0, lc};
Point(3) = {xc+ss, yc+ss, 0, lc};
Point(4) = {xc-ss, yc+ss, 0, lc};
// Circle
Point(5) = {xc-srq, yc-srq, 0, lc};
Point(6) = {xc+srq, yc-srq, 0, lc};
Point(7) = {xc+srq, yc+srq, 0, lc};
Point(8) = {xc-srq, yc+srq, 0, lc};
// Center
Point(9) = {xc, yc, 0, lc};

// Lines
// Outside
Line(10) = {1,4};
Line(11) = {4,3};
Line(12) = {2,3};
Line(13) = {1,2};
//Circle
Circle(14)={8,9,7};
Circle(15)={7,9,6};
Circle(16)={6,9,5};
Circle(17)={5,9,8};
// Diagonals
Line(18)={8,4};
Line(19)={7,3};
Line(20)={6,2};
Line(21)={5,1};

//Theta of the circle
Transfinite Line {11,12,13,10,14,15,16,17} = nt Using Progression 1.0;
// Radial direction
Transfinite Line {18,19,20,21} = nr Using Progression 1.0;

// Line loops
Line Loop(22) = {-18,14,19,-11};
Line Loop(23) = {-10,-21,17,18};
Line Loop(24) = {21,13,-20,16};
Line Loop(25) = {20,12,-19,15};

// Surfaces
Plane Surface(22) = {22};
Plane Surface(23) = {23};
Plane Surface(24) = {24};
Plane Surface(25) = {25};

// Transfinite mapping
Transfinite Surface {22:25};
Recombine Surface {22:25};

// Physical surfaces
Physical Surface(1) = {22:25};

// Physical boundaries
Physical Line(0) = {10, 11, 12, 13}; // outside
Physical Line(1) = {14, 15, 16, 17}; // inside

