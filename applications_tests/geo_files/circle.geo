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
// Center
Point(9) = {xc, yc, 0, lc};
// Circle
Point(27) = {xc-srq, yc-srq, 0, lc};
Point(28) = {xc+srq, yc-srq, 0, lc};
Point(29) = {xc+srq, yc+srq, 0, lc};
Point(30) = {xc-srq, yc+srq, 0, lc};

// Lines
Circle(31)={30,9,29};
Circle(32)={29,9,28};
Circle(33)={28,9,27};
Circle(34)={27,9,30};

//Theta of the circle
Transfinite Line {31,32,33,34} = nt Using Progression 1.0;

// Line loops
Line Loop(26) = {31, 32, 33, 34};

// Surfaces
Plane Surface(26) = {26};

// Transfinite mapping
Transfinite Surface {26};
Recombine Surface {26};

// Physical surfaces
Physical Surface(1) = {26};

// Physical boundaries
Physical Line(2) = {31, 32, 33, 34};

