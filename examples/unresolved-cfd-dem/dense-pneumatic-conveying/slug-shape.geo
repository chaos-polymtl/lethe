// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

//+
SetFactory("OpenCASCADE");
lc = 1;
r = 0.042;
L = 1;
Ls = 0.5;
layer = 0.021;
buffer = 0.01;


// Before slug
Point(0) = {-L/2,      -r+layer,  r, lc};
Point(1) = {-L/2,      -r+layer, -r, lc};
Point(2) = {-(Ls/2+2*r-layer), -r+layer, -r, lc};
Point(3) = {-(Ls/2+2*r-layer), -r+layer,  r, lc};

Line(0)={0,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,0};

Line(100)={2,0};

Line Loop(1) = {1,0,100};
Line Loop(100) = {3,2,100};

// Rear slope of slug
Point(4) = {-(Ls/2), r, -r, lc};
Point(5) = {-(Ls/2), r,  r, lc};

Line(4)={2,4};
Line(5)={4,5};
Line(6)={5,3};

Line(200)={3,4};

Line Loop(2) = {6,5,200};
Line Loop(200) = {2,4,200};

// Top slug
Point(6) = {(Ls/2), r, -r, lc};
Point(7) = {(Ls/2), r,  r, lc};

Line(7)={4,6};
Line(8)={6,7};
Line(9)={7,5};

Line(300)={5,6};

Line Loop(3) = {8,9,300};
Line Loop(300) = {5,7,300};

// Front slope of slug
Point(8) = {(Ls/2+2*r-layer), -r+layer, -r, lc};
Point(9) = {(Ls/2+2*r-layer), -r+layer,  r, lc};

Line(10)={7,9};
Line(11)={9,8};
Line(12)={8,6};

Line(400)={8,7};

Line Loop(4) = {10,11,400};
Line Loop(400) = {8,12,400};

// Front stationnary layer
Point(10) = {L/2,      -r+layer,  r, lc};
Point(11) = {L/2,      -r+layer, -r, lc};

Line(13)={9,10};
Line(14)={10,11};
Line(15)={11,8};

Line(500)={11,9};

Line Loop(5) = {13,14,500};
Line Loop(500) = {11,15,500};

Plane Surface(1) = {1};
Plane Surface(100) = {100};
Plane Surface(2) = {2};
Plane Surface(200) = {200};
Plane Surface(3) = {3};
Plane Surface(300) = {300};
Plane Surface(4) = {4};
Plane Surface(400) = {400};
Plane Surface(5) = {5};
Plane Surface(500) = {500};
Physical Surface(0) = {1,100,2,200,3,300,4,400,5,500};

