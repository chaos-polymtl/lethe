// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Variables
// h : jet height
// d : distance from bottom
// Lin : inlet length (development flow)
// Lout : outlet length
h = 0.01;
d = h/2;
Lin = 8*h;
Lout = 100*h;

// Mesh parameters
cell_size = 1.0;
points_v_Lin = 41; // nodes along x-axis in inlet
points_v_Lout = 501; // nodes along x-axis in outlet
points_h_constricted = 11; // nodes along y-axis in constricted section
points_h_outer = 6; // nodes along y-axis in expanded sections

// Points
Point(1) = {Lin, 0, 0, cell_size};
Point(2) = {Lin+Lout, 0, 0, cell_size};
Point(3) = {Lin+Lout, d, 0, cell_size};
Point(4) = {Lin+Lout, d+h, 0, cell_size};
Point(5) = {Lin+Lout, d+h+d, 0, cell_size};
Point(6) = {Lin, d+h+d, 0, cell_size};
Point(7) = {Lin, d+h, 0, cell_size};
Point(8) = {0, d+h, 0, cell_size};
Point(9) = {0, d, 0, cell_size};
Point(10) = {Lin, d, 0, cell_size};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 1};
Line(11) = {7, 4};
Line(12) = {10, 3};
Line(13) = {7, 10};

// Surfaces
Curve Loop(1) = {1, 2, -12, 10};
Plane Surface(1) = {1};
Curve Loop(2) = {9, -13, 7, 8};
Plane Surface(2) = {2};
Curve Loop(3) = {12, 3, -11, 13};
Plane Surface(3) = {3};
Curve Loop(4) = {11, 4, 5, 6};
Plane Surface(4) = {4};
Recombine Surface {1, 2, 3, 4};

// Quad meshing (transfinite)
Transfinite Surface{1} = {1, 2, 3, 10};
Transfinite Surface{2} = {9, 10, 7, 8};
Transfinite Surface{3} = {10, 3, 4, 7};
Transfinite Surface{4} = {7, 4, 5, 6};
Transfinite Curve{10, 2} = points_h_outer Using Progression 1;
Transfinite Curve{8, 13, 3} = points_h_constricted Using Progression 1;
Transfinite Curve{6, 4} = 5 Using Progression 1;
Transfinite Curve{1, 12, 11, 5} = points_v_Lout Using Progression 1;
Transfinite Curve{7, 9} = points_v_Lin Using Progression 1;

// Physical lines
// Boundaries
Physical Line(0) = {8}; // inlet
Physical Line(1) = {2, 3, 4}; // outlet
Physical Line(2) = {9, 10, 1, 5, 6, 7}; // walls

// Surface
Physical Surface(0) = {1, 2, 3, 4};
Recombine Surface {1, 2, 3, 4};






