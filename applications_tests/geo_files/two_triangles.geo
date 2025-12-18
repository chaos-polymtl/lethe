// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Description : This file defines a mesh of two triangles that are sharing an edge.
// The two triangles are not coplanar, nor perpendicular.

// Define a variable
ymin= -1.;
ymax= 0.;
zmax= 1;
zmin= -1;
x0= -0.01;
x1= -0.01;
x2= 0.5;

lc = 4.;

Point(0)= {x0, -2., zmin, lc};
Point(1)= {x1, ymax, zmin, lc};
Point(2)= {x1, ymax, zmax, lc};
Point(3)= {x2, ymin, 0, lc};


Line(0)={0,1};
Line(1)={1,2};
Line(2)={2,0};
Line(3)={2,3};
Line(4)={3,1};

Line Loop(1) = {0,1,2};
Line Loop(2) = {1,3,4};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Surface(0) = {1,2};

