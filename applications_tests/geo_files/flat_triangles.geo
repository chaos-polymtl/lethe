// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Description : This file defines two parallel triangles that are sharing an edge. The two triangles are making a square.

// Define a variable

ymin= 0.;
ymax= 0.;
zmax= 1;
zmin= -1;
x0= -1;
x1= 0;
x2= 1;

lc = 2.;

Point(0)= {x0, ymin, 0, lc};
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

