// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

//+
SetFactory("OpenCASCADE");
lc = 1;
W = 0.2;
L = 0.45;
H = 0.0;


// Left side
Point(0) = {-L, H, 0, lc};
Point(1) = {-L, H, W, lc};
Point(2) = {L,  H, W, lc};
Point(3) = {L,  H, 0, lc};

Line(0)={0,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,0};

Line(100)={2,0};

Line Loop(1) = {1,0,100};
Line Loop(100) = {3,2,100};

Plane Surface(1) = {1};
Plane Surface(100) = {100};

Physical Surface(0) = {1,100};

