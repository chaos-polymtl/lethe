// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

gap = 0.0001;

Point(0) = {gap-0.007, gap, 0, 1.0};
Point(1) = {gap-0.007, 0.1-gap, 0, 1.0};
Point(2) = {-0.05+gap, 0.1-gap, 0, 1.0};
Point(3) = {-0.05+gap, gap, 0, 1.0};

Line(0)={0,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,0};

Line Loop(1) = {0,1,2,3};

Plane Surface(1) = {1} ;

Physical Surface(0) = {1};