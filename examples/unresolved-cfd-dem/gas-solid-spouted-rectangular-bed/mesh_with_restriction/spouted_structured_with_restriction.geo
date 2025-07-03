// SPDX-FileCopyrightText: Copyright (c) 2022 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Define a variable

H=1; // height of bed
W=0.128; // half width of bed from end of channel to end of bed
Z = 0.04; // depth of the bed in z
iW=0.012; // half width of channel
B=0.04; // point at channel inlet
nw=18;  // number of points  in the half width of bed after channel
nwi=3; // number of points in width of channel
nb=6; // number of points in the height of the channel
nh=101; // number of points in the height of the bed without channel
zl = 4; // number of cells in the z direction (depth)
T = 0.2; // height of top part
nt = 21;
L = 0.07; // half width of top part
nl = 10;

Point(2) = {-iW, -B, 0};
Point(3) = {+iW, -B, 0};
Point(11) = {-W-iW, 0, 0};
Point(12) = {-iW, 0, 0};
Point(13) = {iW,0, 0};
Point(14) = {iW+W, 0, 0};
Point(21) = {-W-iW, H, 0};
Point(22) = {-iW, H, 0};
Point(23) = {iW,H, 0};
Point(24) = {iW+W, H, 0};
Point(31) = {-L, T+H, 0};
Point(32) = {L, T+H, 0};
Point(33) = {-iW, T+H, 0};
Point(34) = {+iW, T+H, 0};

Line(1)={11,12};
Line(2)={12,2};
Line(3)={2,3};
Line(4)={3,13};
Line(5)={13,14};
Line(6)={14,24};
Line(7)={24,23};
Line(8)={23,22};
Line(9)={22,21};
Line(10)={21,11};

Line(101)={12,13};
Line(102)={12,22};
Line(103)={13,23};
Line(104)={21,31};
Line(105)={31,33};
Line(106)={33,34};
Line(107)={34,32};
Line(108)={32,24};
Line(109)={22,33};
Line(110)={23,34};


Line Loop(1) = {1,102,9,10};
Line Loop(2) = {101,103,8,-102};
Line Loop(3) = {5,6,7,-103};
Line Loop(4) = {2,3,4,-101};
Line Loop(5) = {104, 105, -109, 9};
Line Loop(6) = {109, 106, -110, 8};
Line Loop(7) = {110, 107, 108, 7};


Plane Surface(1) = {1} ;
Plane Surface(2) = {2} ;
Plane Surface(3) = {3} ;
Plane Surface(4) = {4} ;
Plane Surface(5) = {5} ;
Plane Surface(6) = {6} ;
Plane Surface(7) = {7} ;

Transfinite Surface {1,2,3,4,5,6,7};
Recombine Surface{1,2,3,4,5,6,7};

Extrude {0, 0, Z} {
  Surface{1}; Surface{3}; Surface{2}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Layers{zl}; Recombine;
}


Transfinite Line {10,102,103,6} = Ceil(nh) Using Progression 1;
Transfinite Line {1,5,9,7, 107, 105} = Ceil(nw) Using Progression 1;
Transfinite Line {3,101,8,106} = Ceil(nwi) Using Progression 1;
Transfinite Line {2,4} = Ceil(nb) Using Progression 1;
Transfinite Line {104,108,109,110} = Ceil(nt) Using Progression 1;

Physical Volume(0) = {1:7};

Physical Surface(0) = {1,2,3,4,5,6,7,131,132,145,154,176,185,193,198,207,220,242,259,264}; // sides
Physical Surface(1) = {189}; // channel inlet
Physical Surface(2) = {211,233,255}; // outlet
Physical Surface(3) = {119,141}; // bed inlet (background velocity)