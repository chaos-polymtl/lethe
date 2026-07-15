// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// 2D midplane cross-section of a particle-laden jet impinging on a wall.
// x-axis: flow direction (horizontal), y-axis: height along wall (vertical).
//
// The chamber now extends back to x = -dc so it surrounds the nozzle, modelling
// the ambient fluid entrained around the jet. The nozzle is a channel of width
// wc with solid walls of thickness tc on each side; those walls run from the
// back face (x = -dc) and end at the junction (x = 0). The wall bands are NOT
// meshed (un-meshed solid, no-slip BC). Impingement wall is at x = dp.
//
// y-bands (bottom -> top): outer-bottom | wall | nozzle | wall | outer-top
//   chamber  (0 .. dp):  all five bands are fluid (wall has ended)
//   back     (-dc .. 0): nozzle + two outer bands are fluid; wall bands = void

// Projection zone (main chamber)
hp = 0.2;    // width of chamber
dp = 0.05;   // height of chamber (jet-to-wall distance)
zp = 8;      // cells in flow direction (chamber)
nhp = 15;    // points across each chamber outer half-band (graded toward jet)
exp = 1.05;  // grading progression toward jet axis

// Jet nozzle (charging zone)
wc = 0.0072; // nozzle width
dc = 0.075;  // nozzle length
zc = 21;     // cells in flow direction (nozzle / back region)
nhc = 5;     // points across nozzle half-width
tc = 0.001;  // nozzle wall thickness (solid, un-meshed). 2D: flat wall, no sqrt(2)
ntc = 3;     // points across each nozzle wall band
wo = wc/2 + tc;  // outer-wall offset from axis

// ---- Points: a 3 (x) by 6 (y) grid -------------------------------------------
// y-levels: -hp/2, -wo, -wc/2, wc/2, wo, hp/2

// Back / inlet column (x = -dc)
Point(1) = {-dc, -hp/2, 0};
Point(2) = {-dc, -wo,   0};
Point(3) = {-dc, -wc/2, 0};
Point(4) = {-dc,  wc/2, 0};
Point(5) = {-dc,  wo,   0};
Point(6) = {-dc,  hp/2, 0};

// Junction column (x = 0)
Point(7)  = {0, -hp/2, 0};
Point(8)  = {0, -wo,   0};
Point(9)  = {0, -wc/2, 0};
Point(10) = {0,  wc/2, 0};
Point(11) = {0,  wo,   0};
Point(12) = {0,  hp/2, 0};

// Impingement-wall column (x = dp)
Point(13) = {dp, -hp/2, 0};
Point(14) = {dp, -wo,   0};
Point(15) = {dp, -wc/2, 0};
Point(16) = {dp,  wc/2, 0};
Point(17) = {dp,  wo,   0};
Point(18) = {dp,  hp/2, 0};

// ---- Lines -------------------------------------------------------------------
// Vertical (constant x). At x = -dc only the fluid-facing segments exist; the
// wall-band back faces (P2-P3, P4-P5) are omitted (they bound only solid).
Line(1) = {1, 2};    // anterior face, bottom (x = -dc)
Line(2) = {3, 4};    // inlet (nozzle, x = -dc)
Line(3) = {5, 6};    // anterior face, top (x = -dc)

Line(4) = {7, 8};    // junction, outer-bottom band
Line(5) = {8, 9};    // junction, bottom wall band  -> nozzle-wall rim
Line(6) = {9, 10};   // junction, nozzle band
Line(7) = {10, 11};  // junction, top wall band     -> nozzle-wall rim
Line(8) = {11, 12};  // junction, outer-top band

Line(9)  = {13, 14}; // impingement, outer-bottom band
Line(10) = {14, 15}; // impingement, bottom wall band
Line(11) = {15, 16}; // impingement, nozzle band
Line(12) = {16, 17}; // impingement, top wall band
Line(13) = {17, 18}; // impingement, outer-top band

// Horizontal (constant y). Back region (x in [-dc,0]) then chamber (x in [0,dp]).
Line(14) = {1, 7};   // y = -hp/2  (side wall, back)
Line(15) = {2, 8};   // y = -wo    (outer nozzle wall, bottom)
Line(16) = {3, 9};   // y = -wc/2  (inner nozzle wall, bottom)
Line(17) = {4, 10};  // y =  wc/2  (inner nozzle wall, top)
Line(18) = {5, 11};  // y =  wo    (outer nozzle wall, top)
Line(19) = {6, 12};  // y =  hp/2  (side wall, back)

Line(20) = {7, 13};  // y = -hp/2  (side wall, chamber)
Line(21) = {8, 14};  // y = -wo
Line(22) = {9, 15};  // y = -wc/2
Line(23) = {10, 16}; // y =  wc/2
Line(24) = {11, 17}; // y =  wo
Line(25) = {12, 18}; // y =  hp/2  (side wall, chamber)

// ---- Surfaces ----------------------------------------------------------------
// Back region (x in [-dc, 0]): nozzle + two outer bands (wall bands are void)
Line Loop(1) = {14, 4, -15, -1};   // surrounding fluid, bottom
Line Loop(2) = {16, 6, -17, -2};   // nozzle channel
Line Loop(3) = {18, 8, -19, -3};   // surrounding fluid, top

// Chamber (x in [0, dp]): all five bands are fluid
Line Loop(4) = {20, 9, -21, -4};   // outer-bottom band
Line Loop(5) = {21, 10, -22, -5};  // bottom wall band (fluid here; void behind)
Line Loop(6) = {22, 11, -23, -6};  // nozzle band
Line Loop(7) = {23, 12, -24, -7};  // top wall band (fluid here; void behind)
Line Loop(8) = {24, 13, -25, -8};  // outer-top band

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};

// ---- Transfinite meshing -----------------------------------------------------
// Flow direction (x)
Transfinite Line {14, 15, 16, 17, 18, 19} = zc + 1 Using Progression 1; // back (dc)
Transfinite Line {20, 21, 22, 23, 24, 25} = zp + 1 Using Progression 1; // chamber (dp)

// Across the nozzle channel (y) — uniform
Transfinite Line {2, 6, 11} = nhc Using Progression 1;

// Across the nozzle wall bands (y) — uniform, thin
Transfinite Line {5, 7, 10, 12} = ntc Using Progression 1;

// Chamber outer half-bands (y) — graded fine toward the jet
// Bottom lines run outer(-hp/2)->inner(-wo): reverse so the fine end is inner
Transfinite Line {-1, -4, -9} = nhp Using Progression exp;
// Top lines run inner(wo)->outer(hp/2): forward keeps the fine end at inner
Transfinite Line {3, 8, 13} = nhp Using Progression exp;

Transfinite Surface {1, 2, 3, 4, 5, 6, 7, 8};
Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8};

// ---- Physical groups (boundary IDs match the 3D convention) ------------------
Physical Surface(0) = {1, 2, 3, 4, 5, 6, 7, 8};  // fluid domain

Physical Line(0) = {2};                  // Inlet (nozzle back face, x = -dc)
Physical Line(1) = {1, 3};               // Anterior face (box back wall around nozzle)
Physical Line(2) = {9, 10, 11, 12, 13};  // Impingement wall (x = dp)
Physical Line(3) = {15, 16, 17, 18, 5, 7};  // Nozzle walls: outer(15,18) inner(16,17) rim(5,7)
Physical Line(4) = {14, 20, 19, 25};     // Side walls (top & bottom, full length)


Mesh.ElementOrder = 1;
Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;