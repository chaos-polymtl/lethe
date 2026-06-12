// ============================================================
// 3D rectangular impinging-jet chamber
//
// Jet direction: positive x
//
// Boundary IDs:
//   1 = inlet at x = -Lp
//   2 = outlet at y = -Ly/2
//   3 = outlet at y = +Ly/2
//   4 = all walls
//
// Volume/material ID:
//   10 = computational domain
// ============================================================

SetFactory("Built-in");

// ------------------------------------------------------------
// Geometry dimensions
// ------------------------------------------------------------

D = 1.0;

// Chamber dimensions
Lx = 6.0 * D;
Ly = 10.0 * D;
Lz = 6.0 * D;

// Square inlet-pipe dimensions
Dy = 1.0 * D;
Dz = 1.0 * D;

// Inlet-pipe length
Lp = 2.0 * D;

// Chamber coordinates
x0 = 0.0;
x1 = Lx;

y0 = -Ly / 2.0;
y1 = -Dy / 2.0;
y2 =  Dy / 2.0;
y3 =  Ly / 2.0;

z0 = -Lz / 2.0;
z1 = -Dz / 2.0;
z2 =  Dz / 2.0;
z3 =  Lz / 2.0;

// ------------------------------------------------------------
// Mesh resolution
// ------------------------------------------------------------

// Chamber cells in x
Nx = 4;

// Pipe cells in x
Np = 3;

// Cells in each outer y-region
NySide = 3;

// Cells across the inlet in y
NyJet = 3;

// Cells in each outer z-region
NzSide = 3;

// Cells across the inlet in z
NzJet = 3;

lc = 1.0;

// ------------------------------------------------------------
// Points on the plane x = 0
//
// z = z3    13 ----- 14 ----- 15 ----- 16
//            |        |        |        |
// z = z2     9 ----- 10 ----- 11 ----- 12
//            |        | inlet  |        |
// z = z1     5 ------ 6 ------ 7 ------ 8
//            |        |        |        |
// z = z0     1 ------ 2 ------ 3 ------ 4
//
//           y0       y1       y2       y3
//
// Surface 5 is the pipe-to-chamber opening.
// ------------------------------------------------------------

// Row z = z0
Point(1) = {x0, y0, z0, lc};
Point(2) = {x0, y1, z0, lc};
Point(3) = {x0, y2, z0, lc};
Point(4) = {x0, y3, z0, lc};

// Row z = z1
Point(5) = {x0, y0, z1, lc};
Point(6) = {x0, y1, z1, lc};
Point(7) = {x0, y2, z1, lc};
Point(8) = {x0, y3, z1, lc};

// Row z = z2
Point(9)  = {x0, y0, z2, lc};
Point(10) = {x0, y1, z2, lc};
Point(11) = {x0, y2, z2, lc};
Point(12) = {x0, y3, z2, lc};

// Row z = z3
Point(13) = {x0, y0, z3, lc};
Point(14) = {x0, y1, z3, lc};
Point(15) = {x0, y2, z3, lc};
Point(16) = {x0, y3, z3, lc};

// ------------------------------------------------------------
// Lines in the y-direction
// ------------------------------------------------------------

// z = z0
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};

// z = z1
Line(4) = {5, 6};
Line(5) = {6, 7};
Line(6) = {7, 8};

// z = z2
Line(7) = {9, 10};
Line(8) = {10, 11};
Line(9) = {11, 12};

// z = z3
Line(10) = {13, 14};
Line(11) = {14, 15};
Line(12) = {15, 16};

// ------------------------------------------------------------
// Lines in the z-direction
// ------------------------------------------------------------

// y = y0
Line(13) = {1, 5};
Line(14) = {5, 9};
Line(15) = {9, 13};

// y = y1
Line(16) = {2, 6};
Line(17) = {6, 10};
Line(18) = {10, 14};

// y = y2
Line(19) = {3, 7};
Line(20) = {7, 11};
Line(21) = {11, 15};

// y = y3
Line(22) = {4, 8};
Line(23) = {8, 12};
Line(24) = {12, 16};

// ------------------------------------------------------------
// Nine surfaces on x = 0
// ------------------------------------------------------------

// Bottom z-region
Curve Loop(1) = {1, 16, -4, -13};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 19, -5, -16};
Plane Surface(2) = {2};

Curve Loop(3) = {3, 22, -6, -19};
Plane Surface(3) = {3};

// Middle z-region
Curve Loop(4) = {4, 17, -7, -14};
Plane Surface(4) = {4};

Curve Loop(5) = {5, 20, -8, -17};
Plane Surface(5) = {5}; // Pipe opening

Curve Loop(6) = {6, 23, -9, -20};
Plane Surface(6) = {6};

// Top z-region
Curve Loop(7) = {7, 18, -10, -15};
Plane Surface(7) = {7};

Curve Loop(8) = {8, 21, -11, -18};
Plane Surface(8) = {8};

Curve Loop(9) = {9, 24, -12, -21};
Plane Surface(9) = {9};

// ------------------------------------------------------------
// Structured surface mesh
// ------------------------------------------------------------

// Outer y-regions
Transfinite Curve {
  1, 3,
  4, 6,
  7, 9,
  10, 12
} = NySide + 1;

// Inlet y-region
Transfinite Curve {
  2, 5, 8, 11
} = NyJet + 1;

// Outer z-regions
Transfinite Curve {
  13, 15,
  16, 18,
  19, 21,
  22, 24
} = NzSide + 1;

// Inlet z-region
Transfinite Curve {
  14, 17, 20, 23
} = NzJet + 1;

Transfinite Surface {
  1, 2, 3,
  4, 5, 6,
  7, 8, 9
};

Recombine Surface {
  1, 2, 3,
  4, 5, 6,
  7, 8, 9
};

// ------------------------------------------------------------
// Extrude the chamber in the positive x-direction
// ------------------------------------------------------------

chamber[] = Extrude {Lx, 0, 0}
{
  Surface {
    1, 2, 3,
    4, 5, 6,
    7, 8, 9
  };

  Layers {Nx};
  Recombine;
};

// ------------------------------------------------------------
// Extrude the inlet pipe in the negative x-direction
// ------------------------------------------------------------

pipe[] = Extrude {-Lp, 0, 0}
{
  Surface {5};

  Layers {Np};
  Recombine;
};

// ------------------------------------------------------------
// Identify external surfaces
// ------------------------------------------------------------

eps = 1.0e-6 * D;

// Pipe inlet at x = -Lp
inlet[] = Surface In BoundingBox {
  -Lp - eps, y1 - eps, z1 - eps,
  -Lp + eps, y2 + eps, z2 + eps
};

// Lower-y outlet
outlet1[] = Surface In BoundingBox {
  x0 - eps, y0 - eps, z0 - eps,
  x1 + eps, y0 + eps, z3 + eps
};

// Upper-y outlet
outlet2[] = Surface In BoundingBox {
  x0 - eps, y3 - eps, z0 - eps,
  x1 + eps, y3 + eps, z3 + eps
};

// Impingement wall at x = Lx
rightWall[] = Surface In BoundingBox {
  x1 - eps, y0 - eps, z0 - eps,
  x1 + eps, y3 + eps, z3 + eps
};

// Chamber wall at z = z0
lowerWall[] = Surface In BoundingBox {
  x0 - eps, y0 - eps, z0 - eps,
  x1 + eps, y3 + eps, z0 + eps
};

// Chamber wall at z = z3
upperWall[] = Surface In BoundingBox {
  x0 - eps, y0 - eps, z3 - eps,
  x1 + eps, y3 + eps, z3 + eps
};

// Pipe wall at y = y1
pipeYMinus[] = Surface In BoundingBox {
  -Lp - eps, y1 - eps, z1 - eps,
  x0 + eps,  y1 + eps, z2 + eps
};

// Pipe wall at y = y2
pipeYPlus[] = Surface In BoundingBox {
  -Lp - eps, y2 - eps, z1 - eps,
  x0 + eps,  y2 + eps, z2 + eps
};

// Pipe wall at z = z1
pipeZMinus[] = Surface In BoundingBox {
  -Lp - eps, y1 - eps, z1 - eps,
  x0 + eps,  y2 + eps, z1 + eps
};

// Pipe wall at z = z2
pipeZPlus[] = Surface In BoundingBox {
  -Lp - eps, y1 - eps, z2 - eps,
  x0 + eps,  y2 + eps, z2 + eps
};

// x = 0 chamber wall, excluding the pipe opening Surface 5
leftChamberWalls[] = {
  1, 2, 3,
  4,    6,
  7, 8, 9
};

walls[] = {
  leftChamberWalls[],
  rightWall[],
  lowerWall[],
  upperWall[],
  pipeYMinus[],
  pipeYPlus[],
  pipeZMinus[],
  pipeZPlus[]
};

// ------------------------------------------------------------
// Physical boundary groups
// ------------------------------------------------------------

Physical Surface("inlet", 1) = {
  inlet[]
};

Physical Surface("outlet_1", 2) = {
  outlet1[]
};

Physical Surface("outlet_2", 3) = {
  outlet2[]
};

Physical Surface("walls", 4) = {
  walls[]
};

// ------------------------------------------------------------
// Physical volume group
// ------------------------------------------------------------

allVolumes[] = Volume{:};

Physical Volume("domain", 10) = {
  allVolumes[]
};

// ------------------------------------------------------------
// Mesh output settings
// ------------------------------------------------------------

Mesh.ElementOrder = 1;
Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;
Mesh.SaveAll = 0;