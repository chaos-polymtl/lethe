// BLADES GEOMETRY
// Valérie Bibeau, Polytechnique Montréal
// 2020

SetFactory("OpenCASCADE");

// -------------------------------------------
// Dimensionless geometry variables
// -------------------------------------------

T = 1.000;

//ratioTD = 3;
//  D = T/ratioTD;
D = 0.333;
ratioHT = 1;
  H = T*ratioHT;
ratioTC = 4;
  C = T/ratioTC;
//ratioDW = 5;
//  W = D/ratioDW;
//ratioDW = 10;
//  W_Hub = D/ratioDW;
W = 0.05;
W_Hub = 0.10;
H_blade = D/4;
E = 0.25*W;



theta = Pi/4;
d = W*Cos(theta);
r = d/2;
h = W*Sin(theta);


// -------------------------------------------
// Mesh
// -------------------------------------------
Mesh.CharacteristicLengthMin = 0.001;
Mesh.CharacteristicLengthMax = 0.02;
//Mesh.ElementOrder = 1;
//Mesh.SecondOrderLinear = 0;
//Mesh.HighOrderOptimize = 1;
//Mesh.SubdivisionAlgorithm = 2; // quad

// -------------------------------------------
// Cylinder (pbt)
// -------------------------------------------
Cylinder(2) = {C-H/2, 0, 0, H-C, 0, 0, W/2, 2*Pi};

// -------------------------------------------
// Cylinder (pbt hub)
// -------------------------------------------
Cylinder(3) = {C-H/2, 0, 0, H_blade, 0, 0, W_Hub/2, 2*Pi};

// -------------------------------------------
// Blade 1
// -------------------------------------------
Box (4) = { C-H/2, -E/2, 0, H_blade, E, D/2 };
Rotate {{ 0,0,1 }, {C-H/2+H_blade/2, 0, 0}, theta } {Volume{4};}

// Other blades
Rotate {{ 1,0,0 }, {C-H/2+H_blade/2, 0, 0}, Pi/2 } {Duplicata{Volume{4};}}
Rotate {{ 1,0,0 }, {C-H/2+H_blade/2, 0, 0}, 2*Pi/2 } {Duplicata{Volume{4};}}
Rotate {{ 1,0,0 }, {C-H/2+H_blade/2, 0, 0}, 3*Pi/2 } {Duplicata{Volume{4};}}

BooleanUnion {Volume{2}; Delete;}{Volume{3:7}; Delete;}

Physical Surface("0") = {0:22};
Physical Volume(0) = {1};

