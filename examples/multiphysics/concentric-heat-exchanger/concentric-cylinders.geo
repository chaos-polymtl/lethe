// Generates the hexahedral mesh for the flow within the concentric heat exchanger example
// This is achieved by generating 3 concentric circles in 2D, meshing them and extruding the resulting mesh 
// Two physical volumes are used. Volume 0 indicates the fluid, whereas Volume 1 indicates the solid

SetFactory("Built-in");

lc = 0.1;
lf = 0.1;
le = 0.1;

RO=1.;
RI=2.;
RE=3.;
L=50;
nl=50;

Point(0) = {0, 0, 0, lc};
Point(1) = {RO, 0, 0, lc};
Point(2) = {0, -RO , 0, lc};
Point(3) = {-RO, 0, 0, lc};
Point(4) = {0, RO, 0, lc};

Point(5) = {RI, 0, 0, lf};
Point(6) = {0, -RI , 0, lf};
Point(7) = {-RI, 0, 0, lf};
Point(8) = {0, RI, 0, lf};

Point(9) = {RE, 0, 0, le};
Point(10) = {0, -RE , 0, le};
Point(11) = {-RE, 0, 0, le};
Point(12) = {0, RE, 0, le};

Circle(1)={1,0,2};
Circle(2)={2,0,3};
Circle(3)={3,0,4};
Circle(4)={4,0,1};

Circle(5)={5,0,6};
Circle(6)={6,0,7};
Circle(7)={7,0,8};
Circle(8)={8,0,5};

Circle(9)={9,0,10};
Circle(10)={10,0,11};
Circle(11)={11,0,12};
Circle(12)={12,0,9};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8};
Line Loop(3) = {9,10,11,12};

Plane Surface(1) = {1} ;
Plane Surface(2) = {1,2} ;
Plane Surface(3) = {2,3} ;
Recombine Surface{1,2,3};

// Extrusion
Extrude {0, 0, L} {
  Surface{1,2,3}; Layers{nl}; Recombine;
}

Physical Volume(0) = {1, 3};
Physical Volume(1) = {2};
Physical Surface(4) = {2,76,105,109,113,117};
Physical Surface(0) = {1};
Physical Surface(1) = {34};
Physical Surface(2) = {3};
Physical Surface(3) = {118};




