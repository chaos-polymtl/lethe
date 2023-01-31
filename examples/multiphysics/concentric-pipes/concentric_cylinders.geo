SetFactory("Built-in");

lc = 0.5;
lf = 0.5;
le = 0.5;

RO=2;
RI=5;
RE=7.5;
L=200;
nl=75;

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


// Reconstruction
//Physical Line(1)={1};
//Physical Line(2)={2};
//Physical Line(3)={3,6,9};
//Physical Line(4)={10};
//Physical Line(5)={8};
//Physical Line(6)={15,16};
//Physical Line(7)={18,19};
//Physical Line(8)={14,11,17};
//
//Recombine Surface{1,2,3,4,5,6};
//

//
Physical Volume(0) = {1, 3};
Physical Volume(1) = {2};
Physical Surface(0) = {105,109,113,117};
Physical Surface(1) = {1};
Physical Surface(2) = {34};
Physical Surface(3) = {3};
Physical Surface(4) = {118};




