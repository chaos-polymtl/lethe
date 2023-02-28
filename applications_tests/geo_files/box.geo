// Box geometry with two physical volumes
// This is a primitive geometry which is used for the application test
// that use the solid region mechanism (e.g. conjugate heat transfer)

H=0.1;
H2=0.2;
L=1;
Mesh.CharacteristicLengthMin = 0.001;
Mesh.CharacteristicLengthMax = 0.015;


Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, H, 0};
Point(4) = {0, H, 0};
Point(5) = {L, H2, 0};
Point(6) = {0, H2, 0};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={4,3};
Line(4)={1,4};
Line(5)={3,5};
Line(6)={5,6};
Line(7)={6,4};
Line Loop(1) = {1,2,-3,-4};
Plane Surface(1) = {1} ;

Line Loop(2) = {3,5,6,7};
Plane Surface(2) = {2} ;

Physical Surface(0) = {1};
Physical Surface(1) = {2};
Physical Line(0)={4};
Physical Line(1)={2};
Physical Line(2)={1};
Physical Line(3)={6};
Physical Line(4)={5,7};


