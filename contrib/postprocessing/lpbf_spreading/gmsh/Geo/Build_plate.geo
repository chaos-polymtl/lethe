// Define a variable

y  = 0;
x1 = 0;
x2 = 0.0105;
z1 = 0;
z2 = 0.073;

Point(0) = { x1, y, z1, 1};
Point(1) = { x2, y, z1, 1};
Point(2) = { x2, y, z2, 1};
Point(3) = { x1, y, z2, 1};

Line(0) = {0,1};
Line(1) = {1,2};
Line(2) = {2,0}; 
Line(3) = {0,3};
Line(4) = {3,2};

Line Loop(1) = {0,1,2};
Line Loop(2) = {2,3,4};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Physical Surface(11) = {1,2};




