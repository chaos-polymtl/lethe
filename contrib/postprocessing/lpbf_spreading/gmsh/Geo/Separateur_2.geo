// Define a variable

y1  = 0;
y2  =-0.1;
x1 = 0.0006;
x2 = 0.0036;
x3 = 0.0001;
x4 = 0.;
z1 = 0;
z2 = 0.073;


// Separator
Point(0) = { x1, y1, z1, 1};
Point(1) = { x2, y1, z1, 1};
Point(2) = { x2, y1, z2, 1};
Point(3) = { x1, y1, z2, 1};

Point(4) = { x1, y2, z1, 1};
Point(5) = { x1, y2, z2, 1};
Point(6) = { x2, y2, z1, 1};
Point(7) = { x2, y2, z2, 1};

Line(0)  = {0,1};
Line(1)  = {1,2};
Line(2)  = {2,0}; 
Line(3)  = {0,3};
Line(4)  = {3,2};

Line(5)  = {3,5};
Line(6)  = {5,0};
Line(7)  = {0,4};
Line(8)  = {4,5};

Line(9)  = {2,7};
Line(10) = {7,1};
Line(11) = {1,6};
Line(12) = {6,7};

// Gap
//	Top points
Point(8)  = { x3, y1, z1, 1};
Point(9)  = { x4, y1, z1, 1};
Point(10) = { x4, y1, z2, 1};
Point(11) = { x3, y1, z2, 1};


//	Bottom points
Point(12) = { x3, y2, z1, 1};
Point(13) = { x4, y2, z1, 1};
Point(14) = { x4, y2, z2, 1};
Point(15) = { x3, y2, z2, 1};

// Top of gap
Line(13) = {8,9};
Line(14) = {9,10};
Line(15) = {10,8};
Line(16) = {8,11};
Line(17) = {11,10};

// X_minus gap
Line(18) = {11,15};
Line(19) = {15,8};
Line(20) = {8,12};
Line(21) = {12,15};

// X_plus gap
Line(22) = {10,14};
Line(23) = {14,9};
Line(24) = {9,13};
Line(25) = {13,14};


Line Loop(1)  = {0,1,2};
Line Loop(2)  = {2,3,4};
Line Loop(3)  = {3,5,6};
Line Loop(4)  = {6,7,8};
Line Loop(5)  = {1,9,10};
Line Loop(6)  = {10,11,12};
Line Loop(7)  = {13,14,15};
Line Loop(8)  = {15,16,17};
Line Loop(9)  = {16,18,19};
Line Loop(10) = {19,20,21};
Line Loop(11) = {14,22,23};
Line Loop(12) = {23,24,25};


Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Physical Surface(10) = {1,2,3,4,5,6,7,8,9,10,11,12};




