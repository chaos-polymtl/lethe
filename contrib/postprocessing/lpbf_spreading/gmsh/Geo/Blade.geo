x  = 0.0;
y1 = 0.0;
y2 = 0.1;
z1 = 0.0;
z2 = 0.1;

// Points
Point(0) = { x, y1, z1, 1};
Point(1) = { x, y2, z1, 1};
Point(2) = { x, y2, z2, 1};
Point(3) = { x, y1, z2, 1};

// Line
Line(0) = {0,1};
Line(1) = {1,2};
Line(2) = {2,0}; 
Line(3) = {0,3};
Line(4) = {3,2}; 
Line Loop(1) = {0,1,2};
Line Loop(2) = {2,3,4};


Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Surface(10) = {1,2};




