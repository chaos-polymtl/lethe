// Define a variable

lc = 1;
H=0.5;
W=0.4;
x0=0.6;


Point(0) = {x0, 0, 0, lc};
Point(1) = {x0, 0, W, lc};
Point(2) = {x0, H, W, lc};
Point(3) = {x0, H, 0, lc};

Line(0)={0,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,0};

Line Loop(1) = {0,1,2,3};

Plane Surface(1) = {1} ;

Physical Surface(0) = {1};

