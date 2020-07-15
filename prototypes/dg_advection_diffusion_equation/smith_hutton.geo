n=1;
lc=1;
L=1;

Point(0) = {-L,0, 0, lc};
Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0 , 0, lc};
Point(3) = {L, L, 0, lc};
Point(4) = {0, L , 0, lc};
Point(5) = {-L, L, 0, lc};

Line(0)={0,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,0};
Line(6)={1,4};

Line Loop(1) = {0,6,4,5};
Line Loop(2) = {1,2,3,-6};

Plane Surface(1) = {-1} ;
Plane Surface(2) = {-2} ;

Transfinite Line {0:6} = 1;
Transfinite Surface {1,2};
Recombine Surface{1,2};

Physical Surface(0) = {1,2};
Physical Line(0) = {0,2,3,4,5};
Physical Line(1) = {1};


