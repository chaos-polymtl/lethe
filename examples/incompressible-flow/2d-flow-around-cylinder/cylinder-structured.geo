// Define a variable
L=32;
H=16;
y0=0;
y1=H;

x0 =0.;
x1 =L;
//nr=30;
//nt=30;
//nl=150;
nr=60;
nt=60;
nl=150;

xc=8;
yc=8;
r=0.5;
srq=1.414213562/2.*r;

Ls=10*r;
ss=1.414213562/2.*Ls;

xs0=xc-ss;
xs1=xc+ss;
ys0=yc-ss;
ys1=yc+ss;

lc = 0.20;

Point(0) = {x0, y0, 0, lc};
Point(1) = {x0, y1, 0, lc};
Point(2) = {xs0, y0, 0, lc};
Point(3) = {xs0, y1, 0, lc};
Point(4) = {xs1, y0, 0, lc};
Point(5) = {xs1, y1, 0, lc};
Point(6) = {x1, y0, 0, lc};
Point(7) = {x1, y1, 0, lc};

Point(11) = {x0, yc-ss, 0, lc};
Point(12) = {x0, yc+ss, 0, lc};
Point(13) = {x1, yc-ss, 0, lc};
Point(14) = {x1, yc+ss, 0, lc};

Point(20) = {xc, yc, 0, lc};
Point(21) = {xc-srq, yc-srq, 0, lc};
Point(22) = {xc+srq, yc-srq, 0, lc};
Point(23) = {xc+srq, yc+srq, 0, lc};
Point(24) = {xc-srq, yc+srq, 0, lc};

Point(31) = {xc-ss, yc-ss, 0, lc};
Point(32) = {xc+ss, yc-ss, 0, lc};
Point(33) = {xc+ss, yc+ss, 0, lc};
Point(34) = {xc-ss, yc+ss, 0, lc};

Line(60) = {0,11};
Line(61) = {11,12};
Line(63) = {12,1};
Line(64) = {6,13};
Line(65) = {13,14};
Line(66) = {14,7};
Line(1) = {0,2};
Line(2) = {1,3};
Line(3) = {2,4};
Line(4) = {3,5};
Line(5) = {4,6};
Line(6) = {5,7};

Line(10) = {11,12};
Line(11) = {11,31};
Line(12) = {12,34};
Line(13) = {31,34};

Line(20) = {32,33};
Line(21) = {32,13};
Line(22) = {13,14};
Line(23) = {33,14};

Line(30) = {2,31};
Line(31) = {4,32};
Line(32) = {31,32};
Line(33) = {34,3};
Line(34) = {34,33};
Line(35) = {33,5};

Circle(40)={24,20,23};
Circle(41)={23,20,22};
Circle(42)={22,20,21};
Circle(43)={21,20,24};

Line(50)={24,34};
Line(51)={23,33};
Line(52)={22,32};
Line(53)={21,31};

//Theta of the circle
Transfinite Line {3,22,4,10,34,20,32,13,40,41,42,43} = Ceil(nt) Using Progression 1.0;
// Radial direction
Transfinite Line {50,51,52,53} = Ceil(nr) Using Progression 1.05;
Transfinite Line {5,21,23,6} = Ceil(nl) Using Progression 1.00;

//
//
//
Line Loop(41) = {-50,40,51,-34};
Line Loop(42) = {-13,-53,43,50};
Line Loop(43) = {53,32,-52,42};
Line Loop(44) = {52,20,-51,41};


Line Loop (50) = {-60,1,30,-11};
Line Loop (51) = {-10,11,13,-12};
Line Loop (52) = {-63,12,33,-2};
Line Loop (53) = {-30,3,31,-32};
Line Loop (54) = {-33,34,35,-4};
Line Loop (55) = {-31,5,64,-21};
Line Loop (56) = {-20,21,22,-23};
Line Loop (57) = {-35,23,66,-6};

// Creates a physical entity 1 (i.e. for a BC)
//Physical Point(1) = {1,2} ;
//Physical Line(0)={5,6,7,8};
//Physical Line(1)={1};
//Physical Line(2)={2,4};
//Physical Line(3)={3};


Plane Surface(41) = {41};
Plane Surface(42) = {42};
Plane Surface(43) = {43};
Plane Surface(44) = {44};
Plane Surface(50) = {50};
Plane Surface(51) = {51};
Plane Surface(52) = {52};
Plane Surface(53) = {53};
Plane Surface(54) = {54};
Plane Surface(55) = {55};
Plane Surface(56) = {56};
Plane Surface(57) = {57};
Transfinite Surface {41:44,50:57};
Recombine Surface{41:44,50:57};


//Plane Surface(2) = {2};
//Plane Surface(3) = {3};
Physical Surface(1) = {41:44,50:57};
Physical Line(0)={40,41,42,43};
Physical Line(1)={60,10,63};
Physical Line(2)={1,3,5,2,4,6};
Physical Line(3)={64,22,66};
