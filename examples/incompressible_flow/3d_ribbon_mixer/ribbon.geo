//SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMin = 0.0001;
Mesh.CharacteristicLengthMax = 0.015;
//Geometry.NumSubEdges = 100; // nicer display of curve

rTank=0.1825;
hTank=0.365;
xRect=0.1304;
yRect=0;
zRect=0.0448;
hRect=0.004;
wRect=0.0361;
zMax=0.333;
rtige=0.0127;
htige=0.333;  //pas le choix pour le maillage detre dans le domaine
l=0.0125; //pour le raffinement
l1=0.0125;


xRect2=xRect*-1;
yRect2=yRect*-1;


nPt=60;
nSpin=1;

//===========================================================================
//Ruban 1
//===========================================================================



// We create series of point for each of the lines of the ribbon
x0=xRect;
z0=zRect;
x1=xRect+wRect;
z1=zRect+hRect;

For i In {0:nPt}
  x = Cos(2*Pi*nSpin*i/nPt)*x0;
  y = Sin(2*Pi*nSpin*i/nPt)*x0;
  z = z0+ i/nPt *(zMax-z0);
  Point(1000+i) = {x,y,z,l};
EndFor

For i In {0:nPt}
  x = Cos(2*Pi*nSpin*i/nPt)*x1;
  y = Sin(2*Pi*nSpin*i/nPt)*x1;
  z = z0+ i/nPt *(zMax-z0);
  Point(2000+i) = {x,y,z,l};
EndFor

For i In {0:nPt}
  x = Cos(2*Pi*nSpin*i/nPt)*x0;
  y = Sin(2*Pi*nSpin*i/nPt)*x0;
  z = z1+ i/nPt *(zMax-z0);
  Point(3000+i) = {x,y,z,l};
EndFor

x0=xRect;
For i In {0:nPt}
  x = Cos(2*Pi*nSpin*i/nPt)*x1;
  y = Sin(2*Pi*nSpin*i/nPt)*x1;
  z = z1+ i/nPt *(zMax-z0);
  Point(4000+i) = {x,y,z,l};
EndFor


// we spline the series of point to create the segment
Spline(1) = {1000:1000+nPt};
Spline(2) = {2000:2000+nPt};
Spline(3) = {3000:3000+nPt};
Spline(4) = {4000:4000+nPt};
//Wire(1) = {1};
//Wire(2) = {2};
//Wire(3) = {3};
//Wire(4) = {4};


// End and final lines
Line(5) = {3000,4000};
Line(6) = {4000,2000};
Line(7) = {2000,1000};
Line(8) = {1000,3000};
Line(9)  = {3000+nPt,4000+nPt};
Line(10) = {4000+nPt,2000+nPt};
Line(11) = {2000+nPt,1000+nPt};
Line(12) = {1000+nPt,3000+nPt};

// Create the surface from the line loops
//Line Loop(100) = {1,-11,-2,7};
Line Loop(100) = {4,-9,-3,5};
Line Loop(101) = {1,12,-3,-8};
Line Loop(102) = {2,-10,-4,6};
Line Loop(103) = {1,-11,-2,7};
Line Loop(200) = {8,5,6,7};
Line Loop(201) = {10,11,12,9};

Surface(1)={100};
Surface(2)={101};
Surface(3)={102};
Surface(4)={103};
Surface(5)={200};
Surface(6)={201};

//===========================================================================
//Ruban 2
//===========================================================================

// We create series of point for each of the lines of the ribbon
x02=xRect2;
z02=zRect;
x2=xRect2-wRect;
z2=zRect+hRect;

For i In {0:nPt}
  x = Cos(2*Pi*nSpin*i/nPt)*x02;
  y = Sin(2*Pi*nSpin*i/nPt)*x02;
  z = z02+ i/nPt *(zMax-z02);
  Point(10000+i) = {x,y,z,l};
EndFor

For i In {0:nPt}
  x = Cos(2*Pi*nSpin*i/nPt)*x2;
  y = Sin(2*Pi*nSpin*i/nPt)*x2;
  z = z02+ i/nPt *(zMax-z02);
  Point(20000+i) = {x,y,z,l};
EndFor

For i In {0:nPt}
  x = Cos(2*Pi*nSpin*i/nPt)*x02;
  y = Sin(2*Pi*nSpin*i/nPt)*x02;
  z = z2+ i/nPt *(zMax-z02);
  Point(30000+i) = {x,y,z,l};
EndFor

x0=xRect;
For i In {0:nPt}
  x = Cos(2*Pi*nSpin*i/nPt)*x2;
  y = Sin(2*Pi*nSpin*i/nPt)*x2;
  z = z2+ i/nPt *(zMax-z02);
  Point(40000+i) = {x,y,z,l};
EndFor


// we spline the series of point to create the segment
Spline(21) = {10000:10000+nPt};
Spline(22) = {20000:20000+nPt};
Spline(23) = {30000:30000+nPt};
Spline(24) = {40000:40000+nPt};



// End and final lines
Line(25) = {30000,40000};
Line(26) = {40000,20000};
Line(27) = {20000,10000};
Line(28) = {10000,30000};
Line(29)  = {30000+nPt,40000+nPt};
Line(30) = {40000+nPt,20000+nPt};
Line(31) = {20000+nPt,10000+nPt};
Line(32) = {10000+nPt,30000+nPt};

// Create the surface from the line loops
Line Loop(1000) = {24,-29,-23,25};
Line Loop(1010) = {21,32,-23,-28};
Line Loop(1020) = {22,-30,-24,26};
Line Loop(1030) = {21,-31,-22,27};
Line Loop(2000) = {28,25,26,27};
Line Loop(2010) = {30,31,32,29};

Surface(10)={1000};
Surface(20)={1010};
Surface(30)={1020};
Surface(40)={1030};
Surface(50)={2000};
Surface(60)={2010};

//===========================================================================
//Tige interieur
//===========================================================================


Point (30) = { 0, 0 ,zRect,l};
Point (31) = {rtige, 0 , zRect,l};
Point (32) = {0, rtige , zRect,l};
Point (33) = {-rtige, 0 , zRect,l};
Point (34) = {0, -rtige , zRect,l};
Point (310) = { 0, 0 ,htige,l};
Point (311) = {rtige, 0 ,  htige,l};
Point (312) = {0, rtige ,  htige,l};
Point (313) = {-rtige, 0 , htige,l};
Point (314) = {0, -rtige , htige,l};

Circle(31000) = {31,30,32};
Circle(31001) = {32,30,33};
Circle(31002) = {33,30,34};
Circle(31003) = {34,30,31};
Circle(31004) = {311,310,312};
Circle(31005) = {312,310,313};
Circle(31006) = {313,310,314};
Circle(31007) = {314,310,311};
Line(31011) = {31,311};
Line(31012) = {32,312};
Line(31013) = {33,313};
Line(31014) = {34,314};

Line Loop(31) ={31000:31003};
Line Loop(32) ={31004:31007};

Line Loop(33) = {31003,31011,-31007,-31014};
Line Loop(34) = {31000,31012,-31004,-31011};
Line Loop(35) = {31001,31013,-31005,-31012};
Line Loop(36) = {31002,31014,-31006,-31013};

Plane Surface (3101) = {31};
Plane Surface (3102) = {32};
Surface (3103) = {33};
Surface (3104) = {34};
Surface (3105) = {35};
Surface (3106) = {36};


//===========================================================================
//Cylindre 
//===========================================================================


Point (0) = { 0, 0 ,0,l1};
Point (1) = {rTank, 0 , 0,l1};
Point (2) = {0, rTank , 0,l1};
Point (3) = {-rTank, 0 , 0,l1};
Point (4) = {0, -rTank , 0,l1};
Point (10) = { 0, 0 ,hTank,l1};
Point (11) = {rTank, 0 ,  hTank,l1};
Point (12) = {0, rTank ,  hTank,l1};
Point (13) = {-rTank, 0 , hTank,l1};
Point (14) = {0, -rTank , hTank,l1};

Circle(1000) = {1,0,2};
Circle(1001) = {2,0,3};
Circle(1002) = {3,0,4};
Circle(1003) = {4,0,1};
Circle(1004) = {11,10,12};
Circle(1005) = {12,10,13};
Circle(1006) = {13,10,14};
Circle(1007) = {14,10,11};
Line(1011) = {1,11};
Line(1012) = {2,12};
Line(1013) = {3,13};
Line(1014) = {4,14};

Line Loop(1) ={1000:1003};
Line Loop(2) ={1004:1007};

Line Loop(3) = {1003,1011,-1007,-1014};
Line Loop(4) = {1000,1012,-1004,-1011};
Line Loop(5) = {1001,1013,-1005,-1012};
Line Loop(6) = {1002,1014,-1006,-1013};

Plane Surface (101) = {1};
Plane Surface (102) = {2};
Surface (103) = {3};
Surface (104) = {4};
Surface (105) = {5};
Surface (106) = {6};

//===========================================================================
//Volume final
//===========================================================================

Surface Loop (2) = {101:106,-1,-2,-3,-4,-5,-6,-10,-20,-30,-40,-50,-60,-3101,-3102,-3103,-3104,-3105,-3106};
Volume (2) = {2};

Physical Surface (0) = {101,103:106}; // Cote
Physical Surface (1) = {102}; // Haut
//Physical Surface ("Ruban1") = {1:6};
//Physical Surface ("Ruban2") = {10:60};
Physical Surface (2) = {1:6,10:60,3101:3106}; // Blade
//Physical Surface ("Tige") = {3101:3106};
Physical Volume (0) = {2};
