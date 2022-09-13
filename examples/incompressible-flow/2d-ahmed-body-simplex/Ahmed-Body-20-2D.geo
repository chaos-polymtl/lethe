//===========================================================================
//Parameters
//===========================================================================
unit = 1000; //length unit : 1 -> mm ; 1000 -> m
phi = 20; //angle at rear, variable
esf = 2.0e-1; //element size factor, used in the free quad mesh

Cphi = Cos(phi*Pi/180);
Sphi = Sin(phi*Pi/180);
Tphi = Tan(phi*Pi/180);

//Ahmed body basic geometry
L = 1044/unit;
H = 288/unit;
R = 100/unit;
Hw = 50/unit; //wheel height (height from the road)
Ls = 222/unit; //slope length

//Fluid domain
xmin = -500/unit;
ymin = -Hw;
xmax = 2500/unit;
ymax = 1000/unit;

//===========================================================================
//2D-Ahmed body with a 20 degree-step
//===========================================================================
//Origin point at bottom left
//Straight lines
Point(1) = {R, 0, 0, esf};
Point(2) = {L, 0, 0, esf};
Point(3) = {L, H-Ls*Sphi, 0, esf};
Point(4) = {L-Ls*Cphi, H, 0, esf};
Point(5) = {R, H, 0, esf};
Point(6) = {0, H-R, 0, esf};
Point(7) = {0, R, 0, esf};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {6, 7};

//Circle arcs at the front
//Center points
Point(8) = {R, H-R, 0, esf};
Point(9) = {R, R, 0, esf};
//Points added at mid-arc for the mesh
midArc = R*Cos(45*Pi/180); //xmid=ymid
Point(10) = {R-midArc, R-midArc, 0, esf};
Point(11) = {R-midArc, H-R+midArc, 0, esf};
//Circle arcs
Circle(6) = {5, 8, 11};
Circle(7) = {11, 8, 6};
Circle(8) = {7, 9, 10};
Circle(9) = {10, 9, 1};

//Points for the corners of the fluid domain
Point(20) = {xmin, ymin, 0, esf};
Point(21) = {xmin, ymax, 0, esf};
Point(22) = {xmax, ymax, 0, esf};
Point(23) = {xmax, ymin, 0, esf};

//===========================================================================
//Layers for the Mesh refinement
//===========================================================================
//Bondary layer
//---------------------------------------------------------------------------
Wbl = Hw/2; //bondary layer width
T45 = Tan(45*Pi/180);
C45 = Cos(45*Pi/180); //=Sin(45*Pi/180)

//Points
Point(30) = {R, -Wbl, 0, esf};
Point(31) = {L+Wbl*T45, -Wbl, 0, esf};
Point(32) = {L+Wbl*T45, H+Wbl-Tphi*(Wbl*T45+Ls*Cphi), 0, esf};
Point(33) = {L-Ls*Cphi, H+Wbl, 0, esf};
Point(34) = {R, H+Wbl, 0, esf};
Point(35) = {R-midArc-Wbl*C45, H-R+midArc+Wbl*C45, 0, esf};
Point(36) = {-Wbl, H-R, 0, esf};
Point(37) = {-Wbl, R, 0, esf};
Point(38) = {R-midArc-Wbl*C45, R-midArc-Wbl*C45, 0, esf};

//Contour lines
Line(14) = {30, 31};
Line(15) = {31, 32};
Line(16) = {32, 33};
Line(17) = {33, 34};
Line(18) = {36, 37};
Circle(28) = {34, 8, 35};
Circle(29) = {35, 8, 36};
Circle(30) = {37, 9, 38};
Circle(31) = {38, 9, 30};

//Inner lines
Line(19) = {36, 6};
Line(20) = {37, 7};
Line(21) = {11, 35};
Line(22) = {5, 34};
Line(23) = {4, 33};
Line(24) = {3, 32};
Line(25) = {2, 31};
Line(26) = {1, 30};
Line(27) = {10, 38};

//Bottom layer
//---------------------------------------------------------------------------
Point(40) = {R-midArc-Wbl*C45, ymin, 0, esf};
Point(41) = {R, ymin, 0, esf};
Point(42) = {L+Wbl*T45, ymin, 0, esf};

//Outer layer
//---------------------------------------------------------------------------
Line(50) = {20, 21};
Line(51) = {21, 22};
Line(52) = {22, 23};
Line(53) = {38, 40};
Line(54) = {40, 41};
Line(55) = {41, 30};
Line(56) = {41, 42};
Line(57) = {42, 31};
Line(58) = {20, 40};
Line(59) = {42, 23};

//===========================================================================
//Transfinite mesh definition
//===========================================================================
//Bondary layer
Transfinite Line {5, 18, 6, 28, 7, 29, 8, 9, 30, 31, 54} = 4 Using Progression 1;
Transfinite Line {4, 17} = 15 Using Progression 1;
Transfinite Line {1, 14, 56} = 16 Using Progression 1;
Transfinite Line {3, 16} = 5 Using Progression 1;
Transfinite Line {2, 15} = 3 Using Progression 1;
Transfinite Line {23, 24, 25, 26, 27, 20, 19, 21, 22} = 3 Using Progression 1;

//Road layer
Transfinite Line {53, 55, 57} = 2 Using Progression 1;

//===========================================================================
//Surfaces
//===========================================================================
//Beware that all Line Loops turn counter-clockwise (no black cells should remain on Gmsh)
//Bondary layer
Line Loop(1) = {-6, -21, 28, 22}; Plane Surface(1) = {1};
Line Loop(2) = {-7, 19, 29, 21}; Plane Surface(2) = {2};
Line Loop(3) = {-19, -5, 20, 18}; Plane Surface(3) = {3};
Line Loop(4) = {-20, -8, -27, 30}; Plane Surface(4) = {4};
Line Loop(5) = {27, 31, -26, -9}; Plane Surface(5) = {5};
Line Loop(6) = {-1, -25, 14, 26}; Plane Surface(6) = {6};
Line Loop(7) = {-2, 25, 15, -24}; Plane Surface(7) = {7};
Line Loop(8) = {-3, 24, 16, -23}; Plane Surface(8) = {8};
Line Loop(9) = {-4, 23, 17, -22}; Plane Surface(9) = {9};
Transfinite Surface {1:9};
//Recombine Surface {1:9};

//Road layer
Line Loop(10) = {55, -31, 53, 54}; Plane Surface(10) = {10};
Line Loop(11) = {14, -57, -56, 55}; Plane Surface(11) = {-11};
Transfinite Surface {10,11};
//Recombine Surface {10,11};

//Outer layer
Line Loop(12) = {52, -59, 57, 15, 16, 17, 28, 29, 18, 30, 53, -58, 50, 51};
Plane Surface(12) = {-12};
//Recombine Surface {12};

//===========================================================================
//Bondary lines
//===========================================================================
//Ahmed Body and road : bc = 0, noslip
Physical Line(0) = {1:9, 58,54,56,59};
//Flow inlet : bc = 1, imposed velocity
Physical Line(1) = {50};
//Top : bc = 2, slip condition
Physical Line(2) = {51};
//Flow outlet : implied bondary condition, free region, pressure becomes close to 0
Physical Line(3) = {52};
//All mesh blocks should be in a physical group
Physical Surface(0) = {1:9, 10:12};

//===========================================================================
