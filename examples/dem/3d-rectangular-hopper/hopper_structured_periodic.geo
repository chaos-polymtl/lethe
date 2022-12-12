SetFactory("Built-in");

// Geometry variables
right_side = 0.06272;
left_side = -right_side;
hole_right_side = 0.012544;
hole_left_side = -hole_right_side;
hole_height = 0.0;
angle_height = -0.04204;
top = 0.16244;
bottom = -top;
width = 0.0056;

// Division variables
nx = 1;
ny = 7;

// Upper left
Point(1) = {hole_left_side, hole_height, 0, 1.0};
Point(2) = {left_side, hole_height, 0, 1.0};
Point(3) = {left_side, top, 0, 1.0};
Point(4) = {hole_left_side, top, 0, 1.0};

Line(1) =  {1,2};
Line(2) =  {2,3};
Line(3) =  {4,3};
Line(4) =  {1,4};

Transfinite Line {1,3} = Ceil(2*nx+1) Using Progression 1;
Transfinite Line {2,4} = Ceil(ny+1) Using Progression 1;
Line Loop(1) = {1,2,-3,-4};
Plane Surface(1) = {1} ;

Transfinite Surface {1};
Physical Surface(0) = {1};


// Upper middle
Point(5) = {hole_right_side, hole_height, 0, 1.0};
Point(6) = {hole_right_side, top, 0, 1.0};

Line(5) = {5,1};
Line(6) = {6,4};
Line(7) = {5,6};

Transfinite Line {5,6} = Ceil(nx+1) Using Progression 1;
Transfinite Line {4,7} = Ceil(ny+1) Using Progression 1;
Line Loop(2) = {5,4,-6,-7};
Plane Surface(2) = {2} ;

Transfinite Surface {2};
Physical Surface(1) = {2};


// Upper right
Point(7) = {right_side, hole_height, 0, 1.0};
Point(8) = {right_side, top, 0, 1.0};

Line(8) = {7,5};
Line(9) = {8,6};
Line(10) = {7,8};

Transfinite Line {8,9} = Ceil(2*nx+1) Using Progression 1;
Transfinite Line {7,10} = Ceil(ny+1) Using Progression 1;
Line Loop(3) = {8,7,-9,-10};
Plane Surface(3) = {3} ;

Transfinite Surface {3};
Physical Surface(2) = {3};

// Bottom middle
Point(9) = {hole_right_side, bottom, 0, 1.0};
Point(10) = {hole_left_side, bottom, 0, 1.0};

Line(11) = {9,10};
Line(12) = {10,1};
Line(13) = {9,5};

Transfinite Line {11,5} = Ceil(nx+1) Using Progression 1;
Transfinite Line {12,13} = Ceil(ny+1) Using Progression 1;
Line Loop(4) = {11,12,-5,-13};
Plane Surface(4) = {4} ;

Transfinite Surface {4};
Physical Surface(3) = {4};


// Bottom left
Point(11) = {left_side, bottom, 0, 1.0};
Point(12) = {left_side, angle_height, 0, 1.0};

Line(14) = {10,11};
Line(15) = {11,12};
Line(16) = {1,12};

Transfinite Line {14,16} = Ceil(2*nx+1) Using Progression 1;
Transfinite Line {12,15} = Ceil(ny+1) Using Progression 1;
Line Loop(5) = {14,15,-16,-12};
Plane Surface(5) = {5} ;

Transfinite Surface {5};
Physical Surface(4) = {5};


// Bottom right 
Point(13) = {right_side, bottom, 0, 1.0};
Point(14) = {right_side, angle_height, 0, 1.0};

Line(17) = {13,9};
Line(18) = {14,5};
Line(19) = {13,14};

Transfinite Line {13,19} = Ceil(ny+1) Using Progression 1;
Transfinite Line {18,17} = Ceil(2*nx+1) Using Progression 1;
Line Loop(6) = {17,13,-18,-19};
Plane Surface(6) = {6} ;

Transfinite Surface {6};
Physical Surface(5) = {6};


// Reconstruction
Physical Line(1)={1};
Physical Line(2)={2};
Physical Line(3)={3,6,9};
Physical Line(4)={10};
Physical Line(5)={8};
Physical Line(6)={15,16};
Physical Line(7)={18,19};
Physical Line(8)={14,11,17};

Recombine Surface{1,2,3,4,5,6};

// Extrusion
Extrude {0, 0, width} {
  Surface{1,2,3,4,5,6}; Layers{1}; Recombine;
}

Physical Volume(0) = {1, 3, 2, 6, 4, 5};

