// Variables
// L_out : step to outlet
// L_in : inlet to step
// h = step height
// beta : expansion ratio = h_out / h_in
L_out = 50.0; //
L_in = 15.0;
h = 1.0;
beta = 2.0;
h_out = h/(1-1/beta);

cell_size = 1.0;

// Points
Point(1) = {L_in+L_out, 0, 0, cell_size};
Point(2) = {L_in+L_out, h_out, 0, cell_size};
Point(3) = {L_in, h, 0, cell_size};
Point(4) = {L_in, 0, 0, cell_size};
Point(5) = {0, h, 0, cell_size};
Point(6) = {0, h_out, 0, cell_size};
Point(7) = {L_in, h_out, 0, 1.0};
Point(8) = {L_in+L_out, h, 0, 1.0};

// Lines
Line(1) = {5,6};
Line(2) = {6,7};
Line(3) = {7,2};
Line(4) = {2,8};
Line(5) = {8,1};
Line(6) = {1,4};
Line(7) = {4,3};
Line(8) = {3,5};
Line(9) = {3,7};

// Surfaces
Curve Loop(1) = {7,9,3,4,5,6};
Plane Surface(1) = {1};
Curve Loop(2) = {8, 1, 2, -9};
Plane Surface(2) = {2};
Recombine Surface {2, 1};

// Quad meshing (transfinite)
Transfinite Surface {1} = {7, 2, 1, 4};
Transfinite Surface {2} = {6, 7, 3, 5};
Transfinite Curve {2, 8} = 26 Using Progression 1; // N-1 inlet cells
Transfinite Curve {3, 6} = 76 Using Progression 1; // N-1 outlet cells
Transfinite Curve {1, 9, 7, 4, 5} = 7 Using Progression 1; // N-1 vertical cells

// Physical properties
// Boundaries
Physical Line(0) = {2,3,6,7,8}; // Wall
Physical Line(1) = {1}; // Inlet
Physical Line(2) = {4,5}; // Outlet

// Surface
Physical Surface(0) = {1,2};

Recombine Surface {2, 1};
