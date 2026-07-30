SetFactory("Built-in");

// Total length of the wall
L=1; 
// Number of cells over the wall length
n_cells_L = 4; 
// Size of one cell
d_cell = L/n_cells_L; 
// Number of cells over the wall depth
n_cells_D = 1; 
// Total depth of the wall
D = n_cells_D * d_cell; 
// Number of cells over the wall height
n_cells_H = 1; 
// Total height of the wall
H = n_cells_H * d_cell; 

// First part of the channel
Point(1) = {0, 0, 0, 1.0};
Point(2) = {L/2, 0, 0, 1.0};
Point(3) = {L/2, H, 0, 1.0};
Point(4) = {0, H, 0, 1.0};

Line(1) =  {1,2};
Line(2) =  {2,3};
Line(3) =  {3,4};
Line(4) =  {4,1};

Transfinite Line {1,3} = n_cells_L/2+1 Using Progression 1;
Transfinite Line {2,4} = n_cells_H+1 Using Progression 1;
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1} ;

Transfinite Surface {1};

// Second part of the channel
Point(5) = {L, -d_cell*n_cells_L/2, 0, 1.0};
Point(6) = {L, -d_cell*n_cells_L/2+H, 0, 1.0};

Line(5) =  {2,5};
Line(6) =  {5,6};
Line(7) =  {6,3};
Line(8) =  {3,2};

Transfinite Line {5,7} = n_cells_L/2+1 Using Progression 1;
Transfinite Line {8,6} = n_cells_H+1 Using Progression 1;
Line Loop(2) = {5,6,7,8};
Plane Surface(2) = {2} ;

Transfinite Surface {2};

// Extrusion
out[] = Extrude {0, 0, D} {
  Surface{1}; Layers{n_cells_D}; Recombine;
};

For i In {0:5}
  Printf("out[%g] = %g", i, out[i]);
EndFor


out_2[] = Extrude {0, 0, D} {
  Surface{2}; Layers{n_cells_D}; Recombine;
};

For i In {0:5}
  Printf("out_2[%g] = %g", i, out_2[i]);
EndFor

// Outlet
Physical Surface(0) = {out_2[3]};
// Solid walls
Physical Surface(1) = {out[0], out[1], out[2], out[4], out[5], out_2[0], out_2[1], out_2[2], out_2[4]};

Physical Volume(0) = {1,2};
