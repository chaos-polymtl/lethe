// Gmsh project created on Tue Sep 13 15:49:36 2022
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {-0.02, 0, 0, 0.04, 0, 0, 0.01525, 2*Pi};
//+
Physical Volume("fluid", 4) = {1};
//+
Physical Surface("inlet", 5) = {3};
//+
Physical Surface("outlet", 6) = {2};
//+
Physical Surface("wall", 7) = {1};
//+
Extrude {0.04, 0, 0} {
  Surface{3}; Layers {5}; Recombine;
}
